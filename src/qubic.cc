#include "qubic.h"
#include "charset.h"
#include "edge_list.h"
#include "version.h"
#include "discretize.h"
#include <cmath> // floor, ceil
#include <cstring> // strcmp
#include <cstddef> // size_t
#include <algorithm>
#include <list>

namespace internal {
  static void seed_update(const DiscreteArray& s, std::vector<std::vector<bits16>>& profile) {
    for (std::size_t i = 0; i < s.size(); i++)
      profile[i][s[i]]++;
  }

  /* scan through all columns and identify the set within threshold,
     * "fuzziness" of the block is controlled by TOLERANCE (-c)
     */
  static void scan_block(const DiscreteArrayList& arr_c, const Symbols& symbols, const std::vector<int>& gene_set,
    Block& b, std::vector<std::vector<bits16>>& profile, double TOLERANCE) {
    std::size_t block_rows, cur_rows;
    block_rows = cur_rows = gene_set.size();
    for (std::size_t j = 0; j < profile.size(); j++)
      for (std::size_t k = 0; k < profile[j].size(); k++)
        profile[j][k] = 0;
    for (std::size_t j = 0; j < cur_rows; j++)
      seed_update(arr_c[gene_set[j]], profile);
    int btolerance = static_cast<int>(std::ceil(TOLERANCE * block_rows));
    for (std::size_t j = 0; j < profile.size(); j++) {
      /* See if this column satisfies tolerance */
      /* here i start from 1 because symbols[0]=0 */
      for (std::size_t i = 1; i < symbols.size(); i++) {
        if (profile[j][i] >= btolerance) {
          b.conds.insert(j);
          break;
        }
      }
    }
  }

  /*************************************************************************/

  static void update_colcand(std::list<std::size_t>& colcand, const DiscreteArray& g1, const DiscreteArray& g2) {
    std::list<std::size_t>::iterator it = colcand.begin();
    while (it != colcand.end()) {
      if (g1[*it] != g2[*it]) colcand.erase(it++); // alternatively, it = colcand.erase(it);
      else ++it;
    }
  }

  /*calculate the weight of the edge with two vertices g1 and g2*/
  static int intersect_row(const std::list<std::size_t>& colcand, const DiscreteArray& g1, const DiscreteArray& g2) {
    int cnt = 0;
    for (auto it = colcand.begin(); it != colcand.end(); ++it) if (g1[*it] == g2[*it] && g1[*it] != 0) cnt++;
    return cnt;
  }

  /*calculate the negative correlation between g1 and g2*/
  inline int reverse_row(const std::list<std::size_t>& colcand, const DiscreteArray& g1, const DiscreteArray& g2,
    const Symbols& symbols) {
    int cnt = 0;
    for (auto it = colcand.begin(); it != colcand.end(); ++it) if (symbols[g1[*it]] == -symbols[g2[*it]]) cnt++;
    return cnt;
  }

  /* calculate the coverage of any row to the current consensus
    * cnt = # of valid consensus columns
    */
  static int seed_current_modify(const DiscreteArray& s, std::list<std::size_t>& colcand, const int components,
    std::vector<std::vector<bits16>>& profile, double TOLERANCE) {
    std::size_t n;
    int threshold = static_cast<int>(std::ceil(components * TOLERANCE));
    int cnt = 0;
    for (std::size_t i = 0; i < profile.size(); i++) {
      std::size_t flag = 0;
      short ss = s[i];
      for (std::size_t k = 1; k < profile[i].size(); k++) {
        n = profile[i][k];
        if (static_cast<short>(k) == ss) n++;
        if (static_cast<int>(n) >= threshold) {
          flag = k;
          break;
        }
      }
      if (flag) {
        cnt++;
        colcand.push_back(i);
      }
    }
    return cnt;
  }

  /*check whether current edge can be treat as a seed*/
  static bool check_seed(const Edge* e, const std::vector<Block>& bb, std::size_t rows) {
    int block_id = bb.size();
    int i, b1, b2, b3;
    b1 = b2 = -1;
    for (i = 0; i < block_id; i++)
      if (bb[i].contains(e->gene_one) && bb[i].contains(e->gene_two))
        return false;
    std::vector<int> profiles(rows, 0);
    bool fg = false;
    for (i = 0; i < block_id; i++)
      if (bb[i].contains(e->gene_one)) {
        fg = true;
        break;
      }
    if (fg)
      b1 = i;
    fg = false;
    for (i = 0; i < block_id; i++)
      if (bb[i].contains(e->gene_two)) {
        fg = true;
        break;
      }
    if (fg)
      b2 = i;
    if (b1 == -1 || b2 == -1)
      return true;
    for (std::set<int>::iterator it = bb[b1].genes_order.begin(); it != bb[b1].genes_order.end(); ++it)
      profiles[*it]++;
    for (std::set<int>::iterator it = bb[b1].genes_reverse.begin(); it != bb[b1].genes_reverse.end(); ++it)
      profiles[*it]++;
    for (std::set<int>::iterator it = bb[b2].genes_order.begin(); it != bb[b2].genes_order.end(); ++it)
      profiles[*it]++;
    for (std::set<int>::iterator it = bb[b2].genes_reverse.begin(); it != bb[b2].genes_reverse.end(); ++it)
      profiles[*it]++;
    for (std::size_t index = 0; index < rows; index++)
      if (profiles[index] > 1) return false;
    b3 = std::max(bb[b1].block_cols(), bb[b2].block_cols());
    return !(e->score < b3/* (bb[b1]->block_cols + bb[b2]->block_cols) / 2*/);
  }

  static long double get_pvalue(const continuous& a, const int& b) {
    long double pvalue = 0;
    long double poisson = 1.0 / exp(a);
    for (int i = 0; i < b + 300; i++) {
      if (i > b - 1) pvalue += poisson;
      else poisson *= a / (i + 1.0);
    }
    return pvalue;
  }

  static void block_init(const DiscreteArrayList& arr_c, Block& b,
    std::vector<int>& genes, std::vector<int>& scores,
    std::vector<bool>& candidates, const int& cand_threshold,
    std::size_t& components, std::vector<long double>& pvalues, bool IS_cond, std::size_t COL_WIDTH, bool IS_area) {
    std::size_t rows = arr_c.size();
    std::size_t cols = arr_c[0].size();
    int score, top;
    int cnt = 0, cnt_all = 0, pid = 0;
    continuous cnt_ave = 0, row_all = static_cast<continuous>(rows);
    long double pvalue;
    int max_cnt, max_i;
    std::vector<int> arr_rows(rows), arr_rows_b(rows);
    std::list<std::size_t> colcand;
    DiscreteArray g1, g2;
    g1 = arr_c[genes[0]];
    g2 = arr_c[genes[1]];
    for (std::size_t i = 0; i < cols; i++)
      if (g1[i] == g2[i] && g1[i] != 0)
        colcand.push_back(i);
    for (std::size_t i = 0; i < rows; i++) {
      arr_rows[i] = intersect_row(colcand, arr_c[genes[0]], arr_c[i]);
      arr_rows_b[i] = arr_rows[i];
    }
    /*we just get the largest 100 rows when we initial a bicluster because we believe that
       * the 100 rows can characterize the structure of the bicluster
       * btw, it can reduce the time complexity*/
    if (rows > 100) {
      std::sort(arr_rows_b.begin(), arr_rows_b.end());
      top = arr_rows_b[rows - 100];
      for (std::size_t i = 0; i < rows; i++)
        if (arr_rows[i] < top)
          candidates[i] = false;
    }
    /*calculate the condition low bound for current seed*/
    int cutoff = static_cast<int>(0.05 * rows);
    b.cond_low_bound = arr_rows_b[rows - cutoff - 1];
    while (components < rows) {
      max_cnt = -1;
      max_i = -1;
      components++;
      cnt_all = 0;
      cnt_ave = 0;
      /******************************************************/
      /*add a function of controlling the bicluster by pvalue*/
      /******************************************************/
      for (std::size_t i = 0; i < rows; i++) {
        if (!candidates[i]) continue;
        cnt = intersect_row(colcand, arr_c[genes[0]], arr_c[i]);
        cnt_all += cnt;
        if (cnt < cand_threshold)
          candidates[i] = false;
        if (cnt > max_cnt) {
          max_cnt = cnt;
          max_i = i;
        }
      }
      cnt_ave = cnt_all / row_all;
      pvalue = get_pvalue(cnt_ave, max_cnt);
      if (IS_cond) {
        if (max_cnt < static_cast<int>(COL_WIDTH) || max_i < 0 || max_cnt < b.cond_low_bound) break;
      } else {
        if (max_cnt < static_cast<int>(COL_WIDTH) || max_i < 0) break;
      }
      if (IS_area) score = components * max_cnt;
      else score = std::min(static_cast<int>(components), max_cnt);
      if (score > b.score) b.score = score;
      if (pvalue < b.pvalue) b.pvalue = pvalue;
      genes.push_back(max_i);
      scores.push_back(score);
      pvalues[pid++] = pvalue;
      update_colcand(colcand, arr_c[genes[0]], arr_c[max_i]);
      candidates[max_i] = false;
    }
    /*be sure to free a pointer when you finish using it*/
  }

  /* compare function for qsort, descending by score */
  static bool block_cmpr(const Block& a, const Block& b) {
    return a.score > b.score;
  }

  /************************************************************************/
  static std::vector<Block> report_blocks(std::vector<Block> bb, int RPT_BLOCK, double FILTER) {
    std::vector<Block> output;
    int num = bb.size();
    std::sort(bb.begin(), bb.end(), block_cmpr);
    /*MIN MAX et al functions can be accessed in struct.h*/
    int n = std::min(num, RPT_BLOCK);
    /*double proportion;*/
    /* the major post-processing here, filter overlapping blocks*/
    int i = 0;
    int j = 0;
    while (i < num && j < n) {
      Block& b_ptr = bb[i];
      double cur_rows = b_ptr.block_rows();
      double cur_cols = b_ptr.block_cols();
      bool flag = true;
      int k = 0;
      while (k < j) {
        double inter_rows = count_intersect(output[k].genes_order, b_ptr.genes_order) +
          count_intersect(output[k].genes_order, b_ptr.genes_reverse) +
          count_intersect(output[k].genes_reverse, b_ptr.genes_order) +
          count_intersect(output[k].genes_reverse, b_ptr.genes_reverse);
        double inter_cols = count_intersect(output[k].conds, b_ptr.conds);
        if (inter_rows * inter_cols > FILTER * cur_rows * cur_cols) {
          flag = false;
          break;
        }
        k++;
      }
      i++;
      if (flag) {
        output.push_back(b_ptr);
        j++;
      }
    }
    return output;
  }

  /************************************************************************/
  static std::vector<Block> cluster(const DiscreteArrayListWithSymbols& all, const std::vector<Edge *>& el,
    std::size_t COL_WIDTH, double TOLERANCE, bool IS_cond, bool IS_area,
    bool IS_pvalue, std::size_t SCH_BLOCK, int RPT_BLOCK, double FILTER, double f, bool verbose) {
    std::vector<Block> bb;
    std::size_t rows = all.list.size();
    std::size_t cols = all.list[0].size();
    std::vector<long double> pvalues(rows);
    std::vector<bool> candidates(rows);
    std::set<int> allincluster;
    for (std::vector<Edge *>::const_iterator it = el.begin(); it != el.end(); ++it) {
      const Edge* e = *it;
      /* check if both genes already enumerated in previous blocks */
      bool flag = true;
      /* speed up the program if the rows bigger than 200 */
      if (rows > 250) {
        if (allincluster.find(e->gene_one) != allincluster.end() && allincluster.find(e->gene_two) != allincluster.end()) flag = false;
      } else flag = check_seed(e, bb, rows);
      if (!flag) continue;
      std::vector<std::vector<bits16>> profile(cols, std::vector<bits16>(all.symbols.size(), 0));
      /*you must allocate a struct if you want to use the pointers related to it*/
      Block b;
      /*initial the b->score*/
      b.score = std::min(2, e->score);
      /*initial the b->pvalue*/
      b.pvalue = 1;
      /* initialize the stacks genes and scores */
      std::vector<int> genes_order, genes_reverse, scores;
      genes_order.reserve(rows);
      genes_reverse.reserve(rows);
      scores.reserve(rows);
      genes_order.push_back(e->gene_one);
      genes_order.push_back(e->gene_two);
      scores.push_back(1);
      scores.push_back(b.score);
      /* branch-and-cut condition for seed expansion */
      int cand_threshold = static_cast<int>(std::floor(COL_WIDTH * TOLERANCE));
      if (cand_threshold < 2) cand_threshold = 2;
      /* maintain a candidate list to avoid looping through all rows */
      for (std::size_t j = 0; j < rows; j++) candidates[j] = true;
      candidates[e->gene_one] = candidates[e->gene_two] = false;
      std::size_t components = 2;
      /* expansion step, generate a bicluster without noise */
      block_init(all.list, b, genes_order, scores, candidates, cand_threshold, components, pvalues, IS_cond, COL_WIDTH,
        IS_area);
      /* track back to find the genes by which we get the best score*/
      std::size_t k;
      for (k = 0; k < components; k++) {
        if (IS_pvalue && (pvalues[k] == b.pvalue && k >= 2 && scores[k] != scores[k + 1])) break;
        if (scores[k] == b.score && (k + 1 == scores.size() || scores[k + 1] != b.score)) break;
      }
      components = k + 1;
      std::fill(candidates.begin(), candidates.end(), true);
      for (std::size_t ki = 0; ki < components - 1; ki++) {
        seed_update(all.list[genes_order[ki]], profile);
        candidates[genes_order[ki]] = false;
      }
      candidates[genes_order[k]] = false;
      genes_order.resize(k + 1);
      std::list<std::size_t> colcand;
      /* add columns satisfy the conservative r */
      int cnt = seed_current_modify(all.list[genes_order[k]], colcand, components, profile, TOLERANCE);
      /* add some new possible genes */
      int m_cnt;
      for (std::size_t ki = 0; ki < rows; ki++) {
        m_cnt = intersect_row(colcand, all.list[genes_order[0]], all.list[ki]);
        if (candidates[ki] && m_cnt >= std::floor(cnt * TOLERANCE)) {
          genes_order.push_back(ki);
          components++;
          candidates[ki] = false;
        }
      }
      /* add genes that negative regulated to the consensus */
      for (std::size_t ki = 0; ki < rows; ki++) {
        m_cnt = reverse_row(colcand, all.list[genes_order[0]], all.list[ki], all.symbols);
        if (candidates[ki] && m_cnt >= std::floor(cnt * TOLERANCE)) {
          genes_reverse.push_back(ki);
          components++;
          candidates[ki] = false;
        }
      }
      /* store gene arrays inside block */
      scan_block(all.list, all.symbols, genes_order, b, profile, TOLERANCE);
      if (b.block_cols() == 0) continue;
      if (IS_pvalue) b.score = static_cast<int>(-(100 * log(b.pvalue)));
      else b.score = components * b.block_cols();
      for (std::vector<int>::iterator iterator = genes_order.begin(); iterator != genes_order.end(); ++iterator) {
        b.genes_order.insert(*iterator);
        allincluster.insert(*iterator);
      }
      for (std::vector<int>::iterator iterator = genes_reverse.begin(); iterator != genes_reverse.end(); ++iterator) {
        b.genes_reverse.insert(*iterator);
        allincluster.insert(*iterator);
      }

      if (f && (b.block_cols() <= 1 || b.block_rows() <= 1)) continue;

      /*save the current block b to the block list bb so that we can sort the blocks by their score*/
      bb.push_back(b);
      /* reaching the results number limit */
      if (bb.size() == SCH_BLOCK) break;
      if (verbose) fputc('.', stdout);
    }
    if (verbose) fprintf(stdout, "\n");
    return report_blocks(bb, RPT_BLOCK, FILTER);
  }
}

class qubic {
  static std::vector<Block> init_qubic(DiscreteArrayListWithSymbols& all, const double c, const double f, std::size_t col_width,
    const int o, const Option& option, const CountHelper& count_helper, const bool verbose) {
    if (verbose) fprintf(stdout, "\nQUBIC %s: greedy biclustering\n\n", VER);
    /* ensure enough searching space */
    int SCH_BLOCK = 2 * o;
    /* the file that stores all blocks */
    EdgeList EdgeList(col_width, count_helper, verbose);
    /* bi-clustering */
    if (verbose) fprintf(stdout, "Clustering started");
    return internal::cluster(all, EdgeList.get_edge_list(), col_width, c, option.cond_, option.area_, option.pvalue_, SCH_BLOCK, o, f, option.filter_1xn_nx1, verbose);
  }

public:
  static std::vector<Block> init_qubic(DiscreteArrayListWithSymbols& all, const double c, const double f, std::size_t col_width,
    const int o, const Option& option, const bool verbose) {
    return init_qubic(all, c, f, col_width, o, option, CountHelper(all.list), verbose);
  }

  static std::vector<Block> init_qubic(DiscreteArrayListWithSymbols& all, const double c, const double f, std::size_t col_width,
    const int o, const Option& option, const bool verbose, const std::vector<std::vector<float>>& weights) {
    return init_qubic(all, c, f, col_width, o, option, WeightedCountHelper(all.list, weights), verbose);
  }
};

/* Open a file to write or die */
FILE* mustOpenWrite(const char* fileName) {
  if (!strcmp(fileName, "stdin"))  return stdin;
  if (!strcmp(fileName, "stdout")) return stdout;
  FILE* f = fopen(fileName, "w");
  if (f) return f;
  fprintf(stderr, "[Error] Can't open %s to write.", fileName);
  throw - 1;
}

/**************************************************************************/

static void write_imported(const char* stream_nm, const DiscreteArrayList& arr_c, const std::vector<std::string>& genes,
  const std::vector<std::string>& conds, const DiscreteArray& symbols) {
  FILE* fw = mustOpenWrite(stream_nm);
  fprintf(fw, "o");
  for (std::size_t col = 0; col < conds.size(); col++)
    fprintf(fw, "\t%s", conds[col].c_str());
  fputc('\n', fw);
  for (std::size_t row = 0; row < genes.size(); row++) {
    fprintf(fw, "%s", genes[row].c_str());
    for (std::size_t col = 0; col < conds.size(); col++)
      fprintf(fw, "\t%d", symbols[arr_c[row][col]]);
    fputc('\n', fw);
  }
  fclose(fw);
}

/* Identified clusters are backtraced to the original data, by
* putting the clustered vectors together, identify common column
*/
static void print_bc(FILE* fw, const Block& b, const int& num,
  const DiscreteArrayList& arr_c, const std::vector<std::string>& genes, const std::vector<std::string>& conds,
  const Symbols& symbols) {
  int block_rows, block_cols;
  int num_1 = 0, num_2 = 0;
  /* block height (genes) */
  block_rows = b.block_rows();
  /* block_width (conditions) */
  block_cols = b.block_cols();
  fprintf(fw, "BC%03d\tS=%d\tPvalue:%g \n", num, block_rows * block_cols, static_cast<double>(b.pvalue));
  /* fprintf(fw, "BC%03d\tS=%d\tPvalue:%lf \n", num, block_rows * block_cols, (double)b.pvalue); */
  fprintf(fw, " Genes [%d]: ", block_rows);
  for (std::set<int>::iterator it = b.genes_order.begin(); it != b.genes_order.end(); ++it)
    fprintf(fw, "%s ", genes[*it].c_str());
  for (std::set<int>::iterator it = b.genes_reverse.begin(); it != b.genes_reverse.end(); ++it)
    fprintf(fw, "%s ", genes[*it].c_str());
  fprintf(fw, "\n");
  fprintf(fw, " Conds [%d]: ", block_cols);
  for (std::set<int>::iterator it = b.conds.begin(); it != b.conds.end(); ++it)
    fprintf(fw, "%s ", conds[*it].c_str());
  fprintf(fw, "\n");
  /* the complete block data output */
  for (std::set<int>::iterator it = b.genes_order.begin(); it != b.genes_order.end(); ++it) {
    fprintf(fw, "%10s:", genes[*it].c_str());
    for (std::set<int>::iterator jt = b.conds.begin(); jt != b.conds.end(); ++jt) {
      fprintf(fw, "\t%d", symbols[arr_c[*it][*jt]]);
      if (it == b.genes_order.begin()) {
        if (symbols[arr_c[*it][*jt]] == +1) num_1++;
        if (symbols[arr_c[*it][*jt]] == -1) num_2++;
      }
    }
    fputc('\n', fw);
  }
  fputc('\n', fw);
  for (std::set<int>::iterator it = b.genes_reverse.begin(); it != b.genes_reverse.end(); ++it) {
    fprintf(fw, "%10s:", genes[*it].c_str());
    for (std::set<int>::iterator jt = b.conds.begin(); jt != b.conds.end(); ++jt)
      fprintf(fw, "\t%d", symbols[arr_c[*it][*jt]]);
    fputc('\n', fw);
  }
  /*fprintf(stdout, "BC%03d: #of 1 and -1 are:\t%d\t%d\n",num,num_1,num_2);
  fputc('\n', fw);*/
}

void print_params(FILE* fw, const std::size_t col_width, double filter, double TOLERANCE, int RPT_BLOCK) {
  fprintf(fw, "# QUBIC version %.1f output\n", 1.9);
  fprintf(fw, "# \n");
  fprintf(fw, "# Parameters: -k %d -f %.2f -c %.2f -o %d",
    static_cast<unsigned int>(col_width), filter, TOLERANCE, RPT_BLOCK);
  fprintf(fw, "\n\n");
}

/**************************************************************************/

void write_blocks(const std::string& tfile, const std::vector<std::string>& row_names, const std::vector<std::string>& col_names, const double c, const int o, const double filter, int col_width, DiscreteArrayListWithSymbols all, std::vector<Block> output, const bool verbose) {
  FILE* fw = mustOpenWrite((tfile + ".blocks").c_str());
  print_params(fw, col_width, filter, c, o);
  for (std::size_t i = 0; i < output.size(); i++)
    print_bc(fw, output[i], i, all.list, row_names, col_names, all.symbols);
  /* clean up */
  fclose(fw);
  if (verbose) fprintf(stdout, "%d clusters are written to %s\n", static_cast<unsigned int>(output.size()),
    (tfile + ".blocks").c_str());
}

void write_chars(const std::string& tfile, const std::vector<std::string>& row_names, const std::vector<std::string>& col_names, DiscreteArrayListWithSymbols all, const bool verbose) {
  write_imported((tfile + ".chars").c_str(), all.list, row_names, col_names, all.symbols);
  if (verbose) fprintf(stdout, "Formatted data are written to %s\n", (tfile + ".chars").c_str());
}

std::size_t fix_col_width(const std::vector<std::vector<short>>& x, const int k) {
  return k == 2 ? std::max(x[0].size() / 20, static_cast<std::size_t>(2)) : k;
}

std::vector<Block> r_main(const std::vector<std::vector<short>>& short_matrix, const double c, const int o,
  const double filter, const int k, const Option& option, const bool verbose) {
  std::size_t col_width = fix_col_width(short_matrix, k);
  if (verbose) fprintf(stdout, "Size of matrix is (%lu, %lu)\n", static_cast<unsigned long>(short_matrix.size()), static_cast<unsigned long>(short_matrix[0].size()));
  DiscreteArrayListWithSymbols all = make_charsets_d(short_matrix, verbose);
  return qubic::init_qubic(all, c, filter, col_width, o, option, verbose);
}

std::vector<Block> r_main(const std::vector<std::vector<short>>& short_matrix, const double c, const int o,
  const double filter, const int k, const Option& option, const bool verbose, const std::vector<std::vector<float>>& weight_matrix) {
  std::size_t col_width = fix_col_width(short_matrix, k);
  if (verbose) fprintf(stdout, "Size of matrix is (%lu, %lu)\n", static_cast<unsigned long>(short_matrix.size()), static_cast<unsigned long>(short_matrix[0].size()));
  DiscreteArrayListWithSymbols all = make_charsets_d(short_matrix, verbose);
  return qubic::init_qubic(all, c, filter, col_width, o, option, verbose, weight_matrix);
}

std::vector<Block> main_d(const std::vector<std::vector<short>>& short_matrix, const std::vector<std::string>& row_names,
  const std::vector<std::string>& col_names, const std::string& tfile,
  const double c, const int o, const double filter, const int k, const Option& option, const bool verbose) {
  std::size_t col_width = fix_col_width(short_matrix, k);
  DiscreteArrayListWithSymbols all = make_charsets_d(short_matrix, verbose);
  std::vector<Block> output = qubic::init_qubic(all, c, filter, col_width, o, option, verbose);
  write_chars(tfile, row_names, col_names, all, verbose);
  write_blocks(tfile, row_names, col_names, c, o, filter, col_width, all, output, verbose);
  return output;
}

void write_rules(const std::string& tfile, const std::vector<std::string>& row_names, std::vector<rule> genes_rules, const bool verbose) {
  FILE* fw = mustOpenWrite((tfile + ".rules").c_str());
  for (std::size_t row = 0; row < row_names.size(); row++) {
    fprintf(fw, "row %s :low=%2.5f, up=%2.5f; %d down-regulated,%d up-regulated\n", row_names[row].c_str(),
      genes_rules[row].lower, genes_rules[row].upper, static_cast<unsigned int>(genes_rules[row].cntl),
      static_cast<unsigned int>(genes_rules[row].cntu));
  }
  fclose(fw);
  if (verbose) fprintf(stdout, "Discretization rules are written to %s\n", (tfile + ".rules").c_str());
}

std::vector<Block> main_c(const std::vector<std::vector<float>>& float_matrix, const std::vector<std::string>& row_names,
  const std::vector<std::string>& col_names, const std::string& tfile, const short r, const double q,
  const double c, const int o, const double filter, const int col_width, const Option& option, const bool verbose) {
  std::vector<rule> genes_rules;
  DiscreteArrayList short_matrix = discretize(float_matrix, r, q, genes_rules);
  write_rules(tfile, row_names, genes_rules, verbose);
  return main_d(short_matrix, row_names, col_names, tfile, c, o, filter, col_width, option, verbose);
}