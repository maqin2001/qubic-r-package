#include "qubic.h"

#include <cassert>
#include <cmath> // floor
#include <cstring> // memset
#include <algorithm>

#include "edge_list.h"
#include "version.h"

struct rule {
  float lower;
  float upper;
  size_t cntl;
  size_t cntu;
};

class qubic {
  /* global data */
  std::vector<std::vector<continuous> > arr;
public:
  DiscreteArrayList arr_c;
  std::vector<discrete> symbols;
  size_t COL_WIDTH;
  std::vector<rule> genes_rules;
private:
  bool IS_DISCRETE;
  double FILTER;
  double TOLERANCE;
  int RPT_BLOCK;
  double QUANTILE;
  discrete DIVIDED;
  size_t rows, cols;
  Prog_options po;

  /* must be defined */
  std::vector<std::vector<bits16> > profile;

  int r_puts() {
    puts("\n===================================================================\n"
      "[Usage]\n"
      "qubic(data, [argument list]);\n"
      "like :\n"
      "qubic(data, file = 'rQUBIC', q = 0.06, c = 0.95, f = 1, k = 2, r = 1, o = 100, d = 'F')\n"
      "===================================================================\n"
      "[Input]\n"
      "-file : input file must be one of two tab-delimited formats\n"
      "  A) continuous data (default, use pre-set discretization (see -q and -r))\n"
      "     -------------------------------------\n"
      "     o        cond1    cond2    cond3\n"
      "     gene1      2.4      3.5     -2.4\n"
      "     gene2     -2.1      0.0      1.2\n"
      "     -------------------------------------\n"
      "  B) discrete data with arbitrary classes (turn on -d)\n"
      "     use '0' for missing or insignificant data\n"
      "     -------------------------------------\n"
      "     o        cond1    cond2    cond3\n"
      "     gene1        1        2        2\n"
      "     gene2       -1        2        0\n"
      "     -------------------------------------\n"
      "-q : use quantile discretization for continuous data\n"
      "     default: 0.06 (see details in Method section in paper)\n"
      "-r : the number of ranks as which we treat the up(down)-regulated value\n"
      "     when discretization\n"
      "     default: 1\n"
      "-d : discrete data, where user should send their processed data\n"
      "     to different value classes, see above\n"
      "-C : the flag using the lower bound of condition number (5 percents of the gene number)\n"
      "===================================================================\n"
      "[Output]\n"
      "-o : number of blocks to report, default: 100\n"
      "-f : filtering overlapping blocks,\n"
      "     default: 1 (do not remove any blocks)\n"
      "-k : minimum column width of the block,\n"
      "     default: 5% of columns, minimum 2 columns\n"
      "-c : consistency level of the block (0.5-1.0], the minimum ratio between the\n"
      "     number of identical valid symbols in a column and the total number \n"
      "     of rows in the output\n"
      "     default: 0.95\n"
      "===================================================================\n");
    return 1;
  }

  static discrete charset_add(std::vector<discrete>& ar, const discrete & s, discrete *bb) {
    /*A signed short can hold all the values between SHRT_MIN  and SHRT_MAX inclusive.SHRT_MIN is required to be -32767 or less,SHRT_MAX must be at least 32767*/
    int ps = s + SHRT_MAX;
    if (bb[ps] < 0) {
      bb[ps] = static_cast<discrete>(ar.size());
      ar.push_back(s);
    }
    return bb[ps];
  }

  static continuous quantile_from_sorted_data(const std::vector<continuous> & sorted_data, const size_t n, const double f) {
    /*floor function returns the largest integral value less than or equal to x*/
    int i = static_cast <int> (floor((n - 1) * f));
    continuous delta = static_cast<continuous>((n - 1) * f - i);
    return (1 - delta)*sorted_data[i] + delta*sorted_data[i + 1];
  }

  static discrete dis_value(const float current, const discrete divided, const std::vector<continuous> & small, const int cntl, const std::vector<continuous> & big, const int cntu) {
    continuous d_space = static_cast<continuous>(1.0 / divided);
    for (discrete i = 0; i < divided; i++) {
      if ((cntl > 0) && (current <= quantile_from_sorted_data(small, cntl, static_cast<continuous>(d_space * (i + 1)))))
        return -i - 1;
      if ((cntu > 0) && (current >= quantile_from_sorted_data(big, cntu, static_cast<continuous>(1.0 - d_space * (i + 1)))))
        return i + 1;
    }
    return 0;
  }

  void discretize(const std::vector<std::vector<continuous> > &arr,
    discrete *bb, std::vector<discrete>& symbols, const double f, DiscreteArrayList &arr_c, const discrete divided) {
    size_t row, col;
    std::vector<continuous> rowdata(arr[0].size());
    std::vector<continuous> big(arr[0].size()), small(arr[0].size());
    size_t i, cntu, cntl;
    float f1, f2, f3, upper, lower;
    for (row = 0; row < arr.size(); row++) {
      for (col = 0; col < arr[0].size(); col++)
        rowdata[col] = arr[row][col];
      std::sort(rowdata.begin(), rowdata.end());

      f1 = quantile_from_sorted_data(rowdata, arr[0].size(), 1 - f);
      f2 = quantile_from_sorted_data(rowdata, arr[0].size(), f);
      f3 = quantile_from_sorted_data(rowdata, arr[0].size(), 0.5);
      if ((f1 - f3) >= (f3 - f2)) {
        upper = 2 * f3 - f2;
        lower = f2;
      } else {
        upper = f1;
        lower = 2 * f3 - f1;
      }
      cntu = 0; cntl = 0;
      for (i = 0; i < arr[0].size(); i++) {
        if (rowdata[i] < lower) {
          small[cntl] = rowdata[i];
          cntl++;
        }
        if (rowdata[i] > upper) {
          big[cntu] = rowdata[i];
          cntu++;
        }
      }
      for (col = 0; col < arr[0].size(); col++)
        arr_c[row][col] = charset_add(symbols, dis_value(arr[row][col], divided, small, cntl, big, cntu), bb);

      rule rule;
      rule.lower = lower;
      rule.upper = upper;
      rule.cntl = cntl;
      rule.cntu = cntu;

      genes_rules.push_back(rule);
    }


  }

  void seed_update(const DiscreteArray& s) {
    for (size_t i = 0; i < cols; i++)
      profile[i][s[i]]++;
  }

  /* scan through all columns and identify the set within threshold,
   * "fuzziness" of the block is controlled by TOLERANCE (-c)
   */
  void scan_block(const std::vector<int> & gene_set, Block & b) {
    size_t i, j;
    size_t block_rows, cur_rows;
    block_rows = cur_rows = gene_set.size();

    size_t k;
    for (j = 0; j < cols; j++)
      for (k = 0; k < symbols.size(); k++)
        profile[j][k] = 0;
    for (j = 0; j < cur_rows; j++)
      seed_update(arr_c[gene_set[j]]);

    int btolerance = static_cast<int>(ceil(TOLERANCE* block_rows));
    for (j = 0; j < cols; j++) {
      /* See if this column satisfies tolerance */
      /* here i start from 1 because symbols[0]=0 */
      for (i = 1; i < symbols.size(); i++) {
        if ((profile[j][i] >= btolerance)) {
          bool result = b.conds.insert(j).second;
          assert(result);
          break;
        }
      }
    }
  }
  /*************************************************************************/



  void update_colcand(std::vector<bool> &colcand, const DiscreteArray &g1, const DiscreteArray &g2) {
    size_t i;
    for (i = 0; i < cols; i++)
      if (colcand[i] && (g1[i] != g2[i]))
        colcand[i] = false;
  }

  /*calculate the weight of the edge with two vertices g1 and g2*/
  int intersect_row(const std::vector<bool> &colcand, const DiscreteArray &g1, const DiscreteArray &g2) {
    size_t i;
    int cnt = 0;
    for (i = 0; i < cols; i++)
      if (colcand[i] && (g1[i] == g2[i]) && (g1[i] != 0))
        cnt++;
    return cnt;
  }

  /*calculate the negative correlation between g1 and g2*/
  int reverse_row(const std::vector<bool> &colcand, const DiscreteArray &g1, const DiscreteArray &g2) {
    size_t i;
    int cnt = 0;
    for (i = 0; i < cols; i++) {
      if (colcand[i] && (symbols[g1[i]] == -symbols[g2[i]])) cnt++;
    }
    return cnt;
  }

  /* calculate the coverage of any row to the current consensus
  * cnt = # of valid consensus columns
  */
  int seed_current_modify(const DiscreteArray &s, std::vector<bool> &colcand, const int components) {
    size_t i, k, flag, n;
    int threshold = static_cast <int> (ceil(components * TOLERANCE));
    discrete ss;
    int cnt = 0;
    for (i = 0; i < cols; i++) {
      flag = 0; ss = s[i];
      for (k = 1; k < symbols.size(); k++) {
        n = profile[i][k];
        if (k == ss) n++;
        if (n >= threshold) {
          flag = k;
          break;
        }
      }
      if (flag) {
        cnt++; colcand[i] = true;
      }
    }
    return cnt;
  }

  /*check whether current edge can be treat as a seed*/
  bool check_seed(const Edge *e, const std::vector<Block> &bb) {
    int block_id = bb.size();
    std::vector<int> profiles(rows);
    int i, b1, b2, b3;
    bool fg = false;
    b1 = b2 = -1;
    for (i = 0; i < block_id; i++)
      if ((bb[i].contains(e->gene_one)) && (bb[i].contains(e->gene_two)))
        return false;

    for (i = 0; i < rows; i++) profiles[i] = 0;
    fg = false;
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
    if ((b1 == -1) || (b2 == -1))
      return true;
    else {
      for (std::set<int>::iterator it = bb[b1].genes_order.begin(); it != bb[b1].genes_order.end(); ++it)
        profiles[*it]++;
      for (std::set<int>::iterator it = bb[b1].genes_reverse.begin(); it != bb[b1].genes_reverse.end(); ++it)
        profiles[*it]++;
      for (std::set<int>::iterator it = bb[b2].genes_order.begin(); it != bb[b2].genes_order.end(); ++it)
        profiles[*it]++;
      for (std::set<int>::iterator it = bb[b2].genes_reverse.begin(); it != bb[b2].genes_reverse.end(); ++it)
        profiles[*it]++;

      for (i = 0; i < rows; i++)
        if (profiles[i] > 1) return false;
      b3 = std::max(bb[b1].block_cols(), bb[b2].block_cols());
      return !(e->score < b3/* (bb[b1]->block_cols + bb[b2]->block_cols) / 2*/);
    }
  }

  long double get_pvalue(const continuous & a, const int & b) {
    int i = 0;
    long double one = 1, pvalue = 0;
    long double poisson = one / exp(a);
    for (i = 0; i < b + 300; i++) {
      if (i >(b - 1)) pvalue = pvalue + poisson;
      else poisson = poisson * a / (i + 1);
    }
    return pvalue;
  }

  void block_init(Block & b,
    std::vector<int> &genes, std::vector<int> &scores,
    std::vector<bool> &candidates, const int &cand_threshold,
    size_t &components, std::vector<long double> & pvalues) {
    int i, score, top;
    int cnt = 0, cnt_all = 0, pid = 0;
    continuous cnt_ave = 0, row_all = static_cast<continuous>(rows);
    long double pvalue;
    int max_cnt, max_i;
    std::vector<int> arr_rows(rows), arr_rows_b(rows);
    std::vector<bool> colcand(cols);
    std::fill(colcand.begin(), colcand.end(), false);
    DiscreteArray g1, g2;

    g1 = arr_c[genes[0]];
    g2 = arr_c[genes[1]];
    for (i = 0; i < cols; i++)
      if ((g1[i] == g2[i]) && (g1[i] != 0))
        colcand[i] = true;

    for (i = 0; i < rows; i++) {
      arr_rows[i] = intersect_row(colcand, arr_c[genes[0]], arr_c[i]);
      arr_rows_b[i] = arr_rows[i];
    }
    /*we just get the largest 100 rows when we initial a bicluster because we believe that
     * the 100 rows can characterize the structure of the bicluster
     * btw, it can reduce the time complexity*/
    if (rows > 100) {
      std::sort(arr_rows_b.begin(), arr_rows_b.end());
      top = arr_rows_b[rows - 100];
      for (i = 0; i < rows; i++)
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
      for (i = 0; i < rows; i++) {
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
      if (po.IS_cond) {
        if (max_cnt < COL_WIDTH || max_i < 0 || max_cnt < b.cond_low_bound) break;
      } else {
        if (max_cnt < COL_WIDTH || max_i < 0) break;
      }
      if (po.IS_area)	score = components * max_cnt;
      else score = std::min(static_cast<int>(components), max_cnt); // ggggggggg
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
  static bool block_cmpr(const Block &a, const Block &b) {
    return a.score > b.score;
  }

  /************************************************************************/

  std::vector<Block> report_blocks(std::vector<Block> bb) {
    std::vector<Block> output;

    int num = bb.size();

    std::sort(bb.begin(), bb.end(), block_cmpr);

    int i, j, k;
    /*MIN MAX et al functions can be accessed in struct.h*/
    int n = std::min(num, RPT_BLOCK);
    bool flag;

    double cur_rows, cur_cols;
    double inter_rows, inter_cols;
    /*double proportion;*/

    /* the major post-processing here, filter overlapping blocks*/
    i = 0; j = 0;
    while (i < num && j < n) {
      Block & b_ptr = bb[i];
      cur_rows = b_ptr.block_rows();
      cur_cols = b_ptr.block_cols();

      flag = true;
      k = 0;
      while (k < j) {
        inter_rows = count_intersect(output[k].genes_order, b_ptr.genes_order) +
          count_intersect(output[k].genes_order, b_ptr.genes_reverse) +
          count_intersect(output[k].genes_reverse, b_ptr.genes_order) +
          count_intersect(output[k].genes_reverse, b_ptr.genes_reverse);
        inter_cols = count_intersect(output[k].conds, b_ptr.conds);

        if (inter_rows * inter_cols > FILTER*cur_rows*cur_cols) {
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


  std::vector<Block> cluster(const std::vector<Edge *> &el) {
    std::vector<Block> bb;

    size_t j, k, components;

    profile.resize(cols, std::vector<bits16>(symbols.size()));

    std::vector<long double> pvalues(rows);
    std::vector<bool> candidates(rows);

    std::set<int> allincluster;

    for (std::vector<Edge *>::const_iterator it = el.begin(); it != el.end(); ++it) {
      const Edge *e = *it;
      /* check if both genes already enumerated in previous blocks */
      bool flag = true;
      /* speed up the program if the rows bigger than 200 */
      if (rows > 250) {
        if (allincluster.find(e->gene_one) != allincluster.end() && allincluster.find(e->gene_two) != allincluster.end()) flag = false;
      } else {
        flag = check_seed(e, bb);
      }
      if (!flag) continue;

      for (j = 0; j < cols; j++)
        for (k = 0; k < symbols.size(); k++)
          profile[j][k] = 0;

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
      int cand_threshold = static_cast<int>(floor(COL_WIDTH * TOLERANCE));
      if (cand_threshold < 2)
        cand_threshold = 2;

      /* maintain a candidate list to avoid looping through all rows */
      for (j = 0; j < rows; j++)
        candidates[j] = true;
      candidates[e->gene_one] = candidates[e->gene_two] = false;
      components = 2;

      /* expansion step, generate a bicluster without noise */
      block_init(b, genes_order, scores, candidates, cand_threshold, components, pvalues);

      /* track back to find the genes by which we get the best score*/
      for (k = 0; k < components; k++) {
        if (po.IS_pvalue)
          if ((pvalues[k] == b.pvalue) && (k >= 2) && (scores[k] != scores[k + 1])) break;
        if ((scores[k] == b.score) && ((k + 1 == scores.size()) || (scores[k + 1] != b.score))) break;
      }
      components = k + 1;
      std::fill(candidates.begin(), candidates.end(), true);

      for (int ki = 0; ki < components - 1; ki++) {
        seed_update(arr_c[genes_order[ki]]);
        candidates[genes_order[ki]] = false;
      }
      candidates[genes_order[k]] = false;
      genes_order.resize(k + 1);

      std::vector<bool> colcand(cols);
      std::fill(colcand.begin(), colcand.end(), false);

      /* add columns satisfy the conservative r */
      int cnt = seed_current_modify(arr_c[genes_order[k]], colcand, components);

      /* add some new possible genes */
      int m_cnt;
      for (size_t ki = 0; ki < rows; ki++) {
        m_cnt = intersect_row(colcand, arr_c[genes_order[0]], arr_c[ki]);
        if (candidates[ki] && (m_cnt >= floor(cnt* TOLERANCE))) {
          genes_order.push_back(ki);
          components++;
          candidates[ki] = false;
        }
      }

      /* add genes that negative regulated to the consensus */
      for (size_t ki = 0; ki < rows; ki++) {
        m_cnt = reverse_row(colcand, arr_c[genes_order[0]], arr_c[ki]);
        if (candidates[ki] && (m_cnt >= floor(cnt * TOLERANCE))) {
          genes_reverse.push_back(ki);
          components++;
          candidates[ki] = false;
        }
      }

      /* store gene arrays inside block */

      scan_block(genes_order, b);
      if (b.block_cols() == 0) continue;
      //b.block_rows = components;
      if (po.IS_pvalue) b.score = static_cast<int>(-(100 * log(b.pvalue)));
      else b.score = components * b.block_cols();

      for (std::vector<int>::iterator it = genes_order.begin(); it != genes_order.end(); ++it) {
        bool result = b.genes_order.insert(*it).second;
        assert(result);
        allincluster.insert(*it);
      }

      for (std::vector<int>::iterator it = genes_reverse.begin(); it != genes_reverse.end(); ++it) {
        bool result = b.genes_reverse.insert(*it).second;
        assert(result);
        allincluster.insert(*it);
      }

      /*save the current block b to the block list bb so that we can sort the blocks by their score*/
      bb.push_back(b);
      /* reaching the results number limit */
      if (bb.size() == po.SCH_BLOCK) break;
      fputc('.', stdout);
    }
    fprintf(stdout, "\n");
    return report_blocks(bb);
  }

  /*make_graph subroutine prototypes */

  /* remove a row from the profile */
  void seed_deduct(const discrete *s) {
    for (size_t i = 0; i < cols; i++) {
      profile[i][s[i]]--;
    }
  }

  std::vector<Block> make_graph(const DiscreteArrayList &arr_c, size_t &COL_WIDTH) {
    EdgeList EdgeList(arr_c, COL_WIDTH);

    /* bi-clustering */
    fprintf(stdout, "Clustering started");
    return cluster(EdgeList.get_edge_list());
  }

  static int intersect_rowE(const std::vector<bool> &colcand, std::vector<discrete> &g1, std::vector<discrete> &g2, const int cols) {
    int i, cnt = 0;
    for (i = 0; i < cols; i++)
      if (colcand[i] && (g1[i] == g2[i]) && g1[i] != 0)
        cnt++;
    return cnt;
  }

  int reverse_rowE(const std::vector<bool> &colcand, std::vector<discrete> &g1, std::vector<discrete> &g2, const int cols) {
    int i, cnt = 0;
    for (i = 0; i < cols; i++)
      if (colcand[i] && (symbols[g1[i]] == -symbols[g2[i]]))
        cnt++;
    return cnt;
  }

  void make_charsets(std::vector<discrete> &symbols) {
    discrete bb[USHRT_MAX];
    memset(bb, -1, USHRT_MAX*sizeof(*bb));
    charset_add(symbols, 0, bb);
    if (IS_DISCRETE) {
      for (size_t i = 0; i < arr.size(); i++)
        for (size_t j = 0; j < arr[0].size(); j++) {
          arr_c[i][j] = charset_add(symbols, (discrete)arr[i][j], bb);
        }
      fprintf(stdout, "Discretized data contains %d classes with charset [ ", static_cast<unsigned int>(symbols.size()));
      for (size_t i = 0; i < symbols.size(); i++)
        fprintf(stdout, "%d ", symbols[i]);  fprintf(stdout, "]\n");
    } else {
      for (size_t i = 0; i < arr.size(); i++)
        for (size_t j = 0; j < arr[0].size(); j++)
          arr_c[i][j] = 0;
      discretize(arr, bb, symbols, QUANTILE, arr_c, DIVIDED);
    }
  }

  std::vector<Block> run_qubic(const double rq, const double rc, const double rf, const int rk, const discrete rr, const int ro, const bool rd) {
    arr_c.resize(rows, DiscreteArray(cols));

    fprintf(stdout, "\nQUBIC %s: greedy biclustering\n\n", VER);
    /* get the program options defined in get_options.c */
    /*set memory for the point which is declared in struct.h*/
    //AllocVar(po);
    /*Initialize the point*/
    IS_DISCRETE = rd;
    po.IS_pvalue = false;
    /* case 'P': po.IS_pvalue = true; */
    COL_WIDTH = rk;
    DIVIDED = rr;
    QUANTILE = rq;
    TOLERANCE = rc;
    RPT_BLOCK = ro;
    po.SCH_BLOCK = 2 * RPT_BLOCK;
    /* ensure enough searching space */
    /*if (po.SCH_BLOCK < 1000) po.SCH_BLOCK = 1000;*/
    FILTER = rf;
    /* case 's': po.IS_SWITCH = true; */
    po.IS_area = false;
    /* case 'S': po.IS_area = true; */
    po.IS_cond = false;
    /* case 'C': po.IS_cond = true; */

    make_charsets(symbols);

    /* the file that stores all blocks */
    return make_graph(arr_c, COL_WIDTH);
  }

public:
  std::vector<Block> init_qubic(const double rq = 0.06, const double rc = 0.95, const double rf = 1, const int rk = 2, const discrete rr = 1, const int ro = 100, const bool rd = false) {
    return run_qubic(rq, rc, rf, rk, rr, ro, rd);
  }

  qubic(const std::vector<std::vector<float> > &data) {
    if (data.size() == 0) throw - 1;
    rows = data.size();
    cols = data[0].size();

    arr = data;
  }
};

/**************************************************************************/
/* file-related operations */

/* Strings */
/* strcmp: a zero value indicates that both strings are equal.
* a value greater than zero indicates that the first character that does not match has a greater value in str1 than in str2;
* And a value less than zero indicates the opposite.
*/
#define sameString(a, b) (strcmp((a), (b))==0)
/* Returns TRUE if two strings are same */

FILE *mustOpen(const char *fileName, const char *mode)
/* Open a file or die */
{
  FILE *f;

  if (sameString(fileName, "stdin")) return stdin;
  if (sameString(fileName, "stdout")) return stdout;
  if ((f = fopen(fileName, mode)) == NULL) {
    const char *modeName = "";
    if (mode) {
      if (mode[0] == 'r') modeName = " to read";
      else if (mode[0] == 'w') modeName = " to write";
      else if (mode[0] == 'a') modeName = " to append";
    }
    fprintf(stderr, "[Error] Can't open %s%s: %s", fileName, modeName, 0);
    throw - 1;
  }
  return f;
}

/**************************************************************************/

static void write_imported(const char* stream_nm, const DiscreteArrayList &arr_c, const std::vector<std::string> &genes, const std::vector<std::string> &conds, const std::vector<discrete> &symbols) {
  size_t row, col;
  FILE *fw;
  fw = mustOpen(stream_nm, "w");
  fprintf(fw, "o");
  for (col = 0; col < conds.size(); col++)
    fprintf(fw, "\t%s", conds[col].c_str());
  fputc('\n', fw);
  for (row = 0; row < genes.size(); row++) {
    fprintf(fw, "%s", genes[row].c_str());
    for (col = 0; col < conds.size(); col++)
      fprintf(fw, "\t%d", symbols[arr_c[row][col]]);
    fputc('\n', fw);
  }
  fprintf(stdout, "Formatted data are written to %s\n", stream_nm);
  fclose(fw);
}

/* Identified clusters are backtraced to the original data, by
* putting the clustered vectors together, identify common column
*/
static void print_bc(FILE* fw, const Block & b, const int & num,
  const DiscreteArrayList & arr_c, const std::vector<std::string> &genes, const std::vector<std::string> &conds, const std::vector<discrete> &symbols) {
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
        if (symbols[arr_c[*it][*jt]] == 1) num_1++;
        if (symbols[arr_c[*it][*jt]] == -1) num_2++;
      }
    }
    fputc('\n', fw);
  }
  fputc('\n', fw);
  for (std::set<int>::iterator it = b.genes_reverse.begin(); it != b.genes_reverse.end(); ++it) {
    fprintf(fw, "%10s:", genes[*it].c_str());
    for (std::set<int>::iterator jt = b.conds.begin(); jt != b.conds.end(); ++jt) {
      fprintf(fw, "\t%d", symbols[arr_c[*it][*jt]]);
    }
    fputc('\n', fw);
  }
  /*fprintf(stdout, "BC%03d: #of 1 and -1 are:\t%d\t%d\n",num,num_1,num_2);
  fputc('\n', fw);*/
}

void print_params(FILE *fw, bool IS_DISCRETE, const std::string &FN, const size_t COL_WIDTH, double FILTER, double TOLERANCE,
  int RPT_BLOCK,
  double QUANTILE,
  discrete DIVIDED) {
  std::string filedesc = "continuous";
  if (IS_DISCRETE)
    filedesc = "discrete";
  fprintf(fw, "# QUBIC version %.1f output\n", 1.9);
  fprintf(fw, "# Datafile %s: %s type\n", FN.c_str(), filedesc.c_str());
  fprintf(fw, "# Parameters: -k %d -f %.2f -c %.2f -o %d",
    static_cast<unsigned int>(COL_WIDTH), FILTER, TOLERANCE, RPT_BLOCK);
  if (!IS_DISCRETE)
    fprintf(fw, " -q %.2f -r %d", QUANTILE, DIVIDED);
  fprintf(fw, "\n\n");
}

/**************************************************************************/

std::vector<Block> r_main(const std::vector<std::vector<float> > &data, const std::vector<std::string > &row_names, const std::vector<std::string > &col_names, const std::string & tfile = "rQUBIC", const double rq = 0.06, const double rc = 0.95, const double rf = 1, const int rk = 2, const discrete rr = 1, const int ro = 100, const bool rd = false) {
  qubic qubic(data);
  std::vector<Block> output = qubic.init_qubic(rq, rc, rf, rk, rr, ro, rd);

  {
    FILE *fw = mustOpen((tfile + ".rules").c_str(), "w");
    for (size_t row = 0; row < row_names.size(); row++) {
      fprintf(fw, "row %s :low=%2.5f, up=%2.5f; %d down-regulated,%d up-regulated\n", row_names[row].c_str(), qubic.genes_rules[row].lower, qubic.genes_rules[row].upper, static_cast<unsigned int>(qubic.genes_rules[row].cntl), static_cast<unsigned int>(qubic.genes_rules[row].cntu));
    }
    fprintf(stdout, "Discretization rules are written to %s\n", (tfile + ".rules").c_str());
    fclose(fw);
  }

  write_imported((tfile + ".chars").c_str(), qubic.arr_c, row_names, col_names, qubic.symbols);

  {
    FILE *fw = mustOpen((tfile + ".blocks").c_str(), "w");
    print_params(fw, rd, tfile, qubic.COL_WIDTH, rf, rc, ro, rq, rr);
    for (size_t i = 0; i < output.size(); i++) {
      print_bc(fw, output[i], i, qubic.arr_c, row_names, col_names, qubic.symbols);
    }
    /* clean up */
    fclose(fw);
  }
  fprintf(stdout, "%d clusters are written to %s\n", static_cast<unsigned int>(output.size()), (tfile + ".blocks").c_str());
  return output;
}

std::vector<Block> r_main(const std::vector<std::vector<float> > &data, const std::string &tfile, const double rq, const double rc, const double rf, const int rk, const short rr, const int ro, const bool rd) {
  qubic qubic(data);
  return qubic.init_qubic(rq, rc, rf, rk, rr, ro, rd);
}

std::vector<Block> r_main(const std::vector<std::vector<float> > &x, const short r, const double q, const double c, const int o, const double f, const int rk, const bool rd) {
  qubic qubic(x);
  return qubic.init_qubic(q, c, f, rk, r, o, rd);
}


std::vector<Block> r_main(const std::vector<std::vector<float> > &data) {
  qubic qubic(data);
  return qubic.init_qubic();
}
