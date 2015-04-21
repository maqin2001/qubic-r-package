#include "struct.h"

#include <algorithm>
#include <bitset>
#include <vector>
#include <set>
#include <memory>
#include <string>

#include <assert.h>

#include "edge_list.h"
#include "version.h"

class qubic {
  /* global data */
  std::vector<std::vector<continuous> > arr;
  DiscreteArrayList arr_c;
  std::vector<discrete> symbols;
  std::vector<std::string> genes;
  std::vector<std::string> conds;
  std::vector<std::string> sub_genes;
  std::vector<bool> sublist;
  size_t rows, cols;
  int TFindex;
  int sub_genes_row;

  Prog_options po;

  /* must be defined */
  std::vector<std::vector<bits16> > profile;
  int col_width;
#define MAXC 100000
  
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
      bb[ps] = ar.size();
      ar.push_back(s);
    }
    return bb[ps];
  }

  static continuous quantile_from_sorted_data(const std::vector<continuous> & sorted_data, const size_t & n, const double & f) {
    /*floor function returns the largest integral value less than or equal to x*/
    int i = static_cast <int> (floor((n - 1) * f));
    continuous delta = (n - 1) * f - i;
    return (1 - delta)*sorted_data[i] + delta*sorted_data[i + 1];
  }

  static discrete dis_value(const float & current, const discrete & divided, const std::vector<continuous> & small, const int & cntl, const std::vector<continuous> & big, const int & cntu) {
    continuous d_space = static_cast<continuous>(1.0 / divided);
    for (discrete i = 0; i < divided; i++) {
      if ((cntl > 0) && (current <= quantile_from_sorted_data(small, cntl, static_cast<continuous>(d_space * (i + 1)))))
        return -i - 1;
      if ((cntu > 0) && (current >= quantile_from_sorted_data(big, cntu, static_cast<continuous>(1.0 - d_space * (i + 1)))))
        return i + 1;
    }
    return 0;
  }

  static void discretize(const char* stream_nm, const std::vector<std::vector<continuous> > &arr, 
    discrete *bb, std::vector<discrete>& symbols, const double& f, DiscreteArrayList &arr_c, const discrete divided,
    const std::vector<std::string>& genes) {
    size_t row, col;
    std::vector<continuous> rowdata(arr[0].size());
    std::vector<continuous> big(arr[0].size()), small(arr[0].size());
    size_t i, cntu, cntl;
    float f1, f2, f3, upper, lower;
    FILE *fw = mustOpen(stream_nm, "w");
    for (row = 0; row < arr.size(); row++)
    {
      for (col = 0; col < arr[0].size(); col++)
        rowdata[col] = arr[row][col];
      std::sort(rowdata.begin(), rowdata.end());

      f1 = quantile_from_sorted_data(rowdata, arr[0].size(), 1 - f);
      f2 = quantile_from_sorted_data(rowdata, arr[0].size(), f);
      f3 = quantile_from_sorted_data(rowdata, arr[0].size(), 0.5);
      if ((f1 - f3) >= (f3 - f2)) {
        upper = 2 * f3 - f2;
        lower = f2;
      }
      else {
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
      fprintf(fw, "row %s :low=%2.5f, up=%2.5f; %d down-regulated,%d up-regulated\n", genes[row].c_str(), lower, upper, cntl, cntu);
    }
    progress("Discretization rules are written to %s", stream_nm);
    fclose(fw);
  }

  void write_imported(const char* stream_nm) {
    int row, col;
    FILE *fw;
    fw = mustOpen(stream_nm, "w");
    fprintf(fw, "o");
    for (col = 0; col < cols; col++)
      fprintf(fw, "\t%s", conds[col].c_str());
    fputc('\n', fw);
    for (row = 0; row < rows; row++) {
      fprintf(fw, "%s", genes[row].c_str());
      for (col = 0; col < cols; col++)
        fprintf(fw, "\t%d", symbols[arr_c[row][col]]);
      fputc('\n', fw);
    }
    progress("Formatted data are written to %s", stream_nm);
    fclose(fw);
  }

  void read_list(FILE* fp) {
#define delims "\t\r\n"
    int i = 0, j = 0;
    sub_genes_row = 0;
    char line[MAXC];
    while (fgets(line, MAXC, fp) != NULL) {
      char *atom = strtok(line, delims);
      sub_genes[sub_genes_row] = atom;
      sub_genes_row++;
    }

    /*update the sub_list*/
    sublist.resize(rows);
    for (i = 0; i < rows; i++) sublist[i] = false;
    for (i = 0; i < sub_genes_row; i++)
      for (j = 0; j < rows; j++)
        if (strcmp(sub_genes[i].c_str(), genes[j].c_str()) == 0)
          sublist[j] = true;
  }

  void seed_update(const DiscreteArray& s) {
    for (int i = 0; i < cols; i++)
      profile[i][s[i]]++;
  }

  /* scan through all columns and identify the set within threshold,
   * "fuzziness" of the block is controlled by TOLERANCE (-c)
   */
  void scan_block(const std::vector<int> & gene_set, Block & b) {
    int i, j;
    int block_rows, cur_rows;
    block_rows = cur_rows = gene_set.size();

    int k;
    for (j = 0; j < cols; j++)
      for (k = 0; k < symbols.size(); k++)
        profile[j][k] = 0;
    for (j = 0; j < cur_rows; j++)
      seed_update(arr_c[gene_set[j]]);

    int btolerance = static_cast<int>(ceil(po.TOLERANCE* block_rows));
    for (j = 0; j < cols; j++) {
      /* See if this column satisfies tolerance */
      /* here i start from 1 because symbols[0]=0 */
      for (i = 1; i < symbols.size(); i++)	{
        if ((profile[j][i] >= btolerance)) {
          bool result = b.conds.insert(j).second;
          assert(result);
          break;
        }
      }
    }
  }

  /*************************************************************************/

  /* Identified clusters are backtraced to the original data, by
   * putting the clustered vectors together, identify common column
   */
  void print_bc(FILE* fw, const Block & b, const int & num) {
    int block_rows, block_cols;
    int num_1 = 0, num_2 = 0;
    /* block height (genes) */
    block_rows = b.block_rows();
    /* block_width (conditions) */
    block_cols = b.block_cols();
    fprintf(fw, "BC%03d\tS=%d\tPvalue:%LG \n", num, block_rows * block_cols, b.pvalue);
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
    /*printf ("BC%03d: #of 1 and -1 are:\t%d\t%d\n",num,num_1,num_2);
    fputc('\n', fw);*/
  }

  void update_colcand(std::vector<bool> & colcand, const DiscreteArray& g1, const DiscreteArray& g2) {
    int i;
    for (i = 0; i < cols; i++)
      if (colcand[i] && (g1[i] != g2[i]))
        colcand[i] = false;
  }

  /*calculate the weight of the edge with two vertices g1 and g2*/
  int intersect_row(const std::vector<bool> & colcand, const DiscreteArray& g1, const DiscreteArray& g2) {
    int i;
    int cnt = 0;
    for (i = 0; i < cols; i++)
      if (colcand[i] && (g1[i] == g2[i]) && (g1[i] != 0))
        cnt++;
    return cnt;
  }

  /*calculate the negative correlation between g1 and g2*/
  int reverse_row(const std::vector<bool> & colcand, const DiscreteArray& g1, const DiscreteArray& g2) {
    int i;
    int cnt = 0;
    for (i = 0; i < cols; i++) {
      if (colcand[i] && (symbols[g1[i]] == -symbols[g2[i]])) cnt++;
    }
    return cnt;
  }

  /* calculate the coverage of any row to the current consensus
  * cnt = # of valid consensus columns
  */
  int seed_current_modify(const DiscreteArray& s, std::vector<bool> &colcand, const int &components) {
    int i, k, flag, n;
    int threshold = static_cast <int> (ceil(components * po.TOLERANCE));
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
  bool check_seed(const Edge *e, const std::vector<Block> & bb) {
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
      b3 = MAX(bb[b1].block_cols(), bb[b2].block_cols());
      return !(e->score < b3/* (bb[b1]->block_cols + bb[b2]->block_cols) / 2*/);
    }
  }

  long double get_pvalue(const continuous & a, const int & b) {
    int i = 0;
    long double one = 1, pvalue = 0;
    long double poisson = one / exp(a);
    for (i = 0; i < b + 300; i++) {
      if (i > (b - 1)) pvalue = pvalue + poisson;
      else poisson = poisson * a / (i + 1);
    }
    return pvalue;
  }

  void block_init(Block & b,
    std::vector<int> & genes, std::vector<int> & scores,
    std::vector<bool> & candidates, const int & cand_threshold,
    int & components, std::vector<long double> & pvalues) {
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
    if (rows > 100)	{
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
      for (i = 0; i < rows; i++)
      {
        if (!candidates[i]) continue;
        if (po.IS_list && !sublist[i]) continue;
        cnt = intersect_row(colcand, arr_c[genes[0]], arr_c[i]);
        cnt_all += cnt;
        if (cnt < cand_threshold)
          candidates[i] = false;
        if (cnt > max_cnt)
        {
          max_cnt = cnt;
          max_i = i;
        }
      }
      cnt_ave = cnt_all / row_all;
      pvalue = get_pvalue(cnt_ave, max_cnt);
      if (po.IS_cond) {
        if (max_cnt < po.COL_WIDTH || max_i < 0 || max_cnt < b.cond_low_bound) break;
      }
      else {
        if (max_cnt < po.COL_WIDTH || max_i < 0) break;
      }
      if (po.IS_area)	score = components * max_cnt;
      else score = MIN(components, max_cnt);
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

  void print_params(FILE *fw) {
    std::string filedesc = "continuous";
    if (po.IS_DISCRETE)
      filedesc = "discrete";
    fprintf(fw, "# QUBIC version %.1f output\n", 1.9);
    fprintf(fw, "# Datafile %s: %s type\n", po.FN.c_str(), filedesc.c_str());
    fprintf(fw, "# Parameters: -k %d -f %.2f -c %.2f -o %d",
      po.COL_WIDTH, po.FILTER, po.TOLERANCE, po.RPT_BLOCK);
    if (!po.IS_DISCRETE)
      fprintf(fw, " -q %.2f -r %d", po.QUANTILE, po.DIVIDED);
    fprintf(fw, "\n\n");
  }

  /* compare function for qsort, descending by score */
  static bool block_cmpr(const Block & a, const Block & b) {
    return a.score > b.score;
  }

  /************************************************************************/
  int report_blocks(FILE* fw, std::vector<Block> bb) {
    int num = bb.size();

    print_params(fw);

    std::sort(bb.begin(), bb.end(), block_cmpr);

    int i, j, k;
    /*MIN MAX et al functions can be accessed in struct.h*/
    int n = MIN(num, po.RPT_BLOCK);
    bool flag;

    std::vector<Block> output;

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

        if (inter_rows * inter_cols > po.FILTER*cur_rows*cur_cols) {
          flag = false;
          break;
        }
        k++;
      }
      i++;
      if (flag) {
        print_bc(fw, b_ptr, j++);
        output.push_back(b_ptr);
      }
    }
    return j;
  }

  /************************************************************************/

  int cluster(FILE *fw, const std::vector<Edge *> & el) {
    std::vector<Block> bb;

    int j, k, components;

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
        else if ((po.IS_TFname) && (e->gene_one != TFindex) && (e->gene_two != TFindex)) flag = false;
        else if ((po.IS_list) && (!sublist[e->gene_one] || !sublist[e->gene_two])) flag = false;
      } else {
        flag = check_seed(e, bb);
        if ((po.IS_TFname) && (e->gene_one != TFindex) && (e->gene_two != TFindex)) flag = false;
        if ((po.IS_list) && (!sublist[e->gene_one] || !sublist[e->gene_two])) flag = false;
      }
      if (!flag) continue;

      for (j = 0; j < cols; j++)
        for (k = 0; k < symbols.size(); k++)
          profile[j][k] = 0;

      /*you must allocate a struct if you want to use the pointers related to it*/
      Block b;
      /*initial the b->score*/
      b.score = MIN(2, e->score);
      /*initial the b->pvalue*/
      b.pvalue = 1;

      /* initialize the stacks genes and scores */
      int ii;

      std::vector<int> genes_order, genes_reverse, scores;

      genes_order.reserve(rows);
      genes_reverse.reserve(rows);
      scores.reserve(rows);

      genes_order.push_back(e->gene_one);
      genes_order.push_back(e->gene_two);
      scores.push_back(1);
      scores.push_back(b.score);

      //for (ii = 2; ii < rows; ii++)
      //{
      //  dsPush(genes, -1);
      //  dsPush(scores, -1);
      //}

      /* branch-and-cut condition for seed expansion */
      int cand_threshold = static_cast<int>(floor(po.COL_WIDTH * po.TOLERANCE));
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
          if ((pvalues[k] == b.pvalue) && (k >= 2) && (dsItem(scores, k) != dsItem(scores, k + 1))) break;
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
      for (int ki = 0; ki < rows; ki++) {
        if (po.IS_list && !sublist[ki]) continue;
        m_cnt = intersect_row(colcand, arr_c[genes_order[0]], arr_c[ki]);
        if (candidates[ki] && (m_cnt >= floor(cnt* po.TOLERANCE))) {
          genes_order.push_back(ki);
          components++;
          candidates[ki] = false;
        }
      }

      /* add genes that negative regulated to the consensus */
      for (int ki = 0; ki < rows; ki++) {
        if (po.IS_list && !sublist[ki]) continue;
        m_cnt = reverse_row(colcand, arr_c[genes_order[0]], arr_c[ki]);
        if (candidates[ki] && (m_cnt >= floor(cnt * po.TOLERANCE))) {
          genes_reverse.push_back(ki);
          components++;
          candidates[ki] = false;
        }
      }

      /* store gene arrays inside block */

      scan_block(genes_order, b);
      if (b.block_cols() == 0) continue;
      //b.block_rows = components;
      if (po.IS_pvalue) b.score = -(100 * log(b.pvalue));
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
      verboseDot();
    }
    /* writes character to the current position in the standard output (stdout) and advances the internal file position indicator to the next position.
     * It is equivalent to putc(character,stdout).*/
    putchar('\n');
    return report_blocks(fw, bb);
  }

  /*make_graph subroutine prototypes */

  /* remove a row from the profile */
  void seed_deduct(const discrete *s) {
    for (size_t i = 0; i < cols; i++) {
      profile[i][s[i]]--;
    }
  }
  
  void make_graph(const char *fn, const DiscreteArrayList &arr_c, int &COL_WIDTH) {
    EdgeList EdgeList(arr_c, COL_WIDTH);

    FILE *fw = mustOpen(fn, "w");
    /* bi-clustering */
    progress("Clustering started");
    int n_blocks = 0;
#ifndef MAKE_GRAPH_ONLY
    n_blocks = cluster(fw, EdgeList.get_edge_list());
#endif
    printf("%d", EdgeList.get_edge_list().size());
    printf("%d clusters are written to %s\n", n_blocks, fn);
    /* clean up */
    fclose(fw);
  }

  /* expand subroutine prototypes */
  DiscreteArrayList another_arr_c;
  std::vector<std::string> another_genes;
  std::vector<std::string> another_conds;
  int another_rows;
  int another_cols;

  static int intersect_rowE(const std::vector<bool> & colcand, std::vector<discrete> & g1, std::vector<discrete> & g2, const int & cols) {
    int i, cnt = 0;
    for (i = 0; i < cols; i++)
      if (colcand[i] && (g1[i] == g2[i]) && g1[i] != 0)
        cnt++;
    return cnt;
  }

  int reverse_rowE(const std::vector<bool> & colcand, std::vector<discrete> & g1, std::vector<discrete> & g2, const int & cols) {
    int i, cnt = 0;
    for (i = 0; i < cols; i++)
      if (colcand[i] && (symbols[g1[i]] == -symbols[g2[i]]))
        cnt++;
    return cnt;
  }

  void init_expand() {
    another_genes = genes;
    another_conds = conds;
    another_arr_c = arr_c;
    another_rows = rows;
    another_cols = cols;
  }

  void make_charsets(const std::string &tfile, std::vector<discrete> &symbols) {
    discrete bb[USHRT_MAX];
    memset(bb, -1, USHRT_MAX*sizeof(*bb));
    charset_add(symbols, 0, bb);
    if (po.IS_DISCRETE) {
      for (size_t i = 0; i < arr.size(); i++)
        for (size_t j = 0; j < arr[0].size(); j++) {
          arr_c[i][j] = charset_add(symbols, (discrete)arr[i][j], bb);
        }
      printf("Discretized data contains %d classes with charset [ ", symbols.size());
      for (discrete i = 0; i < symbols.size(); i++)
        printf("%d ", symbols[i]);  printf("]\n");
    }
    else {
      for (size_t i = 0; i < arr.size(); i++)
        for (size_t j = 0; j < arr[0].size(); j++)
          arr_c[i][j] = 0;
      discretize((tfile + ".rules").c_str(), arr, bb, symbols, po.QUANTILE, arr_c, po.DIVIDED, genes);
    }
  }

  int run_qubic(const std::string &tfile, const double & rq, const double & rc, const double & rf, const int & rk, const discrete & rr, const int & ro, const int & rd) {
    int i = 0, j = 0;

    arr_c.resize(rows, DiscreteArray(cols));

    printf("\nQUBIC %s: greedy biclustering\n\n", VER);
    /* get the program options defined in get_options.c */
    /*set memory for the point which is declared in struct.h*/
    //AllocVar(po);
    /*Initialize the point*/
    po.FN = tfile;
    po.BN = " ";
    po.LN = " ";
    /* case 'l': strcpy(po.LN, optarg); po.IS_list =true; */
    po.TFname = " ";
    /* case 'T': strcpy(po.TFname, optarg); po.IS_TFname = true; */
    if (rd == 1) po.IS_DISCRETE = true;
    else po.IS_DISCRETE = false;
    // po.IS_DISCRETE = true;
    po.IS_TFname = false;
    po.IS_pvalue = false;
    /* case 'P': po.IS_pvalue = true; */
    po.COL_WIDTH = rk;
    po.DIVIDED = rr;
    po.QUANTILE = rq;
    po.TOLERANCE = rc;
    po.FP = NULL;
    po.FB = NULL;
    po.RPT_BLOCK = ro;
    po.SCH_BLOCK = 2 * po.RPT_BLOCK;
    /* ensure enough searching space */
    /*if (po.SCH_BLOCK < 1000) po.SCH_BLOCK = 1000;*/
    po.FILTER = rf;
    /* case 's': po.IS_SWITCH = true; */
    po.IS_area = false;
    /* case 'S': po.IS_area = true; */
    po.IS_cond = false;
    /* case 'C': po.IS_cond = true; */
    po.IS_list = false;

    if (po.IS_list)
      po.FL = mustOpen(po.LN.c_str(), "r");

    /* check if there exist a gene name equals to TFname by -T */
    for (int i = 0; i < genes.size(); i++)
      if (strcmp(genes[i].c_str(), po.TFname.c_str()) == 0)
        printf("%d\n", i);

    make_charsets(tfile, symbols);
    /*read in the sub-gene list*/
    if (po.IS_list) {
      sub_genes.resize(rows);
      read_list(po.FL);
    }
    /*we can do expansion by activate po.IS_SWITCH*/
  {
    /* the file that stores all blocks */
    if (po.IS_list)
      make_graph((tfile + ".block").c_str(), arr_c, po.COL_WIDTH);
    else
      make_graph((tfile + ".blocks").c_str(), arr_c, po.COL_WIDTH);
  }

  return 1;
  }

  int cgetbc(double *rbc, int *ro, char **filer, char **filew) {
    FILE *fpr = fopen(filer[0], "r");
    FILE *fpw = fopen((std::string(filew[0]) + ".bc").c_str(), "w");
    int i = 0, n = *ro, num = 0, ntemp = 0, nbc = 0;
    char temp[1000];
    char czero = '0', tempnum[20];
    while (fscanf(fpr, "%s", temp) != EOF) {
      if (strcmp(temp, "Genes") == 0) {
        fprintf(fpw, "BC%d\n", nbc);
        fscanf(fpr, "%s", tempnum);
        num = 0;
        for (i = 0; i < 20; i++) {
          if (tempnum[i] == ']')  break;
          else if (tempnum[i] == ':')  break;
          else if (tempnum[i] == '[')  continue;
          else {
            num *= 10;
            ntemp = tempnum[i] - czero;
            num += ntemp;
          }
        }
        /* fprintf(fpw, "%d\n", num); */
        for (i = 0; i < num; i++) {
          fscanf(fpr, "%s", temp);
          fprintf(fpw, "%s\t", temp);
        }
        fprintf(fpw, "\n");
        rbc[nbc + n] = num;
        fscanf(fpr, "%s", temp);
        if (strcmp(temp, "Conds") == 0) {
          /* fprintf(fpw, "%s\t", temp); */
          fscanf(fpr, "%s", tempnum);
          num = 0;
          for (i = 0; i < 20; i++) {
            if (tempnum[i] == ']')  break;
            else if (tempnum[i] == ':')  break;
            else if (tempnum[i] == '[')  continue;
            else {
              num *= 10;
              ntemp = tempnum[i] - czero;
              num += ntemp;
            }
          }
          /* fprintf(fpw, "%d\n", num); */
          for (i = 0; i < num; i++) {
            fscanf(fpr, "%s", temp);
            fprintf(fpw, "%s\t", temp);
          }
          fprintf(fpw, "\n");
          rbc[nbc + 2 * n] = num;
          rbc[nbc] = nbc;
          nbc++;
        }
      }
    }
    if (nbc < n) {
      nbc = n - nbc;
      for (i = 0; i < nbc; i++)
        rbc[nbc] = -1;
    }
    fclose(fpw);
    fclose(fpr);
    return 1;
  }

  public:
  void init_qubic(const std::string & tfile = "rQUBIC", const double & rq = 0.06, const double & rc = 0.95, const double & rf = 1, const int & rk = 2, const discrete & rr = 1, const int & ro = 100, const int & rd = 'F') {
    run_qubic(tfile, rq, rc, rf, rk, rr, ro, rd);

    //std::vector<std::vector<int> > data;
    //data.resize(rows, std::vector<int>(cols));

    /* formatted file */
    write_imported((tfile + ".chars").c_str());
    //for (int i = 0; i < rows; i++)
    //  for (int j = 0; j < cols; j++)
    //    data[i][j] = (int)symbols[arr_c[i][j]];
  }

  qubic(const std::vector<std::vector<float> > &data, const std::vector<std::string > &row_names, const std::vector<std::string > &col_names) {
    rows = row_names.size();
    cols = col_names.size();
    
    arr = data;
  
    genes = row_names;
    conds = col_names;
  }

  //qubic(float *r_data, const std::vector<std::string> & r_rowsnames, const std::vector<std::string> & r_colsnames, const int & r_rows, const int & r_cols) {
  //  cols = r_cols;
  //  rows = r_rows;
  //  arr.resize(rows, std::vector<continuous>(cols));
  //  for (size_t i = 0; i < rows; i++) {
  //    for (size_t j = 0; j < cols; j++)
  //      arr[i][j] = r_data[i + j*rows];
  //  }
  //  genes = r_rowsnames;
  //  conds = r_colsnames;
  //}
};

int r_main(const std::vector<std::vector<float> > &data, const std::vector<std::string > &row_names, const std::vector<std::string > &col_names, const std::string & tfile = "rQUBIC", const double & rq = 0.06, const double & rc = 0.95, const double & rf = 1, const int & rk = 2, const discrete & rr = 1, const int & ro = 100, const int & rd = 'F') {
  qubic qubic(data, row_names, col_names);
  qubic.init_qubic(tfile, rq, rc, rf, rk, rr, ro, rd);
  return 1;
}
