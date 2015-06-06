#include "qubic.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <list>

#include "matrix_float.h"
#include "fopen_matrix.h"
#include "edge_list.h"
#include "config.h"

const char USAGE[] =
"===================================================================\n"
"[Usage]\n"
"$ ./qubic -i filename [argument list]\n"
"===================================================================\n"
"[Input]\n"
"-i : input file must be one of two tab-delimited formats\n"
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
//"-d : discrete data, where user should send their processed data\n"
//"     to different value classes, see above\n"
//"-b : the file to expand in specific environment\n"
//"-T : to-be-searched TF name, just consider the seeds containing current TF\n"
//"     default format: B1234\n"
//"-P : the flag to enlarge current biclsuter by the pvalue constrain\n"
//"-S : the flag using area as the value of bicluster to determine when stop\n"
//"-C : the flag using the lower bound of condition number (5 percents of the gene number)\n"
//"-l : the list of genes out of the input file on which we do bicluster\n"
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
//"-s : expansion flag\n"
"===================================================================\n";

void run_qubic(const std::vector<std::vector<float> > &data, const std::vector<std::string > &row_names, const std::vector<std::string > &col_names, const std::string & tfile = "rQUBIC", const double & rq = 0.06, const double & rc = 0.95, const double & rf = 1, const int & rk = 2, const short & rr = 1, const int & ro = 100, const int & rd = 'F') {
  r_main(data, row_names, col_names, tfile, rq, rc, rf, rk, rr, ro, rd);
}

void run_qubic(const MatrixFloat &matrix, const std::string & tfile = "rQUBIC", const double & rq = 0.06, const double & rc = 0.95, const double & rf = 1, const int & rk = 2, const short & rr = 1, const int & ro = 100, const int & rd = 'F') {
  run_qubic(matrix.get_data_const(), matrix.get_row_names(), matrix.get_col_names(), tfile, rq, rc, rf, rk, rr, ro, rd);
}

struct S { double d; int i; };
int get_key(const S &s) { return s.i - 1000; }

char *getCmdOption(char **begin, char **end, const std::string &option) {
  char ** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end)
    return *itr;
  return 0;
}

bool cmdOptionExists(char **begin, char **end, const std::string &option) {
  return std::find(begin, end, option) != end;
}

int main(int argc, char *argv[]) {
  if (cmdOptionExists(argv, argv + argc, "-h")) {
    printf(USAGE);
  }

  const char *file_name = cmdOptionExists(argv, argv + argc, "-i") ? getCmdOption(argv, argv + argc, "-i") : DEFAULT_FILENAME;
  double q = cmdOptionExists(argv, argv + argc, "-q") ? std::atof(getCmdOption(argv, argv + argc, "-q")) : 0.06;
  double c = cmdOptionExists(argv, argv + argc, "-c") ? std::atof(getCmdOption(argv, argv + argc, "-c")) : 0.95;
  double f = cmdOptionExists(argv, argv + argc, "-f") ? std::atof(getCmdOption(argv, argv + argc, "-f")) : 1.0;
  int k = cmdOptionExists(argv, argv + argc, "-k") ? std::atoi(getCmdOption(argv, argv + argc, "-k")) : 2;
  short r = cmdOptionExists(argv, argv + argc, "-r") ? std::atoi(getCmdOption(argv, argv + argc, "-r")) : 1;
  int o = cmdOptionExists(argv, argv + argc, "-o") ? std::atoi(getCmdOption(argv, argv + argc, "-o")) : 100;
  char d = cmdOptionExists(argv, argv + argc, "-d") ? 'T' : 'F';

  MatrixFloat matrix = FopenMatrix::load_matrix<float>(file_name);
  //printf("Size of matrix: %d", matrix.get_data_const().size());
#ifndef LOAD_FILE_ONLY
  run_qubic(matrix, file_name, q, c, f, k, r, o, d);
#endif
  return 0;
}
