#include "qubic.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <list>

#include "fopen_matrix.h"
#include "edge_list.h"

const char USAGE[] =
"\n===================================================================\n"
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
"     -------------------------------------\n";

void run_qubic(const std::vector<std::vector<float> > &data, const std::vector<std::string > &row_names, const std::vector<std::string > &col_names) {
  r_main(data, row_names, col_names);
}

void run_qubic(const Matrix &matrix) {
  run_qubic(matrix.get_data_const(), matrix.get_row_names(), matrix.get_col_names());
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
  char default_file_name[] = "example";// arabidopsis - leaf";// ecoli_466_4297"; example";// 

  if (cmdOptionExists(argv, argv + argc, "-h")) {
    printf(USAGE);
  }

  char * file_name = getCmdOption(argv, argv + argc, "-file");

  if (!file_name) {
    file_name = default_file_name;
  }
  
  Matrix matrix = FopenMatrix::load_matrix(file_name);
  //printf("Size of matrix: %d", matrix.get_data_const().size());
#ifndef LOAD_FILE_ONLY
  run_qubic(matrix);
#endif
  return 0;
}
