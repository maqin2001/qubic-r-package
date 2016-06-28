#include "../src/fopen_matrix.h"
#include <algorithm> // std::find
#include "qubic.h"

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
"-d : discrete data, where user should send their processed data\n"
"     to different value classes, see above\n"
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

char* getCmdOption(char** begin, char** end, const std::string& option) {
  char** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end) return *itr;
  return nullptr;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option) {
  return std::find(begin, end, option) != end;
}

double getCommand(int argc, char* argv[], const const char* option, double defaultValue) {
  return cmdOptionExists(argv, argv + argc, option) ? std::atof(getCmdOption(argv, argv + argc, option)) : defaultValue;
}

int getCommand(int argc, char* argv[], const const char* option, int defaultValue) {
  return cmdOptionExists(argv, argv + argc, option) ? std::atoi(getCmdOption(argv, argv + argc, option)) : defaultValue;
}

int main(int argc, char* argv[]) {
  if (cmdOptionExists(argv, argv + argc, "-h")) printf(USAGE);

  const char*      file_name = getCmdOption(argv, argv + argc, "-i");

  if (!file_name) {
    printf("Input file needed.");
    printf(USAGE);
    return-1;
  }

  short  r = static_cast<short>(getCommand(argc, argv, "-r", 1));
  double                    q = getCommand(argc, argv, "-q", 0.06);
  double                    c = getCommand(argc, argv, "-c", 0.95);
  int                       o = getCommand(argc, argv, "-o", 10);
  double               filter = getCommand(argc, argv, "-f", 1.0);
  int               col_width = getCommand(argc, argv, "-k", 2);

  bool                      d = cmdOptionExists(argv, argv + argc, "-d");

  const char*     w_file_name = getCmdOption(argv, argv + argc, "-w");
  const char*     b_file_name = getCmdOption(argv, argv + argc, "-b");

  if (d) {
    Matrix<short> matrix = FopenMatrix::load_matrix<short>(file_name);
    main_d(matrix.get_data(), matrix.get_row_names(), matrix.get_col_names(), file_name, c, o, filter, col_width, Option(), true);
  } else {
    Matrix<float> matrix = FopenMatrix::load_matrix<float>(file_name);
    main_c(matrix.get_data(), matrix.get_row_names(), matrix.get_col_names(), file_name, r, q, c, o, filter, col_width, Option(), true);
  }

  return 0;
}