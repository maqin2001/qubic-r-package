#include "run_qubic.h"
#include "qubic.h"

std::vector<Block> run_qubic_c(const MatrixFloat &matrix, const std::string &tfile /*= "rQUBIC"*/,
  const short &r /*= 1*/, const double &q /*= 0.06*/, const double &c /*= 0.95*/, const int &o /*= 100*/,
  const double &f /*= 1*/, const int &k /*= 2*/) {
  return main_c(matrix.get_data(), matrix.get_row_names(), matrix.get_col_names(), tfile, r, q, c, o, f, k);
}

std::vector<Block> run_qubic_d(const Matrix<short> &matrix, const std::string &tfile /*= "rQUBIC"*/,
  const double &c /*= 0.95*/, const int &o /*= 100*/,
  const double &f /*= 1*/, const int &k /*= 2*/) {
  return main_d(matrix.get_data(), matrix.get_row_names(), matrix.get_col_names(), tfile, c, o, f, k);
}
