#include "run_qubic.h"
#include "qubic.h"

std::vector<Block> run_qubic(const std::vector<std::vector<float> > &data, const std::vector<std::string > &row_names,
  const std::vector<std::string > &col_names, const std::string &tfile = "rQUBIC", const double &rq = 0.06,
  const double  &rc = 0.95, const double &rf = 1, const int &rk = 2, const short &rr = 1, const int &ro = 100,
  const bool &rd = false) {
  return r_main(data, row_names, col_names, tfile, rq, rc, rf, rk, rr, ro, rd);
}

std::vector<Block> run_qubic(const MatrixFloat &matrix, const std::string &tfile, const double &q, const double &c,
  const double &f, const int &k, const short &r, const int &o, const bool &d) {
  return run_qubic(matrix.get_data(), matrix.get_row_names(), matrix.get_col_names(), tfile, q, c, f, k, r, o, d);
}