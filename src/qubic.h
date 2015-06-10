#ifndef QUBIC_H
#define QUBIC_H 

#include <string>
#include <vector>

#include "block.h"

std::vector<Block> r_main(const std::vector<std::vector<float> > &data, const std::vector<std::string > &row_names, const std::vector<std::string > &col_names, const std::string &tfile, const double rq, const double rc, const double rf, const int rk, const short rr, const int ro, const bool rd);
std::vector<Block> r_main(const std::vector<std::vector<float> > &x, const short r = 1, const double q = 0.06, const double c = 0.95, const int o = 100, const double f = 1, const int rk = 2, const bool d = false);
std::vector<Block> r_main(const std::vector<std::vector<float> > &data);

#endif
