#ifndef QUBIC_H
#define QUBIC_H 

#include <string>
#include <vector>

#include "block.h"

std::vector<Block> r_main(const std::vector<std::vector<float> > &data, const std::vector<std::string > &row_names, const std::vector<std::string > &col_names, const std::string & tfile, const double & rq, const double & rc, const double & rf, const int & rk, const short & rr, const int & ro, const bool &rd);
//std::vector<Block> r_main(const std::vector<std::vector<float> > &data, const double & rq = 0.06, const double & rc = 0.95, const double & rf = 1, const int & rk = 2, const short & rr = 1, const int & ro = 100, const bool &rd = false);
std::vector<Block> r_main(const std::vector<std::vector<float> > &data);

#endif
