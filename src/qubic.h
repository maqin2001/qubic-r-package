#ifndef QUBIC_H
#define QUBIC_H

#include <string>
#include <vector>

#include "block.h"

std::vector<Block> main_c(const std::vector<std::vector<float>> &x, const std::vector<std::string> &row_names,
  const std::vector<std::string> &col_names, const std::string &tfile, const short r,
  const double q, const double c, const int o, const double f, const int k);

std::vector<Block> main_d(const std::vector<std::vector<short>> &x, const std::vector<std::string> &row_names,
  const std::vector<std::string> &col_names, const std::string &tfile, const double c, const int o, const double f,
  const int k);

std::vector<Block> r_main_c(const std::vector<std::vector<float>> &x, const short r = 1, const double q = 0.06,
  const double c = 0.95, const int o = 100, const double f = 1, const int k = 2);

std::vector<Block> r_main_d(const std::vector<std::vector<short>> &x,
  const double c = 0.95, const int o = 100, const double f = 1, const int k = 2);

#endif
