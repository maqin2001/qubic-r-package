
#ifndef DISCRETIZE_H
#define DISCRETIZE_H

#include <vector>

#include "discrete.h"
#include "struct.h"

struct rule {
  float lower;
  float upper;
  size_t cntl;
  size_t cntu;
};

std::vector<std::vector<discrete>> discretize(const std::vector<std::vector<continuous>> &arr, const double f,
  const discrete divided, std::vector<rule> &genes_rules);

#endif