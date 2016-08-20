#include "edge_list.h"
#include "block.h"
std::vector<Block> cluster(const DiscreteArrayListWithSymbols& all, const std::vector<Edge *>& el,
  std::size_t COL_WIDTH, double TOLERANCE, bool IS_cond, bool IS_area,
  bool IS_pvalue, std::size_t SCH_BLOCK, int RPT_BLOCK, double FILTER, double f, bool verbose);
