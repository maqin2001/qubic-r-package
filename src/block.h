#ifndef BLOCK_H
#define BLOCK_H

#include <set>
#include <cstddef> // size_t

/* biclustering block */
class Block {
public:
  std::set<int> genes_order;
  std::set<int> genes_reverse;
  std::set<int> conds;
  int score;
  const size_t block_rows() const { return genes_order.size() + genes_reverse.size(); }
  const size_t block_cols() const { return conds.size(); }
  const bool contains(int gene) const {
    return (genes_order.find(gene) != genes_order.end()) || (genes_reverse.find(gene) != genes_reverse.end());
  }
  int cond_low_bound;
  double significance;
  long double pvalue;
};

#endif
