#ifndef EDGE_LIST_H
#define EDGE_LIST_H

#include <vector>

#include "struct.h"

typedef std::vector<discrete> DiscreteArray;
typedef std::vector<DiscreteArray> DiscreteArrayList;


class EdgeList {
private:
  std::vector<Edge *> edge_list; int col_width;
  int get_key(const Edge* s);
public:
  const std::vector<Edge *> &get_edge_list() const;
  EdgeList(const DiscreteArrayList &, size_t&);
  ~EdgeList();
};

#endif
