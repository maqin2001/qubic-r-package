#ifndef EDGE_LIST_H
#define EDGE_LIST_H

#include "discrete.h"
#include <cstddef> // size_t
#include <cassert>

/* edge between two genes */
struct Edge {
  std::size_t gene_one;
  std::size_t gene_two;
  int score;
};

template <typename T>
class AdjMatrix
{
  std::vector<T> matrix_;
  unsigned int size_;

public:
  AdjMatrix(unsigned int size) : size_(size) {
    matrix_.resize(size_* (size_ + 1) / 2);
  }

  std::size_t get_index(unsigned x, unsigned y) const {
    return (y * (2 * size_ - y + 1)) / 2 + (x - y);
  }

  T& operator()(unsigned int x, unsigned int y) {
    assert(x < size_ && y < size_);
    if (x >= y) return matrix_[get_index(x, y)];
    return matrix_[get_index(y, x)];
  }
  const T& operator()(unsigned int x, unsigned int y) const {
    assert(x < size_ && y < size_);
    if (x >= y) return matrix_[get_index(x, y)];
    return matrix_[get_index(y, x)];
  }
};

class CountHelper {
  static int str_intersect_r(const DiscreteArray &s1, const DiscreteArray &s2) {
    assert(s1.size() == s2.size());
    int common_cnt = 0;
    /* s1 and s2 of equal length, so we check s1 only */
    for (std::size_t i = 0; i < s1.size(); i++)
      if ((s1[i] != 0) && (s1[i] == s2[i])) // Changed order by zy26
        common_cnt++;
    return common_cnt;
  }
  const DiscreteArrayList &arr_;
public:
  virtual ~CountHelper()
  {
  }

  explicit CountHelper(const DiscreteArrayList& arr_c) : arr_(arr_c) {}
  std::size_t size() const {
    return arr_.size();
  }
  virtual int operator()(std::size_t i, std::size_t j) const {
    return str_intersect_r(arr_[i], arr_[j]);
  }
};

class WeightedCountHelperUnused : public CountHelper {
  const AdjMatrix<double>& weights_;
public:
  explicit WeightedCountHelperUnused(const DiscreteArrayList& arr_c, const AdjMatrix<double>& weights) : CountHelper(arr_c), weights_(weights) {}

  int operator()(std::size_t i, std::size_t j) const override {
    return CountHelper::operator()(i, j) * weights_(i, j);
  }
};

class WeightedCountHelper : public CountHelper {
  const std::vector<std::vector<float>>& weights_;
public:
  explicit WeightedCountHelper(const DiscreteArrayList& arr_c, const std::vector<std::vector<float>>& weights) : CountHelper(arr_c), weights_(weights) {}

  int operator()(std::size_t i, std::size_t j) const override {
    return CountHelper::operator()(i, j) + weights_[i][j];
  }
};

class EdgeList {
  std::vector<Edge *> edge_list_;
public:
  const std::vector<Edge *> &get_edge_list() const;
  EdgeList(std::size_t&, const CountHelper& countHelper, bool verbose);
  ~EdgeList();
};

#endif
