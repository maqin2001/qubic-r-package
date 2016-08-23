#ifndef EDGE_LIST_H
#define EDGE_LIST_H

#include "discrete.h"
#include <cstddef> // size_t
#include <cassert>
#include <algorithm>

/* edge between two genes */
struct Edge {
  std::size_t gene_one;
  std::size_t gene_two;
  int score;
};

template <typename T>
class AdjMatrix {
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

inline unsigned str_intersect_r(const DiscreteArray &s1, const DiscreteArray &s2) {
  assert(s1.size() == s2.size());
  int common_cnt = 0;
  /* s1 and s2 of equal length, so we check s1 only */
  for (std::size_t i = 0; i < s1.size(); i++)
    if (s1[i] != 0 && s1[i] == s2[i]) // Changed order by zy26
      common_cnt++;
  return common_cnt;
}

class CountHelper {
  const DiscreteArrayList &arr_;
  const std::size_t col_width_;
public:
  std::size_t col_width() const {
    return col_width_;
  }
  virtual ~CountHelper() {
  }

  explicit CountHelper(const DiscreteArrayList& arr_c, std::size_t col_width) : arr_(arr_c), col_width_(col_width){ }

  std::size_t size() const {
    return arr_.size();
  }
  virtual int operator()(std::size_t i, std::size_t j) const {
    assert(i < j);
    return str_intersect_r(arr_[i], arr_[j]);
  }
};

class CountHelperRealTime : public CountHelper {
  explicit CountHelperRealTime(const DiscreteArrayList& arr_c, std::size_t col_width) : CountHelper(arr_c, col_width) { }
};

class CountHelperSaved : public CountHelper {
protected:
  std::vector<unsigned> intersects_;
public:
  virtual ~CountHelperSaved() {
  }

  explicit CountHelperSaved(const DiscreteArrayList& arr_c, std::size_t col_width) : CountHelper(arr_c, col_width), intersects_(arr_c.size() * (arr_c.size() - 1)) {
    for (std::size_t i = 0; i < arr_c.size(); i++)
      for (std::size_t j = i + 1; j < arr_c.size(); j++)
        intersects_[j * (j - 1) / 2 + i] = str_intersect_r(arr_c[i], arr_c[j]);
  }

  int operator()(std::size_t i, std::size_t j) const override {
    assert(i < j);
    return intersects_[j * (j - 1) / 2 + i];
  }
};

class CountHelperRanked : public CountHelperSaved {
  struct mycomparison {
    bool operator() (unsigned* lhs, unsigned* rhs) const {
      return *lhs < *rhs;
    }
  };
public:
  virtual ~CountHelperRanked() {
  }

  explicit CountHelperRanked(const DiscreteArrayList& arr_c, std::size_t col_width) : CountHelperSaved(arr_c, col_width) {
    std::vector<unsigned*> pintArray(intersects_.size());
    for (std::size_t i = 0; i < intersects_.size(); ++i) {
      pintArray[i] = &intersects_[i];
    }

    std::sort(pintArray.begin(), pintArray.end(), mycomparison());

    // Dereference the pointers and assign their sorted position. not deal tie
    for (std::size_t i = 0; i < intersects_.size(); ++i) {
      *pintArray[i] = i + 1;
    }
  }
};

class WeightedCountHelper : public CountHelperRanked {
  const std::vector<std::vector<float>>& weights_;
public:
  explicit WeightedCountHelper(const DiscreteArrayList& arr_c, const std::vector<std::vector<float>>& weights, std::size_t col_width) : CountHelperRanked(arr_c, col_width), weights_(weights) {}

  int operator()(std::size_t i, std::size_t j) const override {
    return CountHelperRanked::operator()(i, j) + weights_[i][j];
  }
};

class EdgeList {
  std::vector<Edge *> edge_list_;
public:
  const std::vector<Edge *> &get_seeds() const;
  EdgeList(const CountHelper& countHelper, bool verbose);
  ~EdgeList();
};

#endif
