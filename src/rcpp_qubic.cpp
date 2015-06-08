#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <list>

#include <Rcpp.h>

#include "qubic.h"
#include "matrix_float.h"
#include "fopen_matrix.h"
#include "edge_list.h"
#include "config.h"

using namespace Rcpp;

template<typename T>
NumericMatrix from_vector(const std::vector<std::vector<T>>& result) {
  int nc = result.size();
  int nr = result[0].size();
  NumericMatrix m(nr, nc);
  for (int j = 0; j < nc; j++) {
    const std::vector<T>& result_j = result[j];
    if (result_j.size() != nr) stop("incompatible size");
    for (int i = 0; i < nr; i++) {
      m(i,j) = result_j[i];
    }
  }
  return m;
}

template<typename T>
std::vector<std::vector<T>> to_vector(const NumericMatrix& matrix) {
  auto nc = matrix.ncol();
  auto nr = matrix.nrow();
  std::vector<std::vector<T>> result(nc);
  for (int j = 0; j < nc; j++) {
    for (int i = 0; i < nr; i++) {
      result[j].push_back(matrix(i,j));
    }
  }
  return result;
}

// [[Rcpp::export]]
NumericMatrix nothing(NumericMatrix matrix) {
  return from_vector<float>(to_vector<float>(matrix));
}

List get_list() {
  return List::create();
}

// [[Rcpp::export]]
List qubic(NumericMatrix matrix) {
  std::vector<Block> result = r_main(to_vector<float>(matrix));
  int number = result.size();
  
  auto x = LogicalMatrix(matrix.nrow(), number);
  auto y = LogicalMatrix(number, matrix.ncol());
  
  for(int i = 0; i < number; i++) {
    for (auto it = result[i].genes_order.begin(); it != result[i].genes_order.end(); it++) {
      x(*it, i) = true;
    }
    for (auto it = result[i].genes_reverse.begin(); it != result[i].genes_reverse.end(); it++) {
      x(*it, i) = true;
    }
    for (auto it = result[i].conds.begin(); it != result[i].conds.end(); it++) {
      y(i, *it) = true;
    }
  }

  return List::create(
    Named("RowxNumber") = x,
    Named("NumberxCol") = y,
    Named("Number") = result.size(),
    Named("info") = get_list());
}
