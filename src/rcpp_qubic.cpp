#include <csignal>
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
  size_t nc = result.size();
  size_t nr = result[0].size();
  NumericMatrix m(nr, nc);
  for (size_t i = 0; i < nr; i++) {
    const std::vector<T>& result_i = result[i];
    if (result_i.size() != nc) stop("incompatible size");
    for (size_t j = 0; j < nc; j++) {
      m(i, j) = result_i[j];
    }
  }
  return m;
}

template<typename T>
std::vector<std::vector<T>> to_vector(const NumericMatrix& matrix) {
  auto nc = matrix.ncol();
  auto nr = matrix.nrow();
  std::vector<std::vector<T>> result(nr);
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) {
      result[i].push_back(matrix(i, j));
    }
  }
  return result;
}

NumericMatrix nothing(NumericMatrix matrix) {
  return from_vector<float>(to_vector<float>(matrix));
}

List get_list() {
  return List::create();
}

extern "C" void my_function_to_handle_aborts(int signal_number) {
  /*Your code goes here. You can output debugging info.
  If you return from this function, and it was called
  because abort() was called, your program will exit or crash anyway
  (with a dialog box on Windows).
  */
  stop("abort()");
}

List from_blocks(const std::vector<Block> &result) {  
    int number = result.size();

    auto x = LogicalMatrix(matrix.nrow(), number);
    auto y = LogicalMatrix(number, matrix.ncol());

    for (int i = 0; i < number; i++) {
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

// [[Rcpp::export]]
List qubic(const NumericMatrix matrix, const short r, const double q, const double c, const int o, const double f) {
  signal(SIGABRT, &my_function_to_handle_aborts); // may treat abort() more friendly, see http://stackoverflow.com/a/3911102
  try {
    std::vector<Block> result = r_main(to_vector<float>(matrix), r, q, c, o, f);
    return from_blocks(result);
  } catch (double) {
    stop("catch");
  }
}
