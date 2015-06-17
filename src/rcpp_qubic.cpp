#include <csignal>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <list>

#include <Rcpp.h>

#include "config.h"
#include "edge_list.h"
#include "fopen_matrix.h"
#include "matrix_float.h"
#include "option.h"
#include "qubic.h"

using namespace Rcpp;

template<typename T>
NumericMatrix from_vector(const std::vector<std::vector<T>> &result) {
  size_t nc = result.size();
  size_t nr = result[0].size();
  NumericMatrix m(nr, nc);
  for (size_t i = 0; i < nr; i++) {
    const std::vector<T> &result_i = result[i];
    if (result_i.size() != nc) stop("incompatible size");
    for (size_t j = 0; j < nc; j++) m(i, j) = result_i[j];
  }
  return m;
}

template<typename T, typename TMatrix>
std::vector<std::vector<T>> to_vector(const TMatrix &matrix) {
  auto nc = matrix.ncol();
  auto nr = matrix.nrow();
  std::vector<std::vector<T>> result(nr);
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) result[i].push_back(matrix(i, j));
  }
  return result;
}

NumericMatrix nothing(NumericMatrix matrix) {
  return from_vector<float>(to_vector<float, NumericMatrix>(matrix));
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

List from_blocks(const std::vector<Block> &r, const size_t nr, const size_t nc) {
  int number = r.size();
  auto x = LogicalMatrix(nr, number);
  auto y = LogicalMatrix(number, nc);
  for (int i = 0; i < number; i++) {
    for (auto it = r[i].genes_order.begin(); it != r[i].genes_order.end(); it++) x(*it, i) = true;
    for (auto it = r[i].genes_reverse.begin(); it != r[i].genes_reverse.end(); it++) x(*it, i) = true;
    for (auto it = r[i].conds.begin(); it != r[i].conds.end(); it++)
      y(i, *it) = true;
  }
  return List::create(
           Named("RowxNumber") = x,
           Named("NumberxCol") = y,
           Named("Number") = r.size(),
           Named("info") = get_list());
}

// [[Rcpp::export]]
List qubic(const NumericMatrix matrix, const short r, const double q,
           const double c, const int o, const double f, const int k, const bool P, const bool S, const bool C, const bool verbose) {
  // may treat abort() more friendly, see http://stackoverflow.com/a/3911102
  signal(SIGABRT,  &my_function_to_handle_aborts);
  try {    
    std::vector<Block> result = r_main_c(to_vector<float, NumericMatrix>(matrix), r, q, c, o, f, k, Option(P, S, C), verbose);
    return from_blocks(result, matrix.nrow(), matrix.ncol());
  } catch (double) {
    stop("catch");
  }
}

// [[Rcpp::export]]
List qubic_d(const IntegerMatrix matrix, const double c, const int o, const double f, const int k, const bool P, const bool S, const bool C, const bool verbose) {
  // may treat abort() more friendly, see http://stackoverflow.com/a/3911102
  signal(SIGABRT,  &my_function_to_handle_aborts);
  try {
    std::vector<Block> result = r_main_d(to_vector<short, IntegerMatrix>(matrix), c, o, f, k, Option(P, S, C), verbose);
    return from_blocks(result, matrix.nrow(), matrix.ncol());
  } catch (double) {
    stop("catch");
  }
}