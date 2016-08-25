#include "catch.hpp"

#include "../src/fopen_matrix.h"
#include "../src/qubic.h"
#include "../src/discretize.h"

#include <ctime>

#define EXAMPLE_FILENAME "testdata/example"
#define TOY_EXAMPLE_FILENAME "testdata/toy_example"
#define ECOLI_FILENAME "testdata/ecoli_466_4297"

TEST_CASE("run_qubic_test") {
  auto matrix = FopenMatrix::load_matrix<float>(EXAMPLE_FILENAME);

  REQUIRE(matrix.row_names.size() == 9);
  REQUIRE(matrix.col_names.size() == 264);

  auto data = matrix.get_data();

  REQUIRE(data.size() == matrix.row_names.size());
  REQUIRE(data[0].size() == matrix.col_names.size());
  REQUIRE(data[data.size() - 1].size() == matrix.col_names.size());

  std::vector<Block> result = main_c(matrix.get_data(), matrix.get_row_names(), matrix.get_col_names(), EXAMPLE_FILENAME, 1, 0.06, 0.95, 10, 1, 2, Option(), true);
  REQUIRE(result.size() == 2);
  REQUIRE(result[0].genes_order.size() == 3);

  DiscreteArrayList vectors = discretize(matrix.get_data(), 1, 0.06);
  std::vector<Block> result2 = r_main(vectors, 0.95, 10, 1, 2, Option(), true);
  REQUIRE(result2.size() == 2);
  REQUIRE(result2[0].genes_order.size() == 3);
}

TEST_CASE("run_qubic_toy_test") {
  auto matrix = FopenMatrix::load_matrix<short>(TOY_EXAMPLE_FILENAME);

  REQUIRE(matrix.row_names.size() == 20);
  REQUIRE(matrix.col_names.size() == 12);

  auto data = matrix.get_data();

  REQUIRE(data.size() == matrix.row_names.size());
  REQUIRE(data[0].size() == matrix.col_names.size());
  REQUIRE(data[data.size() - 1].size() == matrix.col_names.size());

  std::vector<Block> result = main_d(matrix.get_data(), matrix.get_row_names(), matrix.get_col_names(), TOY_EXAMPLE_FILENAME, 0.95, 10, 1, 2, Option(), true);
  REQUIRE(result.size() == 4);
  REQUIRE(result[0].genes_order.size() == 8);

  std::vector<Block> result2 = r_main(matrix.get_data(), 0.95, 100, 1, 2, Option(), true);
  REQUIRE(result2.size() == 4);
  REQUIRE(result2[0].genes_order.size() == 8);
}

TEST_CASE("run_ecoli_test") {
  auto matrix = FopenMatrix::load_matrix<float>(ECOLI_FILENAME);

  REQUIRE(matrix.row_names.size() == 4297);
  REQUIRE(matrix.col_names.size() == 466);

  auto data = matrix.get_data();

  REQUIRE(data.size() == matrix.row_names.size());
  REQUIRE(data[0].size() == matrix.col_names.size());
  REQUIRE(data[data.size() - 1].size() == matrix.col_names.size());

  std::clock_t begin = std::clock();
  DiscreteArrayList vectors = discretize(matrix.get_data(), 1, 0.06);
  std::vector<Block> result2 = r_main(vectors, 0.95, 100, 1.0, 466 / 20, Option(), true);
  std::clock_t end = std::clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

  printf("elapsed_secs: %f\n", elapsed_secs);

  if (elapsed_secs <= 16) WARN("Less than 16 sec.");
  else if (elapsed_secs >= 17) WARN("More than 17 sec.");

  REQUIRE(result2.size() == 100);
  REQUIRE(result2[0].genes_order.size() == 131);
  REQUIRE(result2[0].block_rows() == 437);
  REQUIRE(result2[0].block_cols() == 29);
}