#include "catch.hpp"

#include "../src/fopen_matrix.h"
#include "../src/qubic.h"
#include "../src/discretize.h"

#define EXAMPLE_FILENAME "testdata/example"
#define TOY_EXAMPLE_FILENAME "testdata/toy_example"

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