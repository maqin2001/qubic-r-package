#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <string>

class Matrix {
public:
  const std::vector<std::vector<float> >& get_data_const() const { return data; }
  const std::vector<std::string>& get_row_names() const { return row_names; }
  const std::vector<std::string>& get_col_names() const { return col_names; }
public: // protected: TODO: Fix me
  std::vector<std::string> row_names;
  std::vector<std::string> col_names;
  std::vector<std::vector<float> > data;
  Matrix(size_t reserved_count);
};

#endif
