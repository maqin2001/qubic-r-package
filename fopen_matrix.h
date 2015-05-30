#ifndef FOPEN_MATRIX_H
#define FOPEN_MATRIX_H

#include "matrix.h"

namespace FopenMatrix {
  template<typename T>
  Matrix<T> load_matrix(const char* file_name, size_t reserved_count = 4096);

  namespace internal {
#include <cassert> 
#include <cstdio>
#include <cstring>

#define MAX_LINE 100000
#define LABEL_LEN 64 

    template<typename T>
    Matrix<T> load_matrix_from_file(FILE* fp, size_t reserved_count) {
      Matrix<T> matrix(reserved_count); // RVO

      char line[MAX_LINE];
      if (fgets(line, MAX_LINE, fp) == NULL)  throw - 1;

      const char delims[] = " \t\r\n";
      char *pch = strtok(line, delims);
      if (pch == NULL) throw - 1;
      pch = strtok(NULL, delims); // ignore the first value
      while (pch != NULL) {
        matrix.col_names.push_back(pch);
        pch = strtok(NULL, delims);
      }

      char value[LABEL_LEN];
      while (1 == fscanf(fp, "%s", value)) {
        matrix.row_names.push_back(value);

        std::vector<T> line_data(matrix.col_names.size());

        for (size_t i = 0; i < matrix.col_names.size(); i++) {
          if (std::is_floating_point<T>::value)
            fscanf(fp, "%f", &line_data[i]);
          else if (std::is_integral<T>::value)
            fscanf(fp, "%d", &line_data[i]);
          else
            throw - 1;
        }
        matrix.data.push_back(line_data);
      }
      assert(matrix.row_names.size() == matrix.data.size());
      return matrix; // RVO
    }
  }

  template<typename T>
  Matrix<T> FopenMatrix::load_matrix(const char* file_name, size_t reserved_count) {
    FILE *fp = fopen(file_name, "r");
    if (NULL == fp) {
      printf("Failed to open 'input.txt'");
      throw - 1;
    }
    Matrix<T> matrix = internal::load_matrix_from_file<float>(fp, reserved_count); // RVO
    fclose(fp);
    return matrix; // RVO
  }
}
#endif
