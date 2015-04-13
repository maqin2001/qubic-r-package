#ifndef FOPEN_MATRIX_H
#define FOPEN_MATRIX_H

#include "matrix.h"

class FopenMatrix : public Matrix {
public:
  static FopenMatrix load_matrix(const char* file_name, size_t reserved_count = 4096);
private:
  FopenMatrix(size_t reserved_count);
};

#endif