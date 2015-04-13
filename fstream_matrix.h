#ifndef FSTREAM_MATRIX_H
#define FSTREAM_MATRIX_H

#include "matrix.h"

class FstreamMatrix : public Matrix {
public:
  static FstreamMatrix load_matrix(const char* file_name, size_t reserved_count = 4096);
private:
  FstreamMatrix(size_t reserved_count);
};

#endif