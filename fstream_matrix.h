#ifndef FSTREAM_MATRIX_H
#define FSTREAM_MATRIX_H

#include "matrix.h"

namespace FstreamMatrix {
  Matrix load_matrix(const char* file_name, size_t reserved_count = 4096);
}

#endif
