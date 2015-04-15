#ifndef FOPEN_MATRIX_H
#define FOPEN_MATRIX_H

#include "matrix.h"

namespace FopenMatrix {
  Matrix load_matrix(const char* file_name, size_t reserved_count = 4096);
}

#endif
