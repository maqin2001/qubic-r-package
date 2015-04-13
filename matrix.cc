#include "matrix.h"

Matrix::Matrix(size_t reserved_count) { row_names.reserve(reserved_count); col_names.reserve(reserved_count); data.reserve(reserved_count); }
