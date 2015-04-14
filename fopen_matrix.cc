#include "fopen_matrix.h"

#include <cassert> 
#include <cstdio>
#include <cstring>

#include "matrix.h"

#define MAX_LINE 100000
#define LABEL_LEN 64 

Matrix FopenMatrix::load_matrix(const char* file_name, size_t reserved_count) {
  Matrix matrix(reserved_count);
  FILE *fp = fopen(file_name, "r");
  if (NULL == fp) {
    printf("Failed to open 'input.txt'");
    throw - 1;
  }

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

    std::vector<float> line_data;
    line_data.reserve(matrix.col_names.size());

    for (int i = 0; i < matrix.col_names.size(); i++) {
      float tmp;
      int items_read = fscanf(fp, "%f", &tmp);
      line_data.push_back(tmp);
    }
    assert(matrix.col_names.size() == line_data.size());
    matrix.data.push_back(line_data);
  }

  assert(matrix.row_names.size() == matrix.data.size());
  
  return matrix; // RVO
}
