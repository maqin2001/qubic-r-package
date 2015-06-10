#ifndef _STRUCT_H
#define _STRUCT_H

#include <cstdarg>
#include <set>

/* Two major data types */
typedef float continuous;
typedef short discrete;


/***** Structures *****/

/* edge between two genes */
typedef struct Edge {
  size_t gene_one;
  size_t gene_two;
  int score;
} Edge;

/* holds running options */
typedef struct Prog_options {
  bool IS_pvalue;
  bool IS_area;
  bool IS_cond;
  size_t SCH_BLOCK;
} Prog_options;

typedef unsigned short int bits16;

int count_intersect(const std::set<int> & ds1, const std::set<int> & ds2);

/* File-related operations */

#endif
