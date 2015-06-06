#ifndef _STRUCT_H
#define _STRUCT_H

#ifndef _GNU_SOURCE 
#define _GNU_SOURCE 
#endif

#include <string>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <cmath>
#include <cctype>
#include <cstdarg>
#include <climits>
#include <ctime>

#include "fib.h"

#include <vector>
#include <set>
#include <stack>

/***** Useful macros *****/

/* Compatibility of __attribute__ with non-GNU */
#ifndef __GNUC__
#  define __attribute__(x) /* Nothing */
#endif

#ifndef __cplusplus
/* Pretend that C has boolean type */
#define TRUE 1
#define FALSE 0

#define boolean unsigned char

#ifndef bool
#define bool unsigned char
#endif
#endif

#define MAX(a,b)  ((a)>(b)?(a):(b))
#define MIN(a,b)  ((a)<(b)?(a):(b))

/* Strings */
/* strcmp: a zero value indicates that both strings are equal.
 * a value greater than zero indicates that the first character that does not match has a greater value in str1 than in str2; 
 * And a value less than zero indicates the opposite.
 */
#define sameString(a, b) (strcmp((a), (b))==0)
/* Returns TRUE if two strings are same */

/* Debugging */
#define internalErr() errAbort("Internal error %s %d", __FILE__, __LINE__)

/* Constants */
//#define LABEL_LEN 64 

#ifndef NULL
#define NULL 0
#endif

/* Two major data types */
typedef float continuous;
typedef short discrete;


/***** Structures *****/

/* edge between two genes */
typedef struct Edge{
  size_t gene_one;
  size_t gene_two;
	int score;
} Edge;

/* biclustering block */
typedef struct Block{
  std::set<int> genes_order;
  std::set<int> genes_reverse;
  std::set<int> conds;
	int score;
  const size_t block_rows() const { return genes_order.size() + genes_reverse.size(); }
  const size_t block_cols() const { return conds.size(); }
  const bool contains(int gene) const { 
    return (genes_order.find(gene) != genes_order.end()) || (genes_reverse.find(gene) != genes_reverse.end()); }
	int cond_low_bound;
	double significance;
	long double pvalue;
} Block;

/* holds running options */
typedef struct Prog_options{
	std::string FN;
  std::string BN;
  std::string LN;
	bool IS_DISCRETE;
	bool IS_TFname;
	bool IS_pvalue;
	bool IS_area;
	bool IS_cond;
	bool IS_list;
  size_t COL_WIDTH;
	discrete DIVIDED;
	int SCH_BLOCK;
	int RPT_BLOCK;
	double FILTER;
  double QUANTILE;
	double TOLERANCE;
  std::string TFname;
	FILE* FP;
	FILE* FB;
	FILE* FL;
} Prog_options;

typedef unsigned short int bits16;

/***** Helper functions *****/

void progress(const char *format, ...)
/* Print progress message */
     __attribute__((format(printf, 1, 2)));

void verboseDot();
/* Print "i-am-alive" dot */

void errAbort(const char *format, ...)
/* Print error message to stderr and exit */
     __attribute__((noreturn, format(printf, 1, 2)));

void uglyTime(char *label, ...);
/* Print label and how long it's been since last call.  Call with 
 * a NULL label to initialize. */

#define dsClear(pds) ((pds).resize(0))
/* Remove all the data in the stack */

#define dsSize(pds) ((pds).size())
/* Return the size of the stack */

#define dsItem(pds, j) ((pds)[j])
/* Return the j-th item in the stack */

int count_intersect(const std::set<int> & ds1, const std::set<int> & ds2);

/* File-related operations */

FILE *mustOpen(const char *fileName, const char *mode);
/* Open a file or die */

#endif
