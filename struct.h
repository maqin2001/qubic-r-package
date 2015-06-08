#ifndef _STRUCT_H
#define _STRUCT_H

#ifndef _GNU_SOURCE 
#define _GNU_SOURCE 
#endif

#include <cstdarg>
#include <cstdio> // FILE
#include <set>

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

/* holds running options */
typedef struct Prog_options{
	bool IS_pvalue;
	bool IS_area;
	bool IS_cond;
	int SCH_BLOCK;
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

int count_intersect(const std::set<int> & ds1, const std::set<int> & ds2);

/* File-related operations */

FILE *mustOpen(const char *fileName, const char *mode);
/* Open a file or die */

#endif
