#include "struct.h"

/**************************************************************************/
/* helper functions for error msgs for allocating memory */

void progress(char *format, ...)
/* Print progress message */
{
	va_list args;
	va_start(args, format);
	vfprintf(stdout, format, args);
	fprintf(stdout, "\n");
	va_end(args);
}

void verboseDot()
/* Print "i-am-alive" dot */
{
	putchar('.');
	fflush(stdout);
}

void errAbort(char *format, ...)
/* Print error message and exit */
{
	va_list args;
	va_start(args, format);
	fprintf(stderr, "[Error] ");
	vfprintf(stderr, format, args);
	fprintf(stderr, "\n");
	va_end(args);
	exit(1);
}

int count_intersect(const std::set<int> & ds1, const std::set<int> & ds2)
/* Return the number of common components between two arrays */
{
  int cnt = 0;
  std::set<int>::const_iterator first1 = ds1.begin();
  std::set<int>::const_iterator last1 = ds1.end();
  std::set<int>::const_iterator first2 = ds2.begin();
  std::set<int>::const_iterator last2 = ds2.end();
  while (first1 != last1 && first2 != last2) {
    if (*first1 < *first2) ++first1;
    else if (*first2 < *first1) ++first2;
    else {
      ++cnt; ++first1; ++first2;
    }
  }
  return cnt;
}

/**************************************************************************/
/* file-related operations */

FILE *mustOpen(const char *fileName, char *mode)
/* Open a file or die */
{
    FILE *f;

    if (sameString(fileName, "stdin")) return stdin;
    if (sameString(fileName, "stdout")) return stdout;
    if ((f = fopen(fileName, mode)) == NULL)
    {
        char *modeName = "";
        if (mode)
        {
            if (mode[0] == 'r') modeName = " to read";
            else if (mode[0] == 'w') modeName = " to write";
            else if (mode[0] == 'a') modeName = " to append";
        }
        errAbort("Can't open %s%s: %s", fileName, modeName, strerror(errno));
    }
    return f;
}

/**************************************************************************/
