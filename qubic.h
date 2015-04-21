#ifndef QUBIC_H
#define QUBIC_H 

#include <string>
#include <vector>

int r_main(float *r_data, const std::vector<std::string> & r_rowsnames, const std::vector<std::string> & r_colsnames, const int & r_rows, const int & r_cols, const std::string & tfile, const float & rq, const double & rc, const double & rf, const int & rk, const short & rr, const int & ro, const int & rd);

int r_main(const std::vector<std::vector<float> > &data, const std::vector<std::string > &row_names, const std::vector<std::string > &col_names, const std::string & tfile, const double & rq, const double & rc, const double & rf, const int & rk, const short & rr, const int & ro, const int & rd);
#endif
