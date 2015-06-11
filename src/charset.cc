#include "charset.h"
discrete charset_add(std::vector<discrete> &ar, const discrete &s, discrete *bb) {
  /*A signed short can hold all the values between SHRT_MIN  and SHRT_MAX inclusive.SHRT_MIN is required to be -32767 or less,SHRT_MAX must be at least 32767*/
  int ps = s + SHRT_MAX;
  if (bb[ps] < 0) {
    bb[ps] = static_cast<discrete>(ar.size());
    ar.push_back(s);
  }
  return bb[ps];
}

DiscreteArrayListWithSymbols make_charsets_d(const std::vector<std::vector<discrete>> &arr) {
  DiscreteArrayListWithSymbols all;
  all.list.resize(arr.size(), DiscreteArray(arr[0].size()));
  discrete bb[USHRT_MAX];
  memset(bb, -1, USHRT_MAX * sizeof(*bb));
  charset_add(all.symbols, 0, bb);
  for (size_t i = 0; i < arr.size(); i++)
    for (size_t j = 0; j < arr[0].size(); j++)
      all.list[i][j] = charset_add(all.symbols, (discrete)arr[i][j], bb);
  fprintf(stdout, "Discretized data contains %d classes with charset [ ", static_cast<unsigned int>(all.symbols.size()));
  for (size_t i = 0; i < all.symbols.size(); i++)
    fprintf(stdout, "%d ", all.symbols[i]);
  fprintf(stdout, "]\n");
  return all;
}