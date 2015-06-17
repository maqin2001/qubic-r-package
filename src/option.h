#ifndef OPTION_H
#define OPTION_H

struct Option {
  Option(bool _p, bool _s, bool _c) : p(_p), s(_s), c(_c){}
  Option() : p(false), s(false), c(false) {}

  bool p;
  bool s;
  bool c;
};

struct DOption {
  double c;
  int o;
  double f;
  int k;
  Option option;
  bool verbose;
};

#endif
