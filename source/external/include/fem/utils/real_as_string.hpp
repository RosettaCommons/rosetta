#ifndef FEM_UTILS_REAL_AS_STRING_HPP
#define FEM_UTILS_REAL_AS_STRING_HPP

#include <cstdio>

namespace fem { namespace utils {

  struct float_as_string_list_directed
  {
    char begin[64];
    int n;

    float_as_string_list_directed(
      double const& val) // intentionally forcing promotion to double
    {
#if defined(_MSC_VER)
      char const* cfmte = "%15.7E";
      char const* cfmt = cfmte;
#else
      char const* cfmt = "%14.7E";
#endif
      float av = (val < 0 ? -val : val);
      if (av < 1e3) {
        if (av < 1e1) {
          if      (av < 1e-1) /*pass*/;
          else if (av < 1e0) cfmt = "%10.7f    ";
          else               cfmt = "%10.6f    ";
        }
        else if (av < 1e2) cfmt = "%10.5f    ";
        else               cfmt = "%10.4f    ";
      }
      else if (av < 1e5) {
        if (av < 1e4) cfmt = "%10.3f    ";
        else          cfmt = "%10.2f    ";
      }
      else if (av < 1e7) {
        if (av < 1e6) cfmt = "%10.1f    ";
        else          cfmt = "%9.0f.    ";
      }
      n = std::sprintf(begin, cfmt, val);
#if defined(_MSC_VER)
      if (cfmt == cfmte && n == 15 && begin[12] == '0') {
        begin[12] = begin[13];
        begin[13] = begin[14];
        begin[14] = '\0';
        n = 14;
      }
#endif
    }
  };

  struct double_as_string_list_directed
  {
    char buffer[64];
    char* begin;
    int n;

    double_as_string_list_directed(
      double const& val)
    :
      begin(buffer)
    {
      char const* cfmte = "%23.15E";
      char const* cfmt = cfmte;
      double av = (val < 0 ? -val : val);
      if (av < 1e7) {
        if (av < 1e3) {
          if (av < 1e1) {
            if      (av < 1e-1) /*pass*/;
            else if (av < 1e0) cfmt = "%18.15f     ";
            else               cfmt = "%18.14f     ";
          }
          else if (av < 1e2) cfmt = "%18.13f     ";
          else               cfmt = "%18.12f     ";
        }
        else if (av < 1e5) {
          if (av < 1e4) cfmt = "%18.11f     ";
          else          cfmt = "%18.10f     ";
        }
        else if (av < 1e6) cfmt = "%18.9f     ";
        else               cfmt = "%18.8f     ";
      }
      else if (av < 1e11) {
        if (av < 1e9) {
          if (av < 1e8) cfmt = "%18.7f     ";
          else          cfmt = "%18.6f     ";
        }
        else if (av < 1e10) cfmt = "%18.5f     ";
        else                cfmt = "%18.4f     ";
      }
      else if (av < 1e13) {
        if (av < 1e12) cfmt = "%18.3f     ";
        else           cfmt = "%18.2f     ";
      }
      else if (av < 1e14) cfmt = "%18.1f     ";
      else if (av < 1e15) cfmt = "%17.0f.     ";
      n = std::sprintf(buffer, cfmt, val);
      if (n == 23 && cfmt == cfmte) {
        char c = buffer[20];
        if ((c == '+' || c == '-') && buffer[0] == ' ') {
          begin++;
          buffer[24] = '\0';
          buffer[23] = buffer[22];
          buffer[22] = buffer[21];
          buffer[21] = '0';
        }
      }
    }
  };

}} // namespace fem::utils

#endif // GUARD
