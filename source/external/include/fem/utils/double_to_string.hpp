#ifndef FEM_UTILS_DOUBLE_TO_STRING_HPP
#define FEM_UTILS_DOUBLE_TO_STRING_HPP

#include <fem/utils/char.hpp>
#include <fem/utils/int_types.hpp>

namespace fem { namespace utils {

  struct double_to_string_scientific_notation
  {
    static const int buffer_capacity = 256;
    char buffer[buffer_capacity];

    double_to_string_scientific_notation() {}

    double_to_string_scientific_notation(
      double value,
      int w, // "w" part of Fortran "w.d" format
      int d, // "d" part of Fortran "w.d" format
      int exp_scale, // Fortran "P" scaling
      char e_or_d)
    {
      TBXX_ASSERT(w > 0);
      TBXX_ASSERT(w < buffer_capacity);
      TBXX_ASSERT(d >= 0);
      int nd = d;
      int na = d;
      if (exp_scale == 0) {
        // pass
      }
      else if (exp_scale > 0) {
        nd++;
        na -= exp_scale - 1;
      }
      else if (exp_scale < 0) {
        nd += exp_scale;
      }
      if (value == 0) {
        if (w < 1 + d + 4 || d == 0 || nd <= 0) {
          std::memset(buffer, '*', w);
          return;
        }
        int i_buf = w;
        buffer[--i_buf] = '0';
        buffer[--i_buf] = '0';
        buffer[--i_buf] = '+';
        buffer[--i_buf] = e_or_d;
        for(int i=0;i<d;i++) {
          buffer[--i_buf] = '0';
        }
        buffer[--i_buf] = '.';
        if (i_buf != 0) {
          buffer[--i_buf] = '0';
          while (i_buf != 0) {
            buffer[--i_buf] = ' ';
          }
        }
        return;
      }
      //
      struct pow_tab_entry {
        int i;
        double g;
        double l;
      };
      static const pow_tab_entry pow_tab[] = {
        {256, 1e256, 1e-256},
        {128, 1e128, 1e-128},
        { 64, 1e064, 1e-064},
        { 32, 1e032, 1e-032},
        { 16, 1e016, 1e-016},
        {  8, 1e008, 1e-008},
        {  4, 1e004, 1e-004},
        {  2, 1e002, 1e-002},
        {  1, 1e001, 1e-001}};
      bool v_negative = (value < 0);
      if (v_negative) value = -value;
      int required_w = (v_negative ? 2 : 1) + nd + 4;
      if (w < required_w || d == 0 || nd <= 0 || na < 0) {
        std::memset(buffer, '*', w);
        return;
      }
      int iexp = 0;
      if (value > 1) {
        for(int i_t=0;i_t<9;i_t++) {
          pow_tab_entry const& t = pow_tab[i_t];
          if (value > t.g) {
            iexp += t.i;
            value *= t.l;
          }
        }
        iexp++;
        value *= 0.1;
      }
      else if (value < 1) {
        for(int i_t=0;i_t<9;i_t++) {
          pow_tab_entry const& t = pow_tab[i_t];
          if (value < t.l) {
            iexp -= t.i;
            value *= t.g;
          }
        }
      }
      static const int max_nd_significant = 16;
      static const double multipliers_dbl[] = {
        1e0,  1e1,  1e2,  1e3,  1e4,  1e5,  1e6,  1e7,  1e8, 1e9,
        1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16};
      typedef utils::int64_t i64;
      TBXX_ASSERT(sizeof(i64) >= sizeof(double));
      static i64 multipliers_i64[max_nd_significant+1] = {0};
      if (multipliers_i64[0] == 0) {
        multipliers_i64[0] = 1;
        for(int i=1;i<=max_nd_significant;i++) {
          multipliers_i64[i] = multipliers_i64[i-1] * 10;
        }
      }
      int nd_significant = std::min(nd, max_nd_significant);
      int nd_unknown = nd - nd_significant;
      double multiplier = multipliers_dbl[nd_significant];
      i64 ival = static_cast<i64>(value * multiplier + 0.5);
      if (ival == multipliers_i64[nd_significant]) {
        ival = multipliers_i64[nd_significant-1];
        iexp++;
      }
      iexp -= exp_scale;
      //
      int i_buf = w;
      bool e_negative = (iexp < 0);
      if (e_negative) iexp = -iexp;
      while (iexp != 0) {
        buffer[--i_buf] = int_as_digit(iexp % 10);
        iexp /= 10;
      }
      while (i_buf > w-2) {
        buffer[--i_buf] = '0';
      }
      buffer[--i_buf] = (e_negative ? '-' : '+');
      if (i_buf != w-4) buffer[--i_buf] = e_or_d;
      //
      for(int i=0;i<na;i++) {
        if (nd_unknown != 0) {
          buffer[--i_buf] = '0';
          nd_unknown--;
        }
        else {
          buffer[--i_buf] = int_as_digit(ival % 10);
          ival /= 10;
        }
      }
      buffer[--i_buf] = '.';
      if (ival == 0) {
        if (i_buf != 0) {
          buffer[--i_buf] = '0';
        }
      }
      else {
        do {
          if (nd_unknown != 0) {
            buffer[--i_buf] = '0';
            nd_unknown--;
          }
          else {
            buffer[--i_buf] = int_as_digit(ival % 10);
            ival /= 10;
          }
        }
        while (ival != 0);
      }
      if (i_buf == 0) {
        if (v_negative) buffer[0] = '-';
      }
      else {
        buffer[--i_buf] = (v_negative ? '-' : ' ');
        while (i_buf != 0) {
          buffer[--i_buf] = ' ';
        }
      }
    }
  };

}} // namespace fem::utils

#endif // GUARD
