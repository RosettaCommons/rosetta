#ifndef FEM_UTILS_STRING_TO_DOUBLE_HPP
#define FEM_UTILS_STRING_TO_DOUBLE_HPP

#include <fem/utils/char.hpp>
#include <fem/utils/simple_streams.hpp>
#include <tbxx/optional_copy.hpp>

namespace fem { namespace utils {

  template <size_t Size>
  struct string_to_double_ad_hoc_limts;

  // Simple upper estimates, not guarding against floating-point errors,
  // but guarding against unusual sizeof(double).
  template <>
  struct string_to_double_ad_hoc_limts<8>
  {
    static const unsigned max_digit_count = 16;
    static const unsigned max_exp_int = 308;
  };

  // Ideally in string_to_double_ad_hoc_limts<8>, but it does not
  // seem to be straightforward to do.
  static const double one_e_minus_0_16[] = {
    1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10,
    1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-16};
  static const double one_e_two_pow_0_8[] = {
    1e1, 1e2, 1e4, 1e8, 1e16, 1e32, 1e64, 1e128, 1e256};

  struct string_to_double
  {
    double result;
    tbxx::optional_copy<std::string> error_message;
    bool stream_end;

    string_to_double() {}

    string_to_double(
      simple_istream& inp,
      unsigned decimal_point_substitute=0,
      int exp_scale_substitute=0)
    {
      reset();
      convert(inp, decimal_point_substitute, exp_scale_substitute);
    }

    void
    reset()
    {
      result = 0;
      error_message.release();
      stream_end = false;
    }

    void
    set_error_message(
      int c)
    {
      static const std::string inp_err(
        "Input error while reading floating-point value.");
      static const std::string inp_eoi(
        "End of input while reading floating-point value.");
      static const char* err_inv =
        "Invalid character while reading floating-point value: ";
      if (is_stream_err(c)) {
        error_message = inp_err;
      }
      else if (is_stream_end(c)) {
        error_message = inp_eoi;
        stream_end = true;
      }
      else {
        error_message = err_inv + format_char_for_display(c);
      }
    }

    void
    convert(
      simple_istream& inp,
      unsigned decimal_point_substitute=0, // "d" part of Fortran "w.d" format
      int exp_scale_substitute=0) // Fortran "P" scaling
    {
      typedef string_to_double_ad_hoc_limts<sizeof(double)> ad_hoc_limits;
      static const std::string err_oor(
        "Out-of-range error while reading floating-point value.");
      int c = inp.get();
      while (is_whitespace(c)) {
        c = inp.get();
      }
      bool sign_flag;
      if (c == '-') {
        sign_flag = true;
        c = inp.get();
      }
      else {
        sign_flag = false;
        if (c == '+') {
          c = inp.get();
        }
      }
      bool had_dot = false;
      bool had_digit = false;
      bool had_19 = false;
      int rexp = -1;
      int dexp = 0;
      if (c == '.') {
        had_dot = true;
        c = inp.get();
      }
      while (is_digit(c)) {
        had_digit = true;
        int d = digit_as_int(c);
        if (d != 0) {
          had_19 = true;
          if (dexp <= ad_hoc_limits::max_digit_count) {
            result += d * one_e_minus_0_16[dexp];
          }
        }
        if (had_19) {
          if (!had_dot) rexp++;
          dexp++;
        }
        else if (had_dot) {
          rexp--;
        }
        c = inp.get();
        if (c == '.') {
          if (had_dot) {
            break;
          }
          had_dot = true;
          c = inp.get();
        }
      }
      if (!had_digit) {
        set_error_message(c);
        return;
      }
      if (sign_flag) result *= -1;
      if (!had_dot && decimal_point_substitute > 0) {
        rexp -= decimal_point_substitute;
      }
      int exp_int = 0;
      if (c == 'e' || c == 'E' || c == 'd' || c == 'D') {
        c = inp.get();
        bool exp_sign_flag;
        if (c == '-') {
          exp_sign_flag = true;
          c = inp.get();
        }
        else {
          exp_sign_flag = false;
          if (c == '+') {
            c = inp.get();
          }
        }
        if (!is_digit(c)) {
          result = 0;
          set_error_message(c);
          return;
        }
        int exp_int_limit = rexp;
        if (exp_int_limit < 0) exp_int_limit *= -1;
        exp_int_limit += ad_hoc_limits::max_exp_int;
        exp_int = digit_as_int(c);
        while (true) {
          c = inp.get();
          if (!is_digit(c)) {
            break;
          }
          exp_int *= 10;
          exp_int += digit_as_int(c);
          if (exp_int > exp_int_limit) {
            result = 0;
            error_message = err_oor;
            return;
          }
        }
        if (exp_sign_flag) exp_int *= -1;
      }
      else {
        exp_int = -exp_scale_substitute;
      }
      if (!is_stream_end_or_err(c)) {
        inp.backup();
      }
      exp_int += rexp;
      bool exp_sign_flag = (exp_int < 0);
      if (exp_sign_flag) exp_int *= -1;
      // using ideas found in ruby-1.8.6.tar.gz, missing/strtod.c
      double exp_double = 1;
      for(unsigned i=0; exp_int != 0; i++) {
        if (exp_int & 1) exp_double *= one_e_two_pow_0_8[i];
        exp_int >>= 1;
      }
      if (exp_sign_flag) result /= exp_double;
      else               result *= exp_double;
    }
  };

}} // namespace fem::utils

#endif // GUARD
