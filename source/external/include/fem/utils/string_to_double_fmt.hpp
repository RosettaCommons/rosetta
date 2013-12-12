#ifndef FEM_UTILS_STRING_TO_DOUBLE_FMT_HPP
#define FEM_UTILS_STRING_TO_DOUBLE_FMT_HPP

#include <fem/utils/string_to_double.hpp>

namespace fem { namespace utils {

  struct string_to_double_fmt : string_to_double
  {
    string_to_double_fmt() {}

    string_to_double_fmt(
      simple_istream& inp,
      unsigned w,
      unsigned d,
      bool blanks_zero,
      int exp_scale)
    {
      reset();
      if (w == 0) return;
      std::string buf;
      buf.reserve(w);
      for(unsigned i=0;i<w;i++) {
        int c = inp.get();
        if (is_stream_err(c)) {
          set_error_message(c);
          return;
        }
        if (is_stream_end(c)) {
          if (i == 0) {
            set_error_message(c);
            return;
          }
          break;
        }
        if (is_end_of_line(c)) {
          inp.backup();
          break;
        }
        if (c == ' ') {
          if (blanks_zero) buf += '0';
        }
        else {
          buf += c;
        }
      }
      if (buf.size() == 0) return;
      simple_istream_from_std_string buf_inp(buf);
      convert(buf_inp, d, exp_scale);
    }
  };

}} // namespace fem::utils

#endif // GUARD
