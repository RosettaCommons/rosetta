#ifndef FEM_READ_HPP
#define FEM_READ_HPP

#include <fem/common.hpp>
#include <fem/format.hpp>
#include <fem/star.hpp>
#include <fem/str_arr_ref.hpp>
#include <fem/utils/misc.hpp>
#include <fem/utils/string_to_double_fmt.hpp>

#define IOSTAT_OK     0
#define IOSTAT_ERROR  1
#define IOSTAT_END   -1

namespace fem {

  class read_loop // TODO copy-constructor potential performance problem
  {
    private:
      utils::slick_ptr<utils::simple_istream> inp;
      bool first_inp_get;
      format::token_loop fmt_loop;
      bool blanks_zero;
      int exp_scale;
      io_modes io_mode;
      int *iostat_ptr;

    public:

      read_loop(
        common& cmn,
        int const& unit,
        unformatted_type const&)
      :
        inp(cmn.io.simple_istream(unit)),
        first_inp_get(true),
        blanks_zero(false),
        exp_scale(0),
        io_mode(io_unformatted),
        iostat_ptr(0)
      {}

      read_loop(
        common& cmn,
        int const& unit,
        star_type const&)
      :
        inp(cmn.io.simple_istream(unit)),
        first_inp_get(true),
        blanks_zero(false),
        exp_scale(0),
        io_mode(io_list_directed),
        iostat_ptr(0)
      {}

      read_loop(
        common& cmn,
        int const& unit,
        str_cref fmt)
      :
        inp(cmn.io.simple_istream(unit)),
        first_inp_get(true),
        fmt_loop(fmt),
        blanks_zero(false),
        exp_scale(0),
        io_mode(io_formatted),
        iostat_ptr(0)
      {}

      read_loop(
        str_cref const& internal_file,
        star_type const&)
      :
        inp(utils::slick_ptr<utils::simple_istream>(new
          utils::simple_istream_from_char_ptr_and_size(
            internal_file.elems(), internal_file.len()))),
        first_inp_get(true),
        blanks_zero(false),
        exp_scale(0),
        io_mode(io_list_directed),
        iostat_ptr(0)
      {}

      read_loop(
        str_cref const& internal_file,
        str_cref fmt)
      :
        inp(utils::slick_ptr<utils::simple_istream>(new
          utils::simple_istream_from_char_ptr_and_size(
            internal_file.elems(), internal_file.len()))),
        first_inp_get(true),
        fmt_loop(fmt),
        blanks_zero(false),
        exp_scale(0),
        io_mode(io_formatted),
        iostat_ptr(0)
      {}

      read_loop&
      rec(
        int const&)
      {
        inp.reset();
        throw TBXX_NOT_IMPLEMENTED();
      }

      read_loop&
      iostat(
        int& iostat)
      {
          this -> iostat_ptr = &iostat;
          iostat = IOSTAT_OK;
          return *this;
      }

      std::string const&
      next_edit_descriptor()
      {
        while (true) {
          utils::token const* t = fmt_loop.next_executable_token();
          std::string const& tv = t->value;
          if (t->type == "string") {
            inp.reset();
            throw TBXX_NOT_IMPLEMENTED();
          }
          else if (t->type == "op") {
            if (tv[0] == ':') {
              // ignored
            }
            else if (tv[0] == '/') {
              skip_to_end_of_line();
              int c = inp_get();
              if (utils::is_stream_end(c)) {
                inp.reset();
                if(this -> iostat_ptr != 0) *iostat_ptr = IOSTAT_END;
                throw read_end("End of input while reading string");
              }
              inp -> backup();
            }
            else if (tv[0] == '$') {
              inp.reset();
              throw TBXX_NOT_IMPLEMENTED();
            }
            else {
              inp.reset();
              throw TBXX_UNREACHABLE_ERROR();
            }
          }
          else if (t->type == "format") {
            if (utils::ends_with_char(tv, 'x')) {
              if (tv.size() == 1) {
                process_fmt_x(1);
              }
              else {
                process_fmt_x(utils::signed_integer_value(
                  tv.data(), 0, tv.size()-1));
              }
            }
            else if (std::strchr("adefgilz", tv[0]) != 0) {
              return tv;
            }
            else if (utils::ends_with_char(tv, 'p')) {
              if (tv.size() == 1) {
                exp_scale = 1;
              }
              else {
                exp_scale = utils::signed_integer_value(
                  tv.data(), 0, tv.size()-1);
              }
            }
            else if (tv[0] == 't') {
              inp.reset();
              throw TBXX_NOT_IMPLEMENTED();
            }
            else if (tv[0] == 's') {
              inp.reset();
              throw TBXX_NOT_IMPLEMENTED();
            }
            else if (tv[0] == 'b') {
              blanks_zero = (tv[1] == 'z');
            }
            else {
              inp.reset();
              throw TBXX_UNREACHABLE_ERROR();
            }
          }
          else {
            inp.reset();
            throw TBXX_UNREACHABLE_ERROR();
          }
        }
      }

      int
      inp_get()
      {
        int result = inp->get();
        if (utils::is_stream_err(result)) {
          inp.reset();
          if(this -> iostat_ptr != 0) *iostat_ptr = IOSTAT_ERROR;
          throw io_err("Error during read");
        }
        if (first_inp_get || io_mode == io_unformatted) {
          first_inp_get = false;
          if (utils::is_stream_end(result)) {
            inp.reset();
            if(this -> iostat_ptr != 0) *iostat_ptr = IOSTAT_END;
            throw read_end("End of input during read");
          }
        }
        if (io_mode == io_formatted && result == '\r') {
          int next_char = inp->get();
          if (next_char == '\n') {
            result = '\n';
          }
          else {
            inp->backup();
          }
        }
        return result;
      }

      void
      process_fmt_x(
        unsigned n)
      {
        for(unsigned i=0;i<n;i++) {
          int c = inp_get();
          if (c == utils::stream_end) {
            return;
          }
          if (utils::is_end_of_line(c)) {
            inp->backup();
            return;
          }
        }
      }

      read_loop&
      operator,(
        char& val)
      {
        inp.reset();
        throw TBXX_NOT_IMPLEMENTED();
        return *this;
      }

      read_loop&
      operator,(
        bool& val)
      {
        if (io_mode == io_unformatted) {
          from_stream_unformatted(
            reinterpret_cast<char*>(&val),
            sizeof(bool));
        }
        else if (io_mode == io_list_directed) {
          inp.reset();
          throw TBXX_NOT_IMPLEMENTED();
        }
        else {
          inp.reset();
          throw TBXX_NOT_IMPLEMENTED();
        }
        return *this;
      }

      read_loop&
      operator,(
        integer_star_1& val)
      {
        if (io_mode == io_unformatted) {
          inp.reset();
          throw TBXX_NOT_IMPLEMENTED();
        }
        else if (io_mode == io_list_directed) {
          inp.reset();
          throw TBXX_NOT_IMPLEMENTED();
        }
        else {
          inp.reset();
          throw TBXX_NOT_IMPLEMENTED();
        }
        return *this;
      }

      read_loop&
      operator,(
        integer_star_2& val)
      {
        if (io_mode == io_unformatted) {
          from_stream_unformatted(
            reinterpret_cast<char*>(&val),
            sizeof(integer_star_2));
        }
        else if (io_mode == io_list_directed) {
          inp.reset();
          throw TBXX_NOT_IMPLEMENTED();
        }
        else {
          inp.reset();
          throw TBXX_NOT_IMPLEMENTED();
        }
        return *this;
      }

      read_loop&
      operator,(
        integer_star_4& val)
      {
        if (io_mode == io_unformatted) {
          from_stream_unformatted(
            reinterpret_cast<char*>(&val),
            sizeof(integer_star_4));
        }
        else if (io_mode == io_list_directed) {
          val = static_cast<int>(read_star_long());
        }
        else {
          std::string const& ed = next_edit_descriptor();
          int n = ed.size();
          if (ed[0] == 'i' && n > 1) {
            n = utils::unsigned_integer_value(ed.data(), 1, n);
            val = static_cast<int>(read_fmt_long(n));
          }
          else {
            val = static_cast<int>(read_star_long());
          }
        }
        return *this;
      }

      read_loop&
      operator,(
        integer_star_8& val)
      {
        if (io_mode == io_unformatted) {
          from_stream_unformatted(
            reinterpret_cast<char*>(&val),
            sizeof(integer_star_8));
        }
        else if (io_mode == io_list_directed) {
          inp.reset();
          throw TBXX_NOT_IMPLEMENTED();
        }
        else {
          inp.reset();
          throw TBXX_NOT_IMPLEMENTED();
        }
        return *this;
      }

      read_loop&
      operator,(
        float& val)
      {
        if (io_mode == io_unformatted) {
          from_stream_unformatted(
            reinterpret_cast<char*>(&val),
            sizeof(float));
        }
        else {
          val = static_cast<float>(
            (io_mode == io_formatted
              ? read_fmt_double()
              : read_star_double()));
        }
        return *this;
      }

      read_loop&
      operator,(
        double& val)
      {
        if (io_mode == io_unformatted) {
          from_stream_unformatted(
            reinterpret_cast<char*>(&val),
            sizeof(double));
        }
        else {
          val = (io_mode == io_formatted
            ? read_fmt_double()
            : read_star_double());
        }
        return *this;
      }

      read_loop&
      operator,(
        std::complex<float>& val)
      {
        if (io_mode == io_unformatted) {
          float re, im;
          from_stream_unformatted(
            reinterpret_cast<char*>(&re),
            sizeof(float));
          from_stream_unformatted(
            reinterpret_cast<char*>(&im),
            sizeof(float));
          val = std::complex<float>(re, im);
        }
        else if (io_mode == io_list_directed) {
          inp.reset();
          throw TBXX_NOT_IMPLEMENTED();
        }
        else {
          inp.reset();
          throw TBXX_NOT_IMPLEMENTED();
        }
        return *this;
      }

      read_loop&
      operator,(
        std::complex<double>& val)
      {
        if (io_mode == io_unformatted) {
          double re, im;
          from_stream_unformatted(
            reinterpret_cast<char*>(&re),
            sizeof(double));
          from_stream_unformatted(
            reinterpret_cast<char*>(&im),
            sizeof(double));
          val = std::complex<double>(re, im);
        }
        else if (io_mode == io_list_directed) {
          inp.reset();
          throw TBXX_NOT_IMPLEMENTED();
        }
        else {
          inp.reset();
          throw TBXX_NOT_IMPLEMENTED();
        }
        return *this;
      }

      read_loop&
      operator,(
        str_ref const& val)
      {
        if (io_mode == io_unformatted) {
          from_stream_unformatted(val.elems(), val.len());
        }
        else {
          int vl = val.len();
          int n = vl;
          if (io_mode == io_formatted) {
            std::string const& ed = next_edit_descriptor();
            if (ed[0] == 'a' && ed.size() > 1) {
              n = utils::unsigned_integer_value(ed.data(), 1, ed.size());
            }
          }
          for(int i=0;i<n-vl;i++) {
            int c = inp_get();
            if (utils::is_stream_end(c)) {
              if (i == 0) {
                inp.reset();
                if(this -> iostat_ptr != 0) *iostat_ptr = IOSTAT_END;
                throw read_end("End of input while reading string");
              }
              break;
            }
            if (utils::is_end_of_line(c)) {
              inp->backup();
            }
          }
          for(int i=0;i<vl;i++) {
            int c = inp_get();
            if (utils::is_stream_end(c)) {
              if (i == 0) {
                inp.reset();
                if(this -> iostat_ptr != 0) *iostat_ptr = IOSTAT_END;
                throw read_end("End of input while reading string");
              }
              val[i] = ' ';
            }
            else if (utils::is_end_of_line(c)) {
              inp->backup();
              val[i] = ' ';
            }
            else
            {
                val[i] = c;
            }
          }
        }
        return *this;
      }

      template <typename T, size_t Ndims>
      read_loop&
      operator,(
        arr_ref<T, Ndims> const& val)
      {
        T* v = val.begin();
        size_t n = val.size_1d();
        for(size_t i=0;i<n;i++) {
          (*this), v[i];
        }
        return *this;
      }

      template <size_t Ndims>
      read_loop&
      operator,(
        str_arr_ref<Ndims> const& val)
      {
        size_t n = val.size_1d();
        int l = val.len();
        char* val_begin = val.begin();
        for(size_t i=0;i<n;i++) {
          (*this), str_ref(&val_begin[i*l], l);
        }
        return *this;
      }

      void
      skip_to_end_of_line()
      {
        while (true) {
          int c = inp_get();
          if (   utils::is_stream_end(c)
              || utils::is_end_of_line(c)) {
            break;
          }
        }
      }

      ~read_loop()
      {
        if (inp.get() == 0) return;
        if (io_mode == io_unformatted) {
          skip_to_end_of_unformatted_record();
        }
        else {
          skip_to_end_of_line();
        }
      }

      long
      read_fmt_long(
        unsigned n)
      {
        bool had_non_blank = false;
        bool negative = false;
        long result = 0;
        for(unsigned i=0;i<n;i++) {
          int c = inp_get();
          if (utils::is_stream_end(c)) {
            break;
          }
          if (c == ',') {
            break;
          }
          if (utils::is_end_of_line(c)) {
            inp->backup();
            break;
          }
          if (c != ' ') {
            if (!had_non_blank) {
              had_non_blank = true;
              if (c == '-') {
                negative = true;
                continue;
              }
              if (c == '+') {
                continue;
              }
            }
            if (!utils::is_digit(c)) {
              inp.reset();
              if(this -> iostat_ptr != 0) *iostat_ptr = IOSTAT_ERROR;
              throw io_err(
                "Invalid character while reading integer value.");
            }
            result *= 10;
            result += utils::digit_as_int(c);
          }
        }
        if (negative) result *= -1;
        return result;
      }

      long
      read_star_long()
      {
        while (true) { // loop scanning for first non-whitespace
          int c = inp_get();
          if (utils::is_stream_end(c)) {
            inp.reset();
            if(this -> iostat_ptr != 0) *iostat_ptr = IOSTAT_END;
            throw read_end(
              "End of input while reading integer value");
          }
          if (!utils::is_whitespace(c)) {
            bool negative = (c == '-');
            if (negative || c == '+') {
              c = inp_get();
              if (utils::is_stream_end(c)) {
                throw read_end(
                  "End of input while reading integer value");
              }
            }
            long result = 0;
            while (true) { // loop collecting digits
              if (!utils::is_digit(c)) {
                throw io_err(
                  "Invalid character while reading integer value.");
              }
              result *= 10;
              result += utils::digit_as_int(c);
              c = inp_get();
              if (   utils::is_stream_end(c)
                  || utils::is_whitespace(c)
                  || c == ',') {
                if (negative) result *= -1;
                if (utils::is_end_of_line(c)) inp->backup();
                return result;
              }
            }
          }
        }
      }

      void
      throw_if_conv_error_message(
        utils::string_to_double const& conv)
      {
        if (conv.error_message) {
          inp.reset();
          if (conv.stream_end) {
                if(this -> iostat_ptr != 0) *iostat_ptr = IOSTAT_END;
            throw read_end(*conv.error_message);
          }
          if(this -> iostat_ptr != 0) *iostat_ptr = IOSTAT_ERROR;
          throw io_err(*conv.error_message);
        }
      }

      double
      read_fmt_double()
      {
        std::string const& ed = next_edit_descriptor();
        int n = ed.size();
        if (n < 2 || std::strchr("defg", ed[0]) == 0) {
          return read_star_double();
        }
        int iw = utils::unsigned_integer_scan(ed.data(), 1, ed.size());
        int w = utils::unsigned_integer_value(ed.data(), 1, iw);
        int d = 0;
        if (iw+1 != ed.size()) {
          d = utils::unsigned_integer_value(ed.data(), iw+1, ed.size());
        }
        utils::string_to_double_fmt conv(*inp, w, d, blanks_zero, exp_scale);
        throw_if_conv_error_message(conv);
        first_inp_get = false;
        return conv.result;
      }

      double
      read_star_double()
      {
        utils::string_to_double conv(*inp);
        throw_if_conv_error_message(conv);
        first_inp_get = false;
        int c = inp_get();
        if (   utils::is_stream_end(c)
            || utils::is_whitespace(c)
            || c == ',') {
          if (utils::is_end_of_line(c)) inp->backup();
          return conv.result;
        }
        inp.reset();
        if(this -> iostat_ptr != 0) *iostat_ptr = IOSTAT_ERROR;
        throw io_err(
          "Invalid character while reading floating-point value: "
          + utils::format_char_for_display(c));
      }

      void
      from_stream_unformatted(
        char* target,
        unsigned target_size)
      {
        for(unsigned i=0;i<target_size;i++) {
          int ic = inp_get();
          char c = static_cast<char>(ic);
          if (c == end_of_unformatted_record) {
            if (inp_get() != ic) {
              inp.reset();
              if(this -> iostat_ptr != 0) *iostat_ptr = IOSTAT_END;
              throw read_end("End of record during unformatted read");
            }
          }
          target[i] = c;
        }
      }

      void
      skip_to_end_of_unformatted_record()
      {
        while (true) {
          char c = static_cast<char>(inp_get());
          if (c == end_of_unformatted_record) {
            if (inp_get() == 0) {
              break;
            }
          }
        }
      }
  };

  struct read_from_string : public read_loop
  {
    read_from_string(
      std::string const& internal_file,
      str_cref fmt)
    :
      read_loop(
        str_cref(
          internal_file.data(),
          static_cast<int>(internal_file.size())),
        fmt)
    {}
  };

  struct common_read
  {
    common& cmn;

    common_read(
      common& cmn_)
    :
      cmn(cmn_)
    {}

    read_loop
    operator()(
      int unit,
      unformatted_type const&)
    {
      read_loop result(cmn, unit, unformatted);
      return result;
    }

    read_loop
    operator()(
      int unit,
      star_type const&)
    {
      read_loop result(cmn, unit, star);
      return result;
    }

    read_loop
    operator()(
      int const& unit,
      str_cref fmt)
    {
      read_loop result(cmn, unit, fmt);
      return result;
    }

    read_loop
    operator()(
      str_cref const& internal_file,
      star_type const&)
    {
      read_loop result(internal_file, star);
      return result;
    }

    read_loop
    operator()(
      str_cref const& internal_file,
      str_cref fmt)
    {
      read_loop result(internal_file, fmt);
      return result;
    }
  };

} // namespace fem

#endif // GUARD
