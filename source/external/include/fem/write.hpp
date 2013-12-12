#ifndef FEM_WRITE_HPP
#define FEM_WRITE_HPP

#include <fem/common.hpp>
#include <fem/format.hpp>
#include <fem/star.hpp>
#include <fem/str_arr_ref.hpp>
#include <fem/utils/double_to_string.hpp>
#include <fem/utils/misc.hpp>
#include <fem/utils/real_as_string.hpp>

#if defined(_MSC_VER)
# define FEM_WRITE_CRLF true
#else
# define FEM_WRITE_CRLF false
#endif

namespace fem {

  struct write_loop_base
  {
    bool write_crlf;
    unsigned pos;
    bool prev_was_string;
    int exp_scale;
    unsigned number_of_x_held;
    bool suppress_new_line_at_end;

    write_loop_base(
      bool write_crlf_=false)
    :
      write_crlf(write_crlf_),
      pos(0),
      prev_was_string(false),
      exp_scale(0),
      number_of_x_held(0),
      suppress_new_line_at_end(false)
    {}
  };

  class write_loop : write_loop_base
    // TODO copy-constructor potential performance problem
  {
    private:
      utils::slick_ptr<utils::simple_ostream> out;
      int internal_file_len;
      io_modes io_mode;
      format::token_loop fmt_loop;
      bool terminated_by_colon;

    public:

      write_loop(
        common& cmn,
        int const& unit,
        unformatted_type const&)
      :
        write_loop_base(FEM_WRITE_CRLF && !is_std_io_unit(unit)),
        out(cmn.io.simple_ostream(unit)),
        internal_file_len(-1),
        io_mode(io_unformatted),
        terminated_by_colon(false)
      {}

      write_loop(
        common& cmn,
        int const& unit,
        star_type const&)
      :
        write_loop_base(FEM_WRITE_CRLF && !is_std_io_unit(unit)),
        out(cmn.io.simple_ostream(unit)),
        internal_file_len(-1),
        io_mode(io_list_directed),
        terminated_by_colon(false)
      {}

      write_loop(
        common& cmn,
        int const& unit,
        str_cref fmt)
      :
        write_loop_base(FEM_WRITE_CRLF && !is_std_io_unit(unit)),
        out(cmn.io.simple_ostream(unit)),
        internal_file_len(-1),
        io_mode(io_formatted),
        fmt_loop(fmt),
        terminated_by_colon(false)
      {}

      write_loop(
        str_ref const& internal_file,
        star_type const&)
      :
        out(utils::slick_ptr<utils::simple_ostream>(new
          utils::simple_ostream_to_char_ptr_and_size(
            internal_file.elems(), internal_file.len()))),
        internal_file_len(internal_file.len()),
        io_mode(io_list_directed),
        terminated_by_colon(false)
      {}

      write_loop(
        str_ref const& internal_file,
        str_cref fmt)
      :
        out(utils::slick_ptr<utils::simple_ostream>(new
          utils::simple_ostream_to_char_ptr_and_size(
            internal_file.elems(), internal_file.len()))),
        internal_file_len(internal_file.len()),
        io_mode(io_formatted),
        fmt_loop(fmt),
        terminated_by_colon(false)
      {}

      template <size_t Ndims>
      write_loop(
        str_arr_ref<Ndims> const& internal_file,
        star_type const&)
      :
        out(utils::slick_ptr<utils::simple_ostream>(new
          utils::simple_ostream_to_char_ptr_and_size(
            internal_file.begin(), internal_file.len()))),
        internal_file_len(internal_file.len()),
        io_mode(io_list_directed),
        terminated_by_colon(false)
      {}

      template <size_t Ndims>
      write_loop(
        str_arr_ref<Ndims> const& internal_file,
        str_cref fmt)
      :
        out(utils::slick_ptr<utils::simple_ostream>(new
          utils::simple_ostream_to_char_ptr_and_size(
            internal_file.begin(), internal_file.len()))),
        internal_file_len(internal_file.len()),
        io_mode(io_formatted),
        fmt_loop(fmt),
        terminated_by_colon(false)
      {}

      std::string const&
      next_edit_descriptor(
        bool final=false)
      {
        while (true) {
          if(terminated_by_colon) {
            static const std::string empty("");
            return empty;
          }
          utils::token const* t = fmt_loop.next_executable_token(final);
          if (t == 0) {
            static const std::string empty("");
            return empty;
          }
          std::string const& tv = t->value;
          if (t->type == "string") {
            to_stream_fmt(tv.data(), tv.size());
          }
          else if (t->type == "op") {
            if (tv[0] == ':') {
              if(final) terminated_by_colon = true;
            }
            else if (tv[0] == '/') {
              if (write_crlf) {
                to_stream_fmt("\r\n", 2);
              }
              else {
                to_stream_fmt("\n", 1);
              }
            }
            else if (tv[0] == '$') {
              suppress_new_line_at_end = true;
            }
            else {
              out.reset();
              throw TBXX_UNREACHABLE_ERROR();
            }
          }
          else if (t->type == "format") {
            if (utils::ends_with_char(tv, 'x')) {
              unsigned n = tv.size();
              if (n != 1) n = utils::unsigned_integer_value(
                tv.data(), n-1);
              number_of_x_held += n;
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
              out.reset();
              throw TBXX_NOT_IMPLEMENTED();
            }
            else if (tv[0] == 's') {
              out.reset();
              throw TBXX_NOT_IMPLEMENTED();
            }
            else if (tv[0] == 'b') {
              out.reset();
              throw TBXX_NOT_IMPLEMENTED();
            }
            else {
              out.reset();
              throw TBXX_UNREACHABLE_ERROR();
            }
          }
          else {
            out.reset();
            throw TBXX_UNREACHABLE_ERROR();
          }
        }
      }

      write_loop&
      operator,(
        char const& val)
      {
        if (io_mode == io_unformatted) {
          to_stream_unformatted(&val, 1);
        }
        else if (io_mode == io_list_directed) {
          to_stream(&val, 1, /*space*/ !prev_was_string);
          prev_was_string = true;
        }
        else {
          out.reset();
          throw TBXX_NOT_IMPLEMENTED();
        }
        return *this;
      }

      write_loop&
      operator,(
        char const* val)
      {
        if (io_mode == io_unformatted) {
          to_stream_unformatted(val, std::strlen(val));
        }
        else if (io_mode == io_list_directed) {
          to_stream_star(val, std::strlen(val), /*space*/ !prev_was_string);
          prev_was_string = true;
        }
        else {
          std::string const& ed = next_edit_descriptor();
          int l = std::strlen(val);
          if (ed[0] == 'a') {
            int n = ed.size();
            if (n > 1) {
              n = utils::unsigned_integer_value(ed.data()+1, n-1);
              int b = n - l;
              if (b < 0) {
                l = n;
              }
              else {
                for(int i=0;i<b;i++) to_stream_fmt(" ", 1);
              }
            }
            to_stream_fmt(val, l);
          }
          else {
            to_stream_fmt(val, l);
          }
        }
        return *this;
      }

      write_loop&
      operator,(
        bool const& val)
      {
        if (io_mode == io_unformatted) {
          to_stream_unformatted(
            reinterpret_cast<char const*>(&val),
            sizeof(bool));
        }
        else if (io_mode == io_list_directed) {
          to_stream_star((val ? "T" : "F"), 1);
          prev_was_string = false;
        }
        else {
          std::string const& ed = next_edit_descriptor();
          int n = ed.size();
          if (ed[0] == 'l' && n > 1) {
            n = utils::unsigned_integer_value(ed.data()+1, n-1);
          }
          else {
            n = 1;
          }
          for(int i=1;i<n;i++) to_stream_fmt(" ", 1);
          to_stream_fmt((val ? "T" : "F"), 1);
        }
        return *this;
      }

      write_loop&
      operator,(
        integer_star_1 const& val)
      {
        if (io_mode == io_unformatted) {
          to_stream_unformatted(
            reinterpret_cast<char const*>(&val),
            sizeof(integer_star_1));
        }
        else if (io_mode == io_list_directed) {
          char buf[64];
          int n = std::sprintf(buf, "%4d", static_cast<int>(val));
          to_stream_star(buf, n);
          prev_was_string = false;
        }
        else {
          out.reset();
          throw TBXX_NOT_IMPLEMENTED();
        }
        return *this;
      }

      write_loop&
      operator,(
        integer_star_2 const& val)
      {
        if (io_mode == io_unformatted) {
          to_stream_unformatted(
            reinterpret_cast<char const*>(&val),
            sizeof(integer_star_2));
        }
        else if (io_mode == io_list_directed) {
          char buf[64];
          int n = std::sprintf(buf, "%6d", static_cast<int>(val));
          to_stream_star(buf, n);
          prev_was_string = false;
        }
        else {
          out.reset();
          throw TBXX_NOT_IMPLEMENTED();
        }
        return *this;
      }

      write_loop&
      operator,(
        integer_star_4 const& val)
      {
        if (io_mode == io_unformatted) {
          to_stream_unformatted(
            reinterpret_cast<char const*>(&val),
            sizeof(integer_star_4));
        }
        else if (io_mode == io_list_directed) {
          char buf[64];
          int n = std::sprintf(buf, "%11d", val);
          to_stream_star(buf, n);
          prev_was_string = false;
        }
        else {
          std::string const& ed = next_edit_descriptor();
          if (ed[0] == 'i') {
            int n = ed.size();
            TBXX_ASSERT(n+2 < 64);
            char fmt[64];
            fmt[0] = '%';
            std::strncpy(fmt+1, ed.data()+1, n-1);
            fmt[n] = 'd';
            fmt[n+1] = '\0';
            char buf[64];
            n = std::sprintf(buf, fmt, val);
            to_stream_fmt(buf, n);
          }
          else {
            char buf[64];
            int n = std::sprintf(buf, " %d", val);
            to_stream_fmt(buf, n);
          }
        }
        return *this;
      }

      write_loop&
      operator,(
        integer_star_8 const& val)
      {
        if (io_mode == io_unformatted) {
          to_stream_unformatted(
            reinterpret_cast<char const*>(&val),
            sizeof(integer_star_8));
        }
        else if (io_mode == io_list_directed) {
          // TODO faster implementation
          std::ostringstream o;
          o.width(21);
          o << val;
          std::string s = o.str();
          to_stream_star(s.data(), s.size());
          prev_was_string = false;
        }
        else {
          out.reset();
          throw TBXX_NOT_IMPLEMENTED();
        }
        return *this;
      }

      protected: // implementation detail
        void
        to_stream_fmt_double(
          double const& val)
        {
          std::string const& ed = next_edit_descriptor();
          if (ed[0] == 'f') {
            int n = ed.size();
            TBXX_ASSERT(n+2 < 64);
            char fmt[64];
            fmt[0] = '%';
            std::strncpy(fmt+1, ed.data()+1, n-1);
            fmt[n] = 'f';
            fmt[n+1] = '\0';
            char buf[64];
            n = std::sprintf(buf, fmt, val);
            to_stream_fmt(buf, n);
          }
          else if ((ed[0] == 'd' || ed[0] == 'e') && ed.size() > 1) {
            int es = ed.size();
            int nw = utils::unsigned_integer_scan(ed.data(), 1, es);
            TBXX_ASSERT(nw > 0);
            int w = utils::unsigned_integer_value(ed.data(), 1, nw);
            int d = 0;
            if (nw != es) {
              TBXX_ASSERT(ed[nw] == '.');
              TBXX_ASSERT(nw+1 != es);
              d = utils::unsigned_integer_value(ed.data(), nw+1, es);
            }
            utils::double_to_string_scientific_notation conv(
              val, w, d, exp_scale, utils::to_upper(ed[0]));
            to_stream_fmt(conv.buffer, w);
          }
          else {
            char buf[64];
            int n = std::sprintf(buf, " %.6g", val);
            to_stream_fmt(buf, n);
          }
        }
        public:

      write_loop&
      operator,(
        float const& val)
      {
        if (io_mode == io_unformatted) {
          to_stream_unformatted(
            reinterpret_cast<char const*>(&val),
            sizeof(float));
        }
        else if (io_mode == io_list_directed) {
          utils::float_as_string_list_directed conv(val);
          to_stream(conv.begin, conv.n);
          prev_was_string = false;
        }
        else {
          to_stream_fmt_double(static_cast<double>(val));
        }
        return *this;
      }

      write_loop&
      operator,(
        double const& val)
      {
        if (io_mode == io_unformatted) {
          to_stream_unformatted(
            reinterpret_cast<char const*>(&val),
            sizeof(double));
        }
        else if (io_mode == io_list_directed) {
          utils::double_as_string_list_directed conv(val);
          to_stream(conv.begin, conv.n);
          prev_was_string = false;
        }
        else {
          to_stream_fmt_double(val);
        }
        return *this;
      }

      write_loop&
      operator,(
        long double const& val)
      {
        if (io_mode == io_unformatted) {
          to_stream_unformatted(
            reinterpret_cast<char const*>(&val),
            actual_sizeof_long_double);
        }
        else if (io_mode == io_list_directed) {
          out.reset();
          throw TBXX_NOT_IMPLEMENTED();
        }
        else {
          out.reset();
          throw TBXX_NOT_IMPLEMENTED();
        }
        return *this;
      }

      write_loop&
      operator,(
        std::complex<float> const& val)
      {
        if (io_mode == io_unformatted) {
          float re = val.real();
          float im = val.imag();
          to_stream_unformatted(
            reinterpret_cast<char const*>(&re),
            sizeof(float));
          to_stream_unformatted(
            reinterpret_cast<char const*>(&im),
            sizeof(float));
        }
        else if (io_mode == io_list_directed) {
          utils::float_as_string_list_directed conv_re(val.real());
          utils::float_as_string_list_directed conv_im(val.imag());
          to_stream_star_complex(
            conv_re.begin, conv_re.n,
            conv_im.begin, conv_im.n);
        }
        else {
          out.reset();
          throw TBXX_NOT_IMPLEMENTED();
        }
        return *this;
      }

      write_loop&
      operator,(
        std::complex<double> const& val)
      {
        if (io_mode == io_unformatted) {
          double re = val.real();
          double im = val.imag();
          to_stream_unformatted(
            reinterpret_cast<char const*>(&re),
            sizeof(double));
          to_stream_unformatted(
            reinterpret_cast<char const*>(&im),
            sizeof(double));
        }
        else if (io_mode == io_list_directed) {
          utils::double_as_string_list_directed conv_re(val.real());
          utils::double_as_string_list_directed conv_im(val.imag());
          to_stream_star_complex(
            conv_re.begin, conv_re.n,
            conv_im.begin, conv_im.n);
        }
        else {
          out.reset();
          throw TBXX_NOT_IMPLEMENTED();
        }
        return *this;
      }

      write_loop&
      operator,(
        std::complex<long double> const& val)
      {
        if (io_mode == io_unformatted) {
          long double re = val.real();
          long double im = val.imag();
          to_stream_unformatted(
            reinterpret_cast<char const*>(&re),
            actual_sizeof_long_double);
          to_stream_unformatted(
            reinterpret_cast<char const*>(&im),
            actual_sizeof_long_double);
        }
        else if (io_mode == io_list_directed) {
          out.reset();
          throw TBXX_NOT_IMPLEMENTED();
        }
        else {
          out.reset();
          throw TBXX_NOT_IMPLEMENTED();
        }
        return *this;
      }

      write_loop&
      operator,(
        str_cref const& val)
      {
        if (io_mode == io_unformatted) {
          to_stream_unformatted(val.elems(), val.len());
        }
        else if (io_mode == io_list_directed) {
          to_stream(val.elems(), val.len(), /*space*/ !prev_was_string);
          prev_was_string = true;
        }
        else {
          std::string const& ed = next_edit_descriptor();
          int n = ed.size();
          if (ed[0] == 'a' && n > 1) {
            n = utils::unsigned_integer_value(ed.data()+1, n-1);
            to_stream(val.elems(), std::min(val.len(), n));
            for(int i=val.len();i<n;i++) {
              to_stream(" ", 1);
            }
          }
          else {
            to_stream(val.elems(), val.len());
          }
        }
        return *this;
      }

      write_loop&
      operator,(
        str_addends const& val)
      {
        int ll = val.lhs.len();
        int rl = val.rhs.len();
        utils::simple_buffer<char> buffer(ll + rl);
        char* b = buffer.space;
        std::memcpy(b, val.lhs.elems(), ll);
        std::memcpy(b+ll, val.rhs.elems(), rl);
        (*this), str_cref(b, ll+rl);
        return *this;
      }

      template <typename T, size_t Ndims>
      write_loop&
      operator,(
        arr_cref<T, Ndims> const& val)
      {
        size_t n = val.size_1d();
        T const* val_begin = val.begin();
        for(size_t i=0;i<n;i++) {
          (*this), val_begin[i];
        }
        return *this;
      }

      template <size_t Ndims>
      write_loop&
      operator,(
        str_arr_cref<Ndims> const& val)
      {
        size_t n = val.size_1d();
        int l = val.len();
        char const* val_begin = val.begin();
        for(size_t i=0;i<n;i++) {
          (*this), str_cref(&val_begin[i*l], l);
        }
        return *this;
      }

      ~write_loop()
      {
        if (out.get() == 0) return;
        if (internal_file_len < 0) {
          if (io_mode == io_unformatted) {
            out->put(end_of_unformatted_record);
            out->put('\0');
          }
          else {
            if (io_mode == io_list_directed) {
              if (pos == 0) out->put(' ');
            }
            else {
              next_edit_descriptor(/*final*/ true);
            }
            if (!suppress_new_line_at_end) {
              if (write_crlf) {
                out->put("\r\n", 2);
              }
              else {
                out->put('\n');
              }
            }
          }
          out->flush();
        }
        else {
          if (io_mode == io_unformatted) {
            out.reset();
            throw TBXX_NOT_IMPLEMENTED();
          }
          else {
            if (io_mode == io_list_directed) {
              if (pos == 0) {
                out->put(' ');
                pos++;
              }
            }
            else {
              next_edit_descriptor(/*final*/ true);
            }
            while (pos < internal_file_len) {
              out->put(' ');
              pos++;
            }
          }
        }
      }

    private:

      void
      to_stream(
        char const* buf,
        unsigned n,
        bool space=true)
      {
        switch (io_mode) {
          case io_unformatted:   to_stream_unformatted(buf, n); break;
          case io_list_directed: to_stream_star(buf, n, space); break;
          default:               to_stream_fmt(buf, n);
        }
      }

      void
      to_stream_unformatted(
        char const* buf,
        unsigned n)
      {
        for(unsigned i=0;i<n;i++) {
          char c = buf[i];
          out->put(c);
          if (c == end_of_unformatted_record) {
            out->put(c);
          }
        }
      }

      void
      to_stream_fmt(
        char const* buf,
        unsigned n)
      {
        while (number_of_x_held != 0) {
          out->put(" ", 1);
          number_of_x_held--;
        }
        out->put(buf, n);
      }

      void
      to_stream_star(
        char const* buf,
        unsigned n,
        bool space=true)
      {
        if (pos == 0) {
          out->put(' ');
          pos = 1;
        }
        else if (pos + (space ? 1 : 0) + n > 80) {
          if (write_crlf) {
            out->put("\r\n ", 3);
          }
          else {
            out->put("\n ", 2);
          }
          pos = 1;
        }
        else if (space) {
          out->put(' ');
          pos++;
        }
        out->put(buf, n);
        pos += n;
      }

      void
      to_stream_star_strip_leading_and_trailing_blank_padding(
        char const* buf,
        unsigned n,
        bool space=true)
      {
        size_t_2 indices = utils::find_leading_and_trailing_blank_padding(
          buf, n);
        to_stream_star(
          buf+indices.elems[0],
          indices.elems[1]-indices.elems[0],
          /*space*/ false);
      }

      void
      to_stream_star_complex(
        char const* buf_re, unsigned n_re,
        char const* buf_im, unsigned n_im)
      {
        to_stream_star("(", 1);
        {
          to_stream_star_strip_leading_and_trailing_blank_padding(
            buf_re, n_re, /*space*/ false);
        }
        to_stream_star(",", 1, /*space*/ false);
        {
          to_stream_star_strip_leading_and_trailing_blank_padding(
            buf_im, n_im, /*space*/ false);
        }
        to_stream_star(")", 1, /*space*/ false);
        prev_was_string = false;
      }
  };

  struct common_write
  {
    common& cmn;

    common_write(
      common& cmn_)
    :
      cmn(cmn_)
    {}

    write_loop
    operator()(
      int unit,
      unformatted_type const&)
    {
      write_loop result(cmn, unit, unformatted);
      return result;
    }

    write_loop
    operator()(
      int unit,
      star_type const&)
    {
      write_loop result(cmn, unit, star);
      return result;
    }

    write_loop
    operator()(
      int const& unit,
      str_cref fmt)
    {
      write_loop result(cmn, unit, fmt);
      return result;
    }

    write_loop
    operator()(
      str_ref const& internal_file,
      star_type const&)
    {
      write_loop result(internal_file, star);
      return result;
    }

    write_loop
    operator()(
      str_ref const& internal_file,
      str_cref fmt)
    {
      write_loop result(internal_file, fmt);
      return result;
    }

    template <size_t Ndims>
    write_loop
    operator()(
      str_arr_ref<Ndims> const& internal_file,
      star_type const&)
    {
      write_loop result(internal_file, star);
      return result;
    }

    template <size_t Ndims>
    write_loop
    operator()(
      str_arr_ref<Ndims> const& internal_file,
      str_cref fmt)
    {
      write_loop result(internal_file, fmt);
      return result;
    }
  };

} // namespace fem

#endif // GUARD
