#ifndef FEM_UTILS_SIMPLE_STREAMS_HPP
#define FEM_UTILS_SIMPLE_STREAMS_HPP

#include <fem/size_t.hpp>
#include <string>
#include <cstdio>

#if defined(_MSC_VER) || _MSC_VER <= 1310 // Visual C++ 7.1
#define std__ferror ferror
#endif

namespace fem { namespace utils {

  static const int stream_end = 256;
  static const int stream_err = 257;

  inline
  bool
  is_stream_end(
    int c) { return c == stream_end; }

  inline
  bool
  is_stream_err(
    int c) { return c == stream_err; }

  inline
  bool
  is_stream_end_or_err(
    int c) { return c >= stream_end; }

  struct simple_ostream
  {
    virtual ~simple_ostream() {}
    virtual void put(char c) = 0;
    virtual void put(char const* str, size_t str_sz) = 0;
    virtual void flush() = 0;
    virtual bool err() = 0;
  };

  struct simple_istream
  {
    virtual ~simple_istream() {}
    virtual int get() = 0;
    virtual void backup() = 0;
  };

  struct simple_ostream_to_char_ptr_and_size : simple_ostream
  {
    char* s;
    size_t sz;
    ssize_t i;

    simple_ostream_to_char_ptr_and_size(
      char* string,
      size_t size) : s(string), sz(size), i(0)
    {}

    void
    put(char c)
    {
      if (i < 0) return;
      if (i < sz) s[i++] = c;
      else i = -1;
    }

    void
    put(char const* str, size_t str_sz)
    {
      if (i < 0) return;
      if (i+str_sz <= sz) {
        std::memcpy(s+i, str, str_sz);
        i += str_sz;
      }
      else {
        ssize_t j = sz - i;
        std::memcpy(s+i, str, j);
        i = -1;
      }
    }

    void
    flush() {}

    bool
    err() { return (i < 0); }
  };

  struct simple_ostream_to_std_string : simple_ostream
  {
    std::string& s;

    simple_ostream_to_std_string(
      std::string& string) : s(string)
    {}

    void
    put(char c) { s.push_back(c); }

    void
    put(char const* str, size_t str_sz)
    {
      for(size_t i=0;i<str_sz;i++) {
        s.push_back(str[i]);
      }
    }

    void
    flush() {}

    bool
    err() { return false; }
  };

  struct simple_ostream_to_c_file : simple_ostream
  {
    FILE* f;

    simple_ostream_to_c_file(
      FILE* f_): f(f_)
    {}

    void
    put(char c) { std::fputc(c, f); }

    void
    put(char const* str, size_t str_sz)
    {
      for(size_t i=0;i<str_sz;i++) {
        std::fputc(str[i], f);
      }
    }

    void
    flush() { std::fflush(f); }

    bool
    err() { return std__ferror(f); }
  };

  struct simple_istream_from_c_str : simple_istream
  {
    char const* s;

    simple_istream_from_c_str(
      char const* string) : s(string)
    {}

    int
    get()
    {
      if (*s == '\0') return stream_end;
      return *s++;
    }

    void
    backup() { s--; }
  };

  struct simple_istream_from_char_ptr_and_size : simple_istream
  {
    char const* s;
    size_t sz;
    size_t i;

    simple_istream_from_char_ptr_and_size(
      char const* string,
      size_t size) : s(string), sz(size), i(0)
    {}

    int
    get()
    {
      if (i == sz) return stream_end;
      return s[i++];
    }

    void
    backup() { i--; }
  };

  struct simple_istream_from_std_string : simple_istream
  {
    std::string s;
    size_t i;

    simple_istream_from_std_string(
      std::string const& string) : s(string), i(0)
    {}

    int
    get()
    {
      if (i == s.size()) return stream_end;
      return s[i++];
    }

    void
    backup() { i--; }
  };

  struct simple_istream_from_c_file : simple_istream
  {
    FILE* f;
    int last_get_result;

    simple_istream_from_c_file(
      FILE* f_): f(f_)
    {}

    int
    get()
    {
      last_get_result = std::fgetc(f);
      if (last_get_result == EOF) {
        last_get_result = (std__ferror(f) ? stream_err : stream_end);
      }
      return last_get_result;
    }

    void
    backup()
    {
      std::ungetc(last_get_result, f);
    }
  };

}} // namespace fem::utils

#endif // GUARD
