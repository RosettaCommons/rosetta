#ifndef FEM_FORMAT_HPP
#define FEM_FORMAT_HPP

#include <fem/io_exceptions.hpp>
#include <fem/str_ref.hpp>
#include <fem/utils/char.hpp>
#include <fem/utils/misc.hpp>
#include <fem/utils/token.hpp>
#include <vector>
#include <stdexcept>
#include <cstring>

namespace fem { namespace format {

  struct tokenizer : utils::noncopyable
  {
    protected:
      char* code;
      char* str_buf;
      unsigned stop;
      unsigned i;
      public:

    std::vector<utils::token> tokens;

    tokenizer(
      char const* fmt,
      unsigned fmt_stop)
    {
      utils::simple_buffer<char> buffer(fmt_stop*2);
      code = buffer.space;
      str_buf = code + fmt_stop;
      stop = 0;
      for(i=0;i<fmt_stop;i++) {
        char c = fmt[i];
        if (c == ' ' || c == '\t') continue;
        if (c != '\'' && c != '"') {
          code[stop++] = utils::to_lower(c);
          continue;
        }
        code[stop++] = c;
        char opening_quote = c;
        for(i++;i<fmt_stop;i++) {
          c = fmt[i];
          code[stop++] = c;
          if (c == opening_quote) {
            if (i+1 == fmt_stop) break;
            if (fmt[i+1] != c) break;
            i++;
            code[stop++] = c;
          }
        }
      }
      if (stop == 0) {
        throw std::runtime_error("Invalid FORMAT specification: empty string");
      }
      stop--;
      if (code[0] != '(') raise_invalid();
      if (code[stop] != ')') raise_invalid();
      tokens.reserve(32); // uncritical; avoids reallocation in most cases
      i = 1;
      while (i < stop) process();
      code = 0;
      str_buf = 0;
    }

    protected:
      void
      raise_invalid()
      {
        throw std::runtime_error(
          "Invalid FORMAT specification: " + std::string(code, stop+1));
      }

      bool
      starts_with(
        char const* substr,
        unsigned start)
      {
        return utils::starts_with(code, start, stop, substr);
      }

      int
      unsigned_integer_scan(
        unsigned start)
      {
        return utils::unsigned_integer_scan(code, start, stop);
      }

      void
      add_token(
        char const* type,
        unsigned start)
      {
        tokens.push_back(utils::token(type,
          std::string(&code[start], &code[i])));
      }

      void
      add_token_string()
      {
        unsigned str_size = 0;
        char opening_quote = code[i];
        for(i++;i<stop;i++) {
          char c = code[i];
          if (c == opening_quote) {
            i++;
            if (i == stop || code[i] != opening_quote) {
              tokens.push_back(utils::token("string",
                std::string(str_buf, str_size)));
              return;
            }
          }
          str_buf[str_size++] = c;
        }
        raise_invalid();
      }

      void
      process()
      {
        unsigned i_code = i;
        char c = code[i_code];
        if (c == ',') {
          i++;
          return;
        }
        if (c == 'x') {
          i++;
          add_token("format", i_code);
          return;
        }
        if (std::strchr("():/$", c) != 0) {
          i++;
          add_token("op", i_code);
          return;
        }
        if (c == '\'' || c == '"') {
          add_token_string();
          return;
        }
        if (c == '+' || c == '-') {
          int j = unsigned_integer_scan(i_code+1);
          if (j < 0 || code[j] != 'p') raise_invalid();
          i = j + 1;
          add_token("format", i_code);
          return;
        }
        int j = unsigned_integer_scan(i_code);
        if (j > 0) {
          i = j;
          if (starts_with("h", j)) {
            throw std::runtime_error(
              "FATAL: Not supported: FORMAT Hollerith edit descriptor: "
                + std::string(code, stop+1));
          }
          if (starts_with("x", i) || starts_with("p", i)) {
            i++;
            add_token("format", i_code);
            return;
          }
          add_token("integer", i_code);
          return;
        }
        if (std::strchr("defgiz", c) != 0) {
          i++;
          j = unsigned_integer_scan(i);
          if (j > 0) {
            if (starts_with(".", j)) {
              j = unsigned_integer_scan(j+1);
              if (j < 0) raise_invalid();
            }
            i = j;
          }
          add_token("format", i_code);
          return;
        }
        if (c == 'a' || c == 'l') {
          i++;
          j = unsigned_integer_scan(i_code+1);
          if (j > 0) i = j;
          add_token("format", i_code);
          return;
        }
        if (starts_with("bn", i_code) || starts_with("bz", i_code)) {
          i += 2;
          add_token("format", i_code);
          return;
        }
        if (c == 's') {
          i++;
          if (starts_with("p", i) || starts_with("s", i)) {
            i++;
          }
          add_token("format", i_code);
          return;
        }
        if (c == 't') {
          i++;
          if (starts_with("l", i) || starts_with("r", i)) {
            i++;
          }
          j = unsigned_integer_scan(i);
          if (j < 0) raise_invalid();
          i = j;
          add_token("format", i_code);
          return;
        }
        raise_invalid();
      }
  };

  struct repeat_point
  {
    unsigned i_fmt;
    unsigned n;
    bool wait_for_closing_parenthesis;

    repeat_point(
      unsigned i_fmt_,
      unsigned n_,
      bool wait_for_closing_parenthesis_=false)
    :
      i_fmt(i_fmt_),
      n(n_),
      wait_for_closing_parenthesis(wait_for_closing_parenthesis_)
    {}
  };

  struct token_loop
  {
    std::vector<utils::token> fmt_tokens;
    unsigned i_fmt;
    unsigned i_fmt_wrap;
    unsigned simple_repeat;
    std::vector<format::repeat_point> repeat_points;

    token_loop()
    :
      i_fmt(0),
      i_fmt_wrap(0),
      simple_repeat(0)
    {}

    token_loop(
      str_cref fmt)
    :
      i_fmt(0),
      i_fmt_wrap(0),
      simple_repeat(0)
    {
      format::tokenizer tz(fmt.elems(), fmt.len());
      fmt_tokens.swap(tz.tokens);
      repeat_points.reserve(32); // uncritical; avoids reallocation in most
    }                            // cases

    utils::token const*
    next_executable_token(
      bool final=false)
    {
      if (simple_repeat != 0) {
        simple_repeat--;
        i_fmt--;
        return &fmt_tokens[i_fmt++];
      }
      while (true) {
        if (i_fmt == fmt_tokens.size()) {
          if (final) {
            return 0;
          }
          if (fmt_tokens.size() == 0) {
            throw io_err("Empty format string but data editing requested.");
          }
          i_fmt = i_fmt_wrap;
          static const utils::token op_slash("op", "/");
          return &op_slash;
        }
        utils::token const* t = &fmt_tokens[i_fmt++];
        std::string const& tv = t->value;
        if (t->type == "integer") {
          if (i_fmt == fmt_tokens.size()) {
            throw std::runtime_error(
              "Trailing lone repeat count in format string.");
          }
          unsigned n = utils::unsigned_integer_value(tv.data(), tv.size());
          if (n == 0) {
            throw std::runtime_error(
              "Zero repeat count in format string.");
          }
          t = &fmt_tokens[i_fmt++];
          if (t->type == "op" && t->value == "(") {
            repeat_points.push_back(format::repeat_point(i_fmt, n));
            if (repeat_points.size() == 1) {
              i_fmt_wrap = i_fmt - 2;
            }
          }
          else {
            simple_repeat = n-1;
            return t;
          }
        }
        else {
          if (t->type == "op") {
            char tv0 = tv[0];
            if (tv0 == '(') {
              repeat_points.push_back(format::repeat_point(i_fmt, 1));
              if (repeat_points.size() == 1) {
                i_fmt_wrap = i_fmt - 1;
              }
            }
            else if (tv0 == ')') {
              if (repeat_points.size() == 0) {
                throw std::runtime_error(
                  "Unexpected closing parenthesis in format string.");
              }
              repeat_point& rp = repeat_points.back();
              rp.n--;
              if (rp.n == 0) {
                repeat_points.pop_back();
              }
              else {
                i_fmt = rp.i_fmt;
              }
            }
            else {
              return t;
            }
          }
          else {
            return t;
          }
        }
      }
    }
  };

}} // namespace fem::format

#endif // GUARD
