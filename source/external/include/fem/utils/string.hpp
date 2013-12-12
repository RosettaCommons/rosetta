#ifndef FEM_UTILS_STRING_HPP
#define FEM_UTILS_STRING_HPP

#include <fem/size_t.hpp>
#include <fem/utils/char.hpp>
#include <tbxx/error_utils.hpp>
#include <algorithm>
#include <cstring>

namespace fem { namespace utils {

  inline
  bool
  starts_with(
    char const* str,
    unsigned start,
    unsigned stop,
    char const* substr)
  {
    for(unsigned j=start;j<stop;) {
      if (*substr == '\0') break;
      if (str[j++] != *substr++) return false;
    }
    return (start != stop);
  }

  inline
  bool
  ends_with_char(
    std::string const& str,
    int c)
  {
    unsigned i = str.size();
    if (i == 0) return false;
    return (str[i-1] == c);
  }

  // compare with fable/__init__.py
  inline
  int
  unsigned_integer_scan(
    char const* code,
    unsigned start,
    unsigned stop)
  {
    unsigned i = start;
    for(;i<stop;i++) {
      int c = code[i];
      if (!is_digit(c)) break;
    }
    if (i == start) return -1;
    return i;
  }

  //! Assumes ASCII or similar.
  inline
  unsigned
  unsigned_integer_value(
    char const* str,
    unsigned start,
    unsigned stop)
  {
    unsigned result = 0;
    unsigned i = start;
    for(;i<stop;i++) {
      result *= 10;
      result += (str[i] - '0');
    }
    return result;
  }

  inline
  unsigned
  unsigned_integer_value(
    char const* str,
    unsigned stop)
  {
    return unsigned_integer_value(str, 0, stop);
  }

  inline
  int
  signed_integer_value(
    char const* str,
    unsigned start,
    unsigned stop)
  {
    bool negative;
    if (str[start] == '-') {
      negative = true;
      start++;
    }
    else {
      negative = false;
      if (str[start] == '+') start++;
    }
    int result = unsigned_integer_value(str, start, stop);
    if (negative) result *= -1;
    return result;
  }

  inline
  void
  copy_with_blank_padding(
    char const* src,
    size_t src_size,
    char* dest,
    size_t dest_size)
  {
    if (dest_size < src_size) {
      std::memmove(dest, src, dest_size);
    }
    else {
      std::memmove(dest, src, src_size);
      for (size_t i=src_size;i<dest_size;i++) {
        dest[i] = ' ';
      }
    }
  }

  inline
  void
  copy_with_blank_padding(
    char const* src,
    char* dest,
    size_t dest_size)
  {
    size_t i;
    for (i=0; i < dest_size && *src != '\0'; i++) {
      dest[i] = *src++;
    }
    for (; i < dest_size; i++) {
      dest[i] = ' ';
    }
  }

  inline
  bool
  string_eq(
    char const* lhs,
    size_t lhs_size,
    char const* rhs,
    size_t rhs_size)
  {
    static const char blank = ' ';
    if (lhs_size < rhs_size) {
      return string_eq(rhs, rhs_size, lhs, lhs_size);
    }
    if (std::memcmp(lhs, rhs, rhs_size) != 0) return false;
    for(size_t i=rhs_size;i<lhs_size;i++) {
      if (lhs[i] != blank) return false;
    }
    return true;
  }

  inline
  bool
  string_eq(
    char const* lhs,
    size_t lhs_size,
    char const* rhs)
  {
    static const char blank = ' ';
    for(size_t i=0;i<lhs_size;i++) {
      if (*rhs == '\0') {
        for(;i<lhs_size;i++) {
          if (lhs[i] != blank) return false;
        }
        return true;
      }
      if (*rhs++ != lhs[i]) return false;
    }
    while(*rhs != '\0') {
      if (*rhs++ != blank) return false;
    }
    return true;
  }

  inline
  int
  string_compare_lexical(
    char const* lhs,
    size_t lhs_size,
    char const* rhs,
    size_t rhs_size)
  {
    size_t n = std::max(lhs_size, rhs_size);
    for(size_t i=0;i<n;i++) {
      char l = (i < lhs_size ? lhs[i] : ' ');
      char r = (i < rhs_size ? rhs[i] : ' ');
      if (l < r) return -1;
      if (l > r) return  1;
    }
    return 0;
  }

  inline
  size_t
  find_leading_blank_padding(
    char const* str,
    size_t stop)
  {
    size_t i = 0;
    while (i != stop) {
      if (str[i] != ' ') break;
      i++;
    }
    return i;
  }

  inline
  size_t
  find_trailing_blank_padding(
    char const* str,
    size_t stop)
  {
    size_t i = stop;
    while (i != 0) {
      i--;
      if (str[i] != ' ') {
        i++;
        break;
      }
    }
    return i;
  }

  inline
  size_t_2
  find_leading_and_trailing_blank_padding(
    char const* str,
    size_t stop)
  {
    return size_t_2(
      find_leading_blank_padding(str, stop),
      find_trailing_blank_padding(str, stop));
  }

  inline
  std::string
  strip_leading_and_trailing_blank_padding(
    std::string const& str)
  {
    size_t_2 indices = find_leading_and_trailing_blank_padding(
      str.data(), str.size());
    if (indices.elems[0] == 0 && indices.elems[1] == str.size()) {
      return str;
    }
    return std::string(
      str.data() + indices.elems[0],
      indices.elems[1] - indices.elems[0]);
  }

  inline
  std::string
  to_lower(
    std::string const& str)
  {
    std::string result = str;
    size_t n = str.size();
    for(size_t i=0;i<n;i++) {
      result[i] = to_lower(result[i]);
    }
    return result;
  }

  inline
  int
  keyword_index(
    char const* valid_vals[],
    std::string const& val,
    char const* throw_info=0)
  {
    std::string
      val_norm = to_lower(strip_leading_and_trailing_blank_padding(val));
    for (int i=0; valid_vals[i] != 0; i++) {
      if (std::strcmp(valid_vals[i], val_norm.c_str()) == 0) {
        return i;
      }
    }
    if (throw_info != 0) {
      std::ostringstream o;
      o << throw_info << ": invalid keyword: \"" << val << "\"";
      throw std::runtime_error(o.str());
    }
    return -1;
  }

  //! Assumes ASCII or similar.
  inline
  std::string
  format_char_for_display(
    int c)
  {
    std::ostringstream o;
    bool printable = (c >= ' ' && c <= '~');
    if (printable) {
      if (c == '"') {
        o << "'\"' (double quote, ";
      }
      else if (c == '\'') {
        o << "\"'\" (single quote, ";
      }
      else {
        o << "\"" << static_cast<char>(c) << "\" (";
      }
    }
    o << "ordinal=" << (c < 0 ? c + 256 : c);
    if (printable) o << ")";
    return o.str();
  }

  inline
  void
  string_reverse_in_place(
    char* s,
    size_t s_size)
  {
    if (s_size == 0) return;
    size_t i = 0;
    size_t j = s_size - 1;
    while (i < j) {
      std::swap(s[i], s[j]);
      i++;
      j--;
    }
  }

  //! Assumes ASCII or similar.
  inline
  int
  int_to_string(
    char* buffer,
    int buffer_size,
    int width,
    int value,
    int left_padding_character=' ')
  {
    int i = 0;
    while (value != 0) {
      if (i == buffer_size) return -1;
      buffer[i++] = int_as_digit(value % 10);
      value /= 10;
    }
    while (i < width) buffer[i++] = left_padding_character;
    string_reverse_in_place(buffer, i);
    return i;
  }

  template <typename VectorOfStringType>
  unsigned
  split_comma_separated(
    VectorOfStringType& result,
    char const* c_str)
  {
    for(unsigned i=0;;i++) {
      char c = c_str[i];
      if (c == '\0') return i;
      if (c == ',' || is_whitespace(c)) continue;
      for(unsigned i_start=i++;;i++) {
        char c = c_str[i];
        if (c == '\0' || c == ',' || is_whitespace(c)) {
          result.push_back(std::string(c_str+i_start, i-i_start));
          if (c == '\0') return i;
          break;
        }
      }
    }
  }

}} // namespace fem::utils

#endif // GUARD
