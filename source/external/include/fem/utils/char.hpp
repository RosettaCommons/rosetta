#ifndef FEM_UTILS_CHAR_HPP
#define FEM_UTILS_CHAR_HPP

namespace fem { namespace utils {

  inline
  bool
  is_end_of_line(
    int c)
  {
    return (c == '\n');
  }

  inline
  bool
  is_whitespace(
    int c)
  {
    return (c == ' '
         || c == '\t'
         || c == '\r'
         || is_end_of_line(c));
  }

  //! Assumes ASCII or similar.
  inline
  bool
  is_digit(
    int c)
  {
    return (c >= '0' && c <= '9');
  }

  //! Assumes ASCII or similar.
  inline
  int
  digit_as_int(
    int c)
  {
    return c - '0';
  }

  inline
  char
  int_as_digit(
    int i)
  {
    return "0123456789"[i];
  }

  //! To avoid locale environment surprises (assumes ASCII or similar).
  inline
  int
  to_lower(
    int c)
  {
    if (c < 'A') return c;
    if (c > 'Z') return c;
    return c + ('a' - 'A');
  }

  //! To avoid locale environment surprises (assumes ASCII or similar).
  inline
  int
  to_upper(
    int c)
  {
    if (c < 'a') return c;
    if (c > 'z') return c;
    return c - ('a' - 'A');
  }

  //! To avoid locale environment surprises (assumes ASCII or similar).
  inline
  bool
  is_upper_a_through_z(
    int c)
  {
    return (c >= 'A' && c <= 'Z');
  }

  //! To avoid locale environment surprises (assumes ASCII or similar).
  inline
  bool
  is_lower_a_through_z(
    int c)
  {
    return (c >= 'a' && c <= 'z');
  }

  //! To avoid locale environment surprises (assumes ASCII or similar).
  inline
  bool
  is_a_through_z(
    int c)
  {
    return (is_upper_a_through_z(c) || is_lower_a_through_z(c));
  }

}} // namespace fem::utils

#endif // GUARD
