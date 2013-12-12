#ifndef FEM_ZERO_HPP
#define FEM_ZERO_HPP

#include <fem/data_type_star.hpp>
#include <complex>

namespace fem {

  /*! Herb Sutter: Why Not Specialize Function Templates?
      C/C++ Users Journal, 19(7), July 2001.
      http://www.gotw.ca/publications/mill17.htm
   */
  template <typename T>
  struct zero_impl;

  template <typename T>
  T
  zero() { return zero_impl<T>::get(); }

#define FEM_ZERO_IMPL(T, ZERO) \
  template <> \
  struct zero_impl<T > \
  { \
    static T get() { return ZERO; } \
  };

  FEM_ZERO_IMPL(char, '\0')
  FEM_ZERO_IMPL(logical_star_1, false)
  FEM_ZERO_IMPL(integer_star_1, 0)
  FEM_ZERO_IMPL(integer_star_2, 0)
  FEM_ZERO_IMPL(integer_star_4, 0)
  FEM_ZERO_IMPL(integer_star_8, 0)
  FEM_ZERO_IMPL(real_star_4, 0.f)
  FEM_ZERO_IMPL(real_star_8, 0.)
  FEM_ZERO_IMPL(real_star_16, 0.L)
#define FEM_LOC std::complex<float>
  FEM_ZERO_IMPL(FEM_LOC, 0.f)
#undef FEM_LOC
#define FEM_LOC std::complex<double>
  FEM_ZERO_IMPL(FEM_LOC, 0.)
#undef FEM_LOC
#define FEM_LOC std::complex<long double>
  FEM_ZERO_IMPL(FEM_LOC, 0.L)
#undef FEM_LOC

  static const char byte0 = zero<char>();
  static const char char0 = zero<char>();
  static const bool bool0 = zero<bool>();
  static const int int0 = zero<int>();
  static const float float0 = zero<float>();
  static const double double0 = zero<double>();
  static const long double long_double0 = zero<long double>();
  static const std::complex<float> cfloat0;
  static const std::complex<double> cdouble0;

  enum fill0_type { fill0 };
  enum no_fill0_type { no_fill0 };

} // namespace fem

#endif // GUARD
