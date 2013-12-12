#ifndef FEM_INTRINSICS_HPP
#define FEM_INTRINSICS_HPP

#include <fem/str_ref.hpp>
#include <cmath>

namespace fem {

  template <typename T>
  inline
  int
  if_arithmetic(
    T const& value)
  {
    if (value == 0) return 0;
    if (value > 0) return 1;
    return -1;
  }

  template <typename V, typename S>
  inline
  V
  sign(
    V const& value,
    S const& sign_source)
  {
    if (sign_source < 0) {
      if (value > 0) return -value;
    }
    else if (value < 0) {
      return -value;
    }
    return value;
  }

  inline
  double
  dsign(
    double const& value,
    double const& sign_source)
  {
    return sign(value, sign_source);
  }

  inline
  int
  isign(
    int const& value,
    int const& sign_source)
  {
    return sign(value, sign_source);
  }


  template <typename T>
  inline
  float
  real(
    T const& value) { return static_cast<float>(value); }

  template <typename T>
  inline
  float
  sngl(
    T const& value) { return static_cast<float>(value); }


  template <typename T>
  inline
  double
  dble(
    T const& value) { return static_cast<double>(value); }

  inline
  float
  sqrt(
    float const& x)
  {
    return std::sqrt(x);
  }

  inline
  double
  sqrt(
    double const& x)
  {
    return std::sqrt(x);
  }

  inline
  double
  dsqrt(
    double const& x)
  {
    return std::sqrt(x);
  }

  inline
  float
  cos(
    float const& angle)
  {
    return std::cos(angle);
  }

  inline
  double
  cos(
    double const& angle)
  {
    return std::cos(angle);
  }

  inline
  float
  acos(
    float const& angle)
  {
    return std::acos(angle);
  }

  inline
  double
  dacos(
    double const& angle)
  {
    return std::acos(angle);
  }

  inline
  double
  acos(
    double const& angle) { return dacos(angle); }

  inline
  float
  sin(
    float const& angle)
  {
    return std::sin(angle);
  }

  inline
  double
  sin(
    double const& angle)
  {
    return std::sin(angle);
  }

  inline
  float
  asin(
    float const& angle)
  {
    return std::asin(angle);
  }

  inline
  double
  dasin(
    double const& angle)
  {
    return std::asin(angle);
  }

  inline
  double
  asin(
    double const& angle) { return dasin(angle); }

  inline
  float
  tan(
    float const& angle)
  {
    return std::tan(angle);
  }

  inline
  double
  dtan(
    double const& angle)
  {
    return std::tan(angle);
  }

  inline
  double
  tan(
    double const& angle) { return dtan(angle); }

  inline
  float
  atan(
    float const& angle)
  {
    return std::atan(angle);
  }

  inline
  double
  datan(
    double const& angle)
  {
    return std::atan(angle);
  }

  inline
  double
  atan(
    double const& angle) { return datan(angle); }

  inline
  float
  atan2(
    float const& y,
    float const& x)
  {
    return std::atan2(y, x);
  }

  inline
  double
  datan2(
    double const& y,
    double const& x)
  {
    return std::atan2(y, x);
  }

  inline
  double
  atan2(
    double const& y,
    double const& x) { return datan2(y, x); }

  inline
  float
  exp(
    float const& x)
  {
    return std::exp(x);
  }

  inline
  double
  exp(
    double const& x)
  {
    return std::exp(x);
  }

  inline
  double
  dexp(
    double const& x)
  {
    return std::exp(x);
  }

  inline
  float
  alog10(
    float const& x)
  {
    return std::log10(x);
  }

  inline
  double
  alog10(
    double const& x)
  {
    return std::log10(x);
  }

  inline
  float
  alog(
    float const& x)
  {
    return std::log(x);
  }

  inline
  double
  dlog(
    double const& x)
  {
    return std::log(x);
  }

  template <typename T>
  inline
  int
  fint(
    T const& val)
  {
    return static_cast<int>(val);
  }

  template <typename T>
  inline
  int
  aint(
    T const& val)
  {
    return static_cast<int>(val);
  }

  template <typename T>
  inline
  float
  ffloat(
    T const& val)
  {
    return static_cast<float>(val);
  }

  inline
  int
  mod(
    int const& v1,
    int const& v2) { return v1 % v2; }

  inline
  float
  amod(
    float const& v1,
    float const& v2) { return std::fmod(v1, v2); }

  inline
  float
  mod(
    float const& v1,
    float const& v2) { return std::fmod(v1, v2); }

  inline
  double
  dmod(
    double const& v1,
    double const& v2) { return std::fmod(v1, v2); }

  inline
  double
  mod(
    double const& v1,
    double const& v2) { return std::fmod(v1, v2); }

  inline
  int
  iabs(
    int const& v) { return std::abs(v); }

  inline
  double
  dabs(
    double const& v) { return std::abs(v); }

  inline
  int
  min0(
    int const& v1,
    int const& v2) { return std::min(v1, v2); }

  inline
  int
  min(
    int const& v1,
    int const& v2) { return std::min(v1, v2); }

  inline
  float
  amin1(
    float const& v1,
    float const& v2) { return std::min(v1, v2); }

  inline
  float
  amin1(
    float const& v1,
    float const& v2,
    float const& v3) { return amin1(amin1(v1, v2), v3); }

  inline
  float
  amin1(
    float const& v1,
    float const& v2,
    float const& v3,
    float const& v4) { return amin1(amin1(v1, v2, v3), v4); }

  inline
  float
  amin1(
    float const& v1,
    float const& v2,
    float const& v3,
    float const& v4,
    float const& v5) { return amin1(amin1(v1, v2, v3, v4), v5); }

  inline
  float
  min(
    float const& v1,
    float const& v2) { return std::min(v1, v2); }

  inline
  double
  dmin1(
    double const& v1,
    double const& v2) { return std::min(v1, v2); }

  inline
  double
  min(
    double const& v1,
    double const& v2) { return std::min(v1, v2); }

  inline
  float
  amin0(
    int const& v1,
    float const& v2) { return std::min(static_cast<float>(v1), v2); }

  inline
  float
  min(
    int const& v1,
    float const& v2) { return std::min(static_cast<float>(v1), v2); }

  inline
  float
  min1(
    float const& v1,
    int const& v2) { return std::min(v1, static_cast<float>(v2)); }

  inline
  float
  min(
    float const& v1,
    int const& v2) { return std::min(v1, static_cast<float>(v2)); }

  inline
  int
  min(
    int const& v1,
    int const& v2,
    int const& v3)
  {
    return min(v1, min(v2, v3));
  }

  inline
  int
  min(
    int const& v1,
    int const& v2,
    int const& v3,
    int const& v4)
  {
    return min(v1, min(v2, v3, v4));
  }

  inline
  float
  min(
    float const& v1,
    float const& v2,
    float const& v3)
  {
    return min(v1, min(v2, v3));
  }

  inline
  float
  min(
    float const& v1,
    float const& v2,
    float const& v3,
    float const& v4)
  {
    return min(v1, min(v2, v3, v4));
  }

  inline
  double
  min(
    double const& v1,
    double const& v2,
    double const& v3)
  {
    return min(v1, min(v2, v3));
  }

  inline
  double
  min(
    double const& v1,
    double const& v2,
    double const& v3,
    double const& v4)
  {
    return min(v1, min(v2, v3, v4));
  }

  inline
  int
  max0(
    int const& v1,
    int const& v2) { return std::max(v1, v2); }

  inline
  int
  max(
    int const& v1,
    int const& v2) { return std::max(v1, v2); }

  inline
  float
  amax1(
    float const& v1,
    float const& v2) { return std::max(v1, v2); }

  inline
  float
  amax1(
    float const& v1,
    float const& v2,
    float const& v3) { return amax1(amax1(v1, v2), v3); }

  inline
  float
  amax1(
    float const& v1,
    float const& v2,
    float const& v3,
    float const& v4) { return amax1(amax1(v1, v2, v3), v4); }

  inline
  float
  amax1(
    float const& v1,
    float const& v2,
    float const& v3,
    float const& v4,
    float const& v5) { return amax1(amax1(v1, v2, v3, v4), v5); }

  inline
  float
  max(
    float const& v1,
    float const& v2,
    float const& v3,
    float const& v4) { return amax1(v1, v2, v3, v4); }

  inline
  float
  max(
    float const& v1,
    float const& v2) { return std::max(v1, v2); }

  inline
  double
  dmax1(
    double const& v1,
    double const& v2) { return std::max(v1, v2); }

  inline
  double
  max(
    double const& v1,
    double const& v2) { return std::max(v1, v2); }

  inline
  float
  amax0(
    int const& v1,
    float const& v2) { return std::max(static_cast<float>(v1), v2); }

  inline
  float
  max(
    int const& v1,
    float const& v2) { return std::max(static_cast<float>(v1), v2); }

  inline
  float
  max1(
    float const& v1,
    int const& v2) { return std::max(v1, static_cast<float>(v2)); }

  inline
  float
  max(
    float const& v1,
    int const& v2) { return std::max(v1, static_cast<float>(v2)); }

  inline
  double
  max(
    float const& v1,
    double const& v2) { return std::max(static_cast<double>(v1), v2); }

  inline
  int
  max(
    int const& v1,
    int const& v2,
    int const& v3) { return max(max(v1, v2), v3); }

  inline
  float
  max(
    float const& v1,
    float const& v2,
    float const& v3) { return max(max(v1, v2), v3); }

  inline
  double
  max(
    double const& v1,
    double const& v2,
    double const& v3) { return max(max(v1, v2), v3); }

  template <typename T>
  inline
  T
  pow2(
    T const& base) { return base * base; }

  template <typename T>
  inline
  T
  pow3(
    T const& base) { return base * base * base; }

  template <typename T>
  inline
  T
  pow4(
    T const& base) { T p2 = pow2(base); return p2 * p2; }

  inline
  int
  pow(
    int const& base,
    int const& exponent)
  {
    if (exponent < 0) return 0;
    int result = 1;
    for(int i=0;i<exponent;i++) {
      result *= base;
    }
    return result;
  }

  inline
  float
  pow(
    int const& base,
    float const& exponent)
  {
    return std::pow(static_cast<float>(base), exponent);
  }

  inline
  double
  pow(
    int const& base,
    double const& exponent)
  {
    return std::pow(static_cast<double>(base), exponent);
  }

  inline
  float
  pow(
    float const& base,
    int const& exponent)
  {
    return std::pow(base, exponent);
  }

  inline
  double
  pow(
    double const& base,
    int const& exponent)
  {
    return std::pow(base, exponent);
  }

  inline
  float
  pow(
    float const& base,
    float const& exponent)
  {
    return std::pow(base, exponent);
  }

  inline
  double
  pow(
    double const& base,
    float const& exponent)
  {
    return std::pow(base, static_cast<double>(exponent));
  }

  inline
  double
  pow(
    float const& base,
    double const& exponent)
  {
    return std::pow(static_cast<double>(base), exponent);
  }

  inline
  double
  pow(
    double const& base,
    double const& exponent)
  {
    return std::pow(base, exponent);
  }

  using std::abs;
  using std::log;

  inline
  int
  ichar(
    str_cref c)
  {
    if (c.len() == 0) {
      std::ostringstream o;
      o << "ichar() argument must be a one-character string,"
        << " but actual string length is " << c.len() << ".";
      throw std::runtime_error(o.str());
    }
    return static_cast<int>(c[0]);
  }

  inline
  str<1>
  fchar(
    int i) { return str<1>(static_cast<char>(i)); }

  inline
  int
  len_trim(
    str_cref c)
  {
    return static_cast<int>(
      utils::find_trailing_blank_padding(c.elems(), c.len()));
  }

  inline
  int
  index(
    str_cref sentence,
    str_cref word)
  {
    size_t i = std::string(sentence).find(std::string(word));
    if (i == std::string::npos) return 0;
    return static_cast<int>(i + 1);
  }

  inline
  int
  lnblnk(
    str_cref c) { return len_trim(c); }


  inline
  int
  nint(
    float const& value)
  {
    return static_cast<int>(std::floor(value+0.5));
  }

  inline
  int
  idnint(
    float const& value)
  {
    return static_cast<int>(std::floor(value+0.5));
  }


  inline
  int
  nint(
    double const& value)
  {
    return static_cast<int>(std::floor(value+0.5));
  }

  inline
  std::complex<float>
  cmplx(
    float const& re,
    float const& im)
  {
    return std::complex<float>(re, im);
  }

  inline
  std::complex<double>
  dcmplx(
    double const& re,
    double const& im)
  {
    return std::complex<double>(re, im);
  }

  inline
  std::complex<double>
  cmplx(
    double const& re,
    double const& im) { return dcmplx(re, im); }

  //! Assumes ASCII or similar.
  inline
  bool
  llt(
    str_cref a,
    str_cref b)
  {
    return a < b;
  }

  //! Assumes ASCII or similar.
  inline
  bool
  lle(
    str_cref a,
    str_cref b)
  {
    return a <= b;
  }

  //! Assumes ASCII or similar.
  inline
  bool
  lgt(
    str_cref a,
    str_cref b)
  {
    return a > b;
  }

  //! Assumes ASCII or similar.
  inline
  bool
  lge(
    str_cref a,
    str_cref b)
  {
    return a >= b;
  }

} // namespace fem

#endif // GUARD
