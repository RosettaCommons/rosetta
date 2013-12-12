#ifndef FEM_UTILS_INT_TYPES_HPP
#define FEM_UTILS_INT_TYPES_HPP

namespace fem { namespace utils {

  // int type names as used in boost/cstdint.hpp

  typedef signed char      int8_t;
  typedef signed short     int16_t;
  typedef signed int       int32_t;
#if defined(__i386__) || defined(__ppc__) || defined(_MSC_VER)
  typedef signed long long int64_t;
#else
  typedef signed long      int64_t;
#endif

}} // namespace fem::utils

#endif // GUARD
