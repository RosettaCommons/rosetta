#ifndef FEM_UTILS_RANDOM_HPP
#define FEM_UTILS_RANDOM_HPP

#include <fem/size_t.hpp>
#include <string>
#include <ctime>

namespace fem { namespace utils {

  inline
  std::string
  random_name_simple(
    size_t size)
  {
    static bool first_call = true;
    static size_t random_state;
    static const size_t random_mod = 225150U;
    if (first_call) {
      first_call = false;
      random_state = static_cast<size_t>(std::time(0)) % random_mod;
    }
    std::string result;
    result.reserve(size);
    for(size_t i=0;i<size;i++) {
      random_state = (random_state * 9538U + 50294U) % random_mod;
      size_t j = random_state % (i == 0 ? 26U : 36U);
      result.push_back("abcdefghijklmnopqrstuvwxyz0123456789"[j]);
    }
    return result;
  }

}} // namespace fem::utils

#endif // GUARD
