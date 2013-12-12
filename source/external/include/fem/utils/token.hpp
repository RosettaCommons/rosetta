#ifndef FEM_UTILS_TOKEN_HPP
#define FEM_UTILS_TOKEN_HPP

#include <string>

namespace fem { namespace utils {

  struct token
  {
    std::string type;
    std::string value;

    token(
      std::string type_,
      std::string value_)
    :
      type(type_),
      value(value_)
    {}

    token(
      std::string type_,
      char value_)
    :
      type(type_),
      value(1, value_)
    {}
  };

}} // namespace fem::utils

#endif // GUARD
