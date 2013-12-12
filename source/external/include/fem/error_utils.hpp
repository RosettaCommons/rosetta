#ifndef FEM_ERROR_UTILS_HPP
#define FEM_ERROR_UTILS_HPP

#include <tbxx/error_utils.hpp>

#define FEM_THROW_UNHANDLED(message) \
  throw std::runtime_error( \
    tbxx::error_utils::file_and_line_as_string(__FILE__, __LINE__) \
    + ": UNHANDLED: " + message)

#endif // GUARD
