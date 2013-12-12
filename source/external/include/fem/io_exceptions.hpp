#ifndef FEM_IO_EXCEPTIONS_HPP
#define FEM_IO_EXCEPTIONS_HPP

#include <stdexcept>
#include <string>

namespace fem {

  struct io_err : std::runtime_error
  {
    explicit
    io_err(std::string const& msg) : std::runtime_error(msg) {}
  };

  struct read_end : std::runtime_error
  {
    explicit
    read_end(std::string const& msg) : std::runtime_error(msg) {}
  };

} // namespace fem

#endif // GUARD
