#ifndef FEM_FILE_POSITIONING_CHAIN_HPP
#define FEM_FILE_POSITIONING_CHAIN_HPP

#include <fem/io.hpp>

namespace fem {

  struct file_positioning_chain
  {
    char const* io_function;
    mutable io_unit* u_ptr;
    int* iostat_ptr;

    private:
      file_positioning_chain const& operator=(
        file_positioning_chain const&); // not implemented

      public:

    file_positioning_chain(
      char const* io_function_,
      io_unit* u_ptr_)
    :
      io_function(io_function_),
      u_ptr(u_ptr_),
      iostat_ptr(0)
    {}

    file_positioning_chain(
      file_positioning_chain const& other)
    :
      io_function(other.io_function),
      u_ptr(other.u_ptr),
      iostat_ptr(other.iostat_ptr)
    {
      other.u_ptr = 0;
    }

    ~file_positioning_chain()
    {
      if (u_ptr == 0) {
        throw TBXX_NOT_IMPLEMENTED(); // should open file file
      }
      if (std::strcmp(io_function, "backspace") == 0) {
        u_ptr->backspace(iostat_ptr);
      }
      else if (std::strcmp(io_function, "endfile") == 0) {
        u_ptr->endfile(iostat_ptr);
      }
      else if (std::strcmp(io_function, "rewind") == 0) {
        u_ptr->rewind(iostat_ptr);
      }
      else {
        throw TBXX_UNREACHABLE_ERROR();
      }
    }

    file_positioning_chain&
    iostat(
      int& val)
    {
      iostat_ptr = &val;
      return *this;
    }
  };

  inline
  file_positioning_chain
  io::backspace(
    int unit)
  {
    return file_positioning_chain(
      "backspace", unit_ptr(unit, /*auto_open*/ true));
  }

  inline
  file_positioning_chain
  io::endfile(
    int unit)
  {
    return file_positioning_chain(
      "endfile", unit_ptr(unit, /*auto_open*/ true));
  }

  inline
  file_positioning_chain
  io::rewind(
    int unit)
  {
    return file_positioning_chain(
      "rewind", unit_ptr(unit, /*auto_open*/ true));
  }

} // namespace fem

#endif // GUARD
