#ifndef FEM_CLOSE_CHAIN_HPP
#define FEM_CLOSE_CHAIN_HPP

#include <fem/io.hpp>

namespace fem {

  struct close_chain
  {
    mutable io* io_ptr;
    int unit;
    int* iostat_ptr;
    bool status_delete;

    private:
      close_chain const& operator=(close_chain const&); // not implemented

      public:

    close_chain(
      io* io_ptr_,
      int unit_)
    :
      io_ptr(io_ptr_),
      unit(unit_),
      iostat_ptr(0),
      status_delete(false)
    {}

    close_chain(
      close_chain const& other)
    :
      io_ptr(other.io_ptr),
      unit(other.unit),
      iostat_ptr(other.iostat_ptr),
      status_delete(other.status_delete)
    {
      other.io_ptr = 0;
    }

    ~close_chain()
    {
      if (io_ptr == 0) return;
      io_unit* u_ptr = io_ptr->unit_ptr(unit);
      if (u_ptr != 0) {
        u_ptr->close(iostat_ptr, status_delete);
        io_ptr->units.erase(io_ptr->units.find(unit));
      }
    }

    close_chain&
    status(
      std::string const& val)
    {
      static const char* keywords[] = {"keep", "delete", 0};
      int i = utils::keyword_index(keywords, val, "CLOSE STATUS");
      status_delete = (i == 1);
      return *this;
    }

    close_chain&
    iostat(
      int& val)
    {
      iostat_ptr = &val;
      return *this;
    }
  };

  inline
  close_chain
  io::close(
    int unit)
  {
    return close_chain(this, unit);
  }

} // namespace fem

#endif // GUARD
