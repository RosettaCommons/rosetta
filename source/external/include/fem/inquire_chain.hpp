#ifndef FEM_INQUIRE_CHAIN_HPP
#define FEM_INQUIRE_CHAIN_HPP

#include <fem/io.hpp>
#include <fem/str_ref.hpp>

namespace fem {

  struct inquire_chain
  {
    mutable io* io_ptr;
    int unit;
    std::string file;
    bool have_file;

    private:
      inquire_chain const& operator=(inquire_chain const&); // not implemented

      public:

    inquire_chain(
      io* io_ptr_,
      int unit_)
    :
      io_ptr(io_ptr_),
      unit(unit_),
      have_file(false)
    {}

    inquire_chain(
      io* io_ptr_,
      std::string const& file_)
    :
      io_ptr(io_ptr_),
      unit(0),
      file(utils::strip_leading_and_trailing_blank_padding(file_)),
      have_file(true)
    {}

    inquire_chain(
      inquire_chain const& other)
    :
      io_ptr(other.io_ptr),
      unit(other.unit),
      file(other.file),
      have_file(other.have_file)
    {
      other.io_ptr = 0;
    }

    ~inquire_chain()
    {
      if (io_ptr == 0) return;
      // TODO
    }

    inquire_chain&
    iostat(
      int& val)
    {
      throw TBXX_NOT_IMPLEMENTED();
    }

    inquire_chain&
    exist(
      bool& val)
    {
      // f77_std 12.10.3.3
      if (have_file) {
        val = utils::path::exists(file.c_str());
      }
      else {
        val = (io_ptr->unit_ptr(unit) != 0);
      }
      return *this;
    }

    inquire_chain&
    opened(
      bool& val)
    {
      if (have_file) {
        val = io_ptr->is_opened_simple(file);
      }
      else {
        io_unit* u_ptr = io_ptr->unit_ptr(unit);
        val = (u_ptr != 0);
      }
      return *this;
    }

    inquire_chain&
    number(
      int& val)
    {
      throw TBXX_NOT_IMPLEMENTED();
    }

    inquire_chain&
    named(
      bool&)
    {
      throw TBXX_NOT_IMPLEMENTED();
    }

    inquire_chain&
    name(
      str_ref val)
    {
      if (have_file) {
        val = utils::path::absolute(file.c_str()).c_str();
      }
      else {
        io_unit* u_ptr = io_ptr->unit_ptr(unit);
        if (u_ptr == 0) {
          val = " ";
        }
        else {
          val = u_ptr->file_name.c_str();
        }
      }
      return *this;
    }

    inquire_chain&
    access(
      str_ref val)
    {
      if (have_file) {
        throw TBXX_NOT_IMPLEMENTED();
      }
      else {
        io_unit* u_ptr = io_ptr->unit_ptr(unit);
        if (u_ptr == 0) {
          val = "UNDEFINED";
        }
        else if (u_ptr->access == io_unit::ac_sequential) {
          val = "SEQUENTIAL";
        }
        else if (u_ptr->access == io_unit::ac_direct) {
          val = "DIRECT";
        }
        else {
          val = "UNDEFINED";
        }
      }
      return *this;
    }

    inquire_chain&
    sequential(
      str_ref val)
    {
      if (have_file) {
        throw TBXX_NOT_IMPLEMENTED();
      }
      else {
        io_unit* u_ptr = io_ptr->unit_ptr(unit);
        if (u_ptr != 0 && u_ptr->access == io_unit::ac_sequential) {
          val = "YES";
        }
        else {
          val = "NO";
        }
      }
      return *this;
    }

    inquire_chain&
    direct(
      str_ref val)
    {
      if (have_file) {
        throw TBXX_NOT_IMPLEMENTED();
      }
      else {
        io_unit* u_ptr = io_ptr->unit_ptr(unit);
        if (u_ptr != 0 && u_ptr->access == io_unit::ac_direct) {
          val = "YES";
        }
        else {
          val = "NO";
        }
      }
      return *this;
    }

    inquire_chain&
    form(
      str_ref)
    {
      throw TBXX_NOT_IMPLEMENTED();
    }

    inquire_chain&
    formatted(
      str_ref)
    {
      throw TBXX_NOT_IMPLEMENTED();
    }

    inquire_chain&
    unformatted(
      str_ref)
    {
      throw TBXX_NOT_IMPLEMENTED();
    }

    inquire_chain&
    recl(
      int&)
    {
      throw TBXX_NOT_IMPLEMENTED();
    }

    inquire_chain&
    nextrec(
      int&)
    {
      throw TBXX_NOT_IMPLEMENTED();
    }

    inquire_chain&
    blank(
      str_ref val)
    {
      if (have_file) {
        throw TBXX_NOT_IMPLEMENTED();
      }
      else {
        io_unit* u_ptr = io_ptr->unit_ptr(unit);
        if (u_ptr == 0) {
          val = "UNDEFINED";
        }
        else if (u_ptr->blank == io_unit::bl_null) {
          val = "NULL";
        }
        else if (u_ptr->blank == io_unit::bl_null) {
          val = "ZERO";
        }
        else {
          val = "UNDEFINED";
        }
      }
      return *this;
    }
  };

  inline
  inquire_chain
  io::inquire_unit(
    int unit)
  {
    return inquire_chain(this, unit);
  }

  inline
  inquire_chain
  io::inquire_file(
    std::string file,
    bool blank_padding_removed_already)
  {
    if (!blank_padding_removed_already) {
      file = utils::strip_leading_and_trailing_blank_padding(file);
    }
    return inquire_chain(this, file);
  }

} // namespace fem

#endif // GUARD
