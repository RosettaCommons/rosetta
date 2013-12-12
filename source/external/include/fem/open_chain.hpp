#ifndef FEM_OPEN_CHAIN_HPP
#define FEM_OPEN_CHAIN_HPP

#include <fem/io.hpp>

namespace fem {

  struct open_chain
  {
    mutable io_unit* u_ptr;
    int* iostat_ptr;

    private:
      open_chain const& operator=(open_chain const&); // not implemented

      public:

    open_chain(
      io_unit* u_ptr_)
    :
      u_ptr(u_ptr_),
      iostat_ptr(0)
    {}

    open_chain(
      open_chain const& other)
    :
      u_ptr(other.u_ptr),
      iostat_ptr(other.iostat_ptr)
    {
      other.u_ptr = 0;
    }

    ~open_chain()
    {
      if (u_ptr == 0) return;
      u_ptr->open(iostat_ptr);
    }

    open_chain&
    access(
      std::string const& val)
    {
      int i = utils::keyword_index(
        io_unit::ac_keywords(), val, "OPEN ACCESS");
      u_ptr->access = static_cast<io_unit::access_types>(i);
      return *this;
    }

    open_chain&
    form(
      std::string const& val)
    {
      int i = utils::keyword_index(
        io_unit::fm_keywords(), val, "OPEN FORM");
      u_ptr->form = static_cast<io_unit::form_types>(i);
      return *this;
    }

    open_chain&
    recl(
      int const& val)
    {
      if (val <= 0) {
        throw io_err(
          "Invalid OPEN RECL: value is zero or negative,"
          " but must be greater than zero");
      }
      u_ptr->recl = static_cast<unsigned>(val);
      return *this;
    }

    open_chain&
    blank(
      std::string const& val)
    {
      int i = utils::keyword_index(
        io_unit::bl_keywords(), val, "OPEN BLANK");
      u_ptr->blank = static_cast<io_unit::blank_types>(i);
      return *this;
    }

    open_chain&
    status(
      std::string const& val)
    {
      int i = utils::keyword_index(
        io_unit::st_keywords(), val, "OPEN STATUS");
      u_ptr->status = static_cast<io_unit::status_types>(i);
      return *this;
    }

    open_chain&
    iostat(
      int& val)
    {
      iostat_ptr = &val;
      return *this;
    }
  };

  inline
  open_chain
  io::open(
    int unit,
    std::string file,
    bool blank_padding_removed_already)
  {
    if (!blank_padding_removed_already) {
      file = utils::strip_leading_and_trailing_blank_padding(file);
    }
    typedef std::map<int, io_unit>::iterator it;
    it map_iter = units.find(unit);
    if (map_iter != units.end()) {
      if (file.size() == 0) {
        file = map_iter->second.file_name; // f77_std 12.10.1.1
      }
      else if (file == map_iter->second.file_name) {
        throw TBXX_NOT_IMPLEMENTED(); // f77_std 12.10.1.1
      }
      map_iter->second.close();
      units.erase(map_iter);
    }
    map_iter = units.insert(
      std::make_pair(unit, io_unit(unit, file))).first;
    return open_chain(&(map_iter->second));
  }

} // namespace fem

#endif // GUARD
