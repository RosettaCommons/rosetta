#ifndef FEM_IO_HPP
#define FEM_IO_HPP

#include <fem/io_exceptions.hpp>
#include <fem/utils/misc.hpp>
#include <fem/utils/path.hpp>
#include <fem/utils/random.hpp>
#include <fem/utils/simple_streams.hpp>
#include <fem/utils/string.hpp>
#include <map>
#include <cstdio>

namespace fem {

  static const char file_not_specified[] = "";
  enum unformatted_type { unformatted };

  enum io_modes {
    io_unformatted=0,
    io_list_directed,
    io_formatted
  };

  static const char end_of_unformatted_record = static_cast<char>(0xAAU);

  struct std_file
  {
    std::FILE* ptr;

    std_file(std::FILE* ptr_=0) : ptr(ptr_) {}
  };

  inline
  bool
  is_std_io_unit(
    int unit)
  {
    return (
         unit == 0
      || unit == 5
      || unit == 6);
  }

  struct io_unit
  {
    int number;
    std::string file_name;
    std_file stream;
    bool prev_op_was_write;

    static
    const char**
    ac_keywords()
    {
      static const char* result[] = {"sequential", "direct", 0};
      return result;
    }
    enum access_types { ac_sequential=0, ac_direct, ac_undef };
    access_types access;

    static
    const char**
    fm_keywords()
    {
      static const char* result[] = {"formatted", "unformatted", 0};
      return result;
    }
    enum form_types { fm_formatted=0, fm_unformatted, fm_undef };
    form_types form;

    unsigned recl;

    static
    const char**
    bl_keywords()
    {
      static const char* result[] = {"null", "zero", 0};
      return result;
    }
    enum blank_types { bl_null=0, bl_zero, bl_undef };
    blank_types blank;

    static
    const char**
    st_keywords()
    {
      static const char* result[] = {"old", "new", "scratch", "unknown", 0};
      return result;
    }
    enum status_types { st_old=0, st_new, st_scratch, st_unknown, st_undef};
    status_types status;

    /*! f77_std 12.9.7
        If no error condition or end-of-file condition exists, the
        value of ios is zero. If an error condition exists, the value
        of ios is positive. If an end-of-file condition exists and
        no error condition exists, the value of ios is negative.
        XXX NOT IMPLEMENTED
     */
    int iostat;

    io_unit(
      int number_,
      std::string file_name_=std::string(""),
      std_file stream_=std_file(0))
    :
      number(number_),
      file_name(file_name_),
      stream(stream_),
      prev_op_was_write(false),
      access(ac_undef),
      form(fm_undef),
      recl(0),
      blank(bl_undef),
      status(st_undef),
      iostat(0)
    {}

    std::string
    get_file_name_set_default_if_necessary()
    {
      if (file_name.size() == 0 && !is_std_io_unit(number)) {
        if (status != st_scratch) {
          char buffer[64];
          int actual_width = utils::int_to_string(
            buffer,
            sizeof(buffer),
            /*width*/ 3,
            /*value*/ number,
            /*left_padding_character*/ '0');
          TBXX_ASSERT(actual_width > 0);
          file_name = "fem_io_unit_" + std::string(buffer, actual_width);
        }
        else {
          size_t run_away_counter = 0;
          while (true) {
            file_name = "fem_io_unit_scratch_" + utils::random_name_simple(8);
            if (!utils::path::exists(file_name.c_str())) {
              break;
            }
            run_away_counter++;
            TBXX_ASSERT(run_away_counter < 1000);
          }
        }
      }
      return file_name;
    }

    void
    open(
      int* iostat_ptr)
    {
      if (status == st_undef) {
        status = st_unknown; // f77_std 12.10.1
      }
      if (status == st_scratch) {
        /* f77_std 12.10.1: SCRATCH must not be specified with a named file. */
        /* ifort compile-time warning; file name simply ignored in this
           implementation. */
        file_name.clear();
      }
      if (access == ac_undef) {
        access = ac_sequential; // f77_std 12.10.1
      }
      { // f77_std 12.10.1, end of paragraph under RECL = rl:
        // This specifier [recl] must be given when a file is being
        // connected for direct access; otherwise, it must be omitted.
        if (access == ac_direct) {
          if (recl == 0) {
            std::ostringstream o;
            o << "OPEN error: unit=" << number
              << " connected for DIRECT access but recl=0";
            throw io_err(o.str());
          }
        }
        else {
          if (recl != 0) {
            std::ostringstream o;
            o << "OPEN error: unit=" << number
              << " connected for SEQUENTIAL access but recl=" << recl
              << " (it must be omitted)";
            throw io_err(o.str());
          }
        }
      }
      if (form == fm_undef) {
        if (access == ac_direct) {
          form = fm_unformatted;
        }
        else {
          form = fm_formatted;
        }
      }
      if (blank == bl_undef) {
        blank = bl_null; // f77_std 12.10.1
      }
      bool file_exists = utils::path::exists(
        get_file_name_set_default_if_necessary().c_str());
      if (status == st_old) {
        if (!file_exists) {
          iostat = 1;
          if (iostat_ptr != 0) {
            *iostat_ptr = iostat;
          }
          return;
          std::ostringstream o;
          o << "OPEN OLD error: file does not exist: " << file_name;
          throw io_err(o.str());
        }
      }
      else if (status == st_new) {
        if (file_exists) {
          iostat = 1;
          if (iostat_ptr != 0) {
            *iostat_ptr = iostat;
            return;
          }
          std::ostringstream o;
          o << "OPEN NEW error: file exists already: " << file_name;
          throw io_err(o.str());
        }
      }
      if (access == ac_direct) {
        throw TBXX_NOT_IMPLEMENTED();
      }
      stream.ptr = std::fopen(file_name.c_str(), "ab+");
      if (stream.ptr == 0 || std::fseek(stream.ptr, 0L, SEEK_SET) != 0) {
        if (status == st_new) {
          iostat = 1;
        }
        else {
          stream.ptr = std::fopen(file_name.c_str(), "rb");
          if (stream.ptr == 0 || std::fseek(stream.ptr, 0L, SEEK_SET) != 0) {
            iostat = 1;
          }
        }
        if (iostat == 1) {
          if (iostat_ptr != 0) {
            *iostat_ptr = iostat;
            return;
          }
          throw io_err("Error opening file: " + file_name);
        }
      }
      if (status == st_new) {
        /* f77_std 12.10.1: Successful execution of an OPEN statement with
           NEW specified creates the file and changes the status to OLD. */
        status = st_old;
      }
    }

    void
    close(
      int* iostat_ptr=0,
      bool status_delete=false)
    {
      if (iostat_ptr != 0) *iostat_ptr = 0; // XXX
      if (is_std_io_unit(number)) return;
      if (stream.ptr != 0) {
        std::fclose(stream.ptr);
        stream.ptr = 0;
      }
      if (status == st_scratch || status_delete) {
        std::remove(file_name.c_str());
      }
    }

    void
    backspace(
      int* iostat_ptr)
    {
      throw TBXX_NOT_IMPLEMENTED();
    }

    void
    endfile(
      int* iostat_ptr)
    {
      if (is_std_io_unit(number)) {
        throw TBXX_NOT_IMPLEMENTED();
      }
      if (!utils::path::truncate_file_at_current_position(stream.ptr)) {
        throw io_err("ENDFILE failure: " + file_name);
      }
      prev_op_was_write = false;
    }

    void
    rewind(
      int* iostat_ptr)
    {
      if (stream.ptr == 0 || std::fseek(stream.ptr, 0L, SEEK_SET) != 0) {
        iostat = 1;
        if (iostat_ptr != 0) {
          *iostat_ptr = iostat;
          return;
        }
        throw io_err("Error rewinding file: " + file_name);
      }
      prev_op_was_write = false;
    }
  };

  struct io : utils::noncopyable
  {
    std::map<int, io_unit> units;

    io()
    {
      units.insert(std::make_pair(0, io_unit(0, "", stderr)));
      units.insert(std::make_pair(5, io_unit(5, "", stdin)));
      units.insert(std::make_pair(6, io_unit(6, "", stdout)));
    }

    ~io()
    {
      typedef std::map<int, io_unit>::iterator it;
      it e = units.end();
      for(it i=units.begin();i!=e;i++) {
        i->second.close();
      }
    }

    io_unit*
    unit_ptr(
      int unit,
      bool auto_open=false)
    {
      typedef std::map<int, io_unit>::iterator it;
      it map_iter = units.find(unit);
      if (map_iter == units.end()) {
        if (!auto_open) return 0;
        map_iter = units.insert(
          std::make_pair(unit, io_unit(unit))).first;
        map_iter->second.open(/*iostat_ptr*/ 0);
      }
      return &(map_iter->second);
    }

    //! Easy C++ access.
    std::string
    file_name_of_unit(
      int unit)
    {
      io_unit* u_ptr = unit_ptr(unit);
      if (u_ptr == 0) return "";
      return u_ptr->file_name;
    }

    inline
    struct inquire_chain
    inquire_unit(
      int unit);

    inline
    struct inquire_chain
    inquire_file(
      std::string file,
      bool blank_padding_removed_already=false);

    inline
    struct open_chain
    open(
      int unit,
      std::string file=std::string(),
      bool blank_padding_removed_already=false);

    inline
    struct close_chain
    close(
      int unit);

    inline
    struct file_positioning_chain
    backspace(
      int unit);

    inline
    struct file_positioning_chain
    endfile(
      int unit);

    inline
    struct file_positioning_chain
    rewind(
      int unit);

    inline
    bool
    is_opened_simple(
      std::string const& file_name) const
    {
      typedef std::map<int, io_unit>::const_iterator it;
      it e = units.end();
      for(it i=units.begin();i!=e;i++) {
        if (i->second.file_name == file_name) {
          return true;
        }
      }
      return false;
    }

    utils::slick_ptr<utils::simple_ostream>
    simple_ostream(
      int unit)
    {
      io_unit* u_ptr = unit_ptr(unit, /*auto_open*/ true);
      std_file& sf = u_ptr->stream;
      if (!u_ptr->prev_op_was_write) {
        if (!is_std_io_unit(unit)) {
          if (!utils::path::truncate_file_at_current_position(sf.ptr)) {
            throw io_err(
              "Cannot truncate file for writing: " + u_ptr->file_name);
          }
        }
        u_ptr->prev_op_was_write = true;
      }
      return utils::slick_ptr<utils::simple_ostream>(new
        utils::simple_ostream_to_c_file(sf.ptr));
    }

    utils::slick_ptr<utils::simple_istream>
    simple_istream(
      int unit)
    {
      io_unit* u_ptr = unit_ptr(unit, /*auto_open*/ true);
      u_ptr->prev_op_was_write = false;
      return utils::slick_ptr<utils::simple_istream>(new
        utils::simple_istream_from_c_file(u_ptr->stream.ptr));
    }

    void
    flush(
      int unit)
    {
      io_unit* u_ptr = unit_ptr(unit);
      if (u_ptr != 0) {
        int status = std::fflush(u_ptr->stream.ptr);
        if (status != 0) {
          throw io_err(
            "Error flushing file: " + u_ptr->file_name);
        }
      }
    }
  };

} // namespace fem

#endif // GUARD
