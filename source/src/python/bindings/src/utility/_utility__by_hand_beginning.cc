// (c) Copyright Rosetta Commons Member Institutions.
#include "utility/vector1.hh"
#include "utility/pointer/access_ptr.hh"

#include <utility/excn/Exceptions.hh>
#include <utility/stream_util.hh>
#include "utility/exit.hh"

#include <platform/types.hh>

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>

#include <iostream>
#include <ostream>
#include <istream>
#include <sstream>
#include <set>
#include <map>



// Includes for dummy bindings to simplify import orders
#include <utility/inline_file_provider.hh>


namespace bp = boost::python;


#ifndef BOOST_ADAPTBX_PYTHON_STREAMBUF_H
#define BOOST_ADAPTBX_PYTHON_STREAMBUF_H

/*
*** License agreement ***
cctbx Copyright (c) 2006, The Regents of the University of
California, through Lawrence Berkeley National Laboratory (subject to
receipt of any required approvals from the U.S. Dept. of Energy).  All
rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

(2) Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

(3) Neither the name of the University of California, Lawrence Berkeley
National Laboratory, U.S. Dept. of Energy nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes,
patches, or upgrades to the features, functionality or performance of
the source code ("Enhancements") to anyone; however, if you choose to
make your Enhancements available either publicly, or directly to
Lawrence Berkeley National Laboratory, without imposing a separate
written license agreement for such Enhancements, then you hereby grant
the following license: a  non-exclusive, royalty-free perpetual license
to install, use, modify, prepare derivative works, incorporate into
other computer software, distribute, and sublicense such enhancements or
derivative works thereof, in binary and source code form.
*/

#include <boost/python/object.hpp>
#include <boost/python/str.hpp>
#include <boost/python/extract.hpp>

#include <boost/optional.hpp>
#include <boost/utility/typed_in_place_factory.hpp>

#include <streambuf>
#include <iostream>

namespace boost_adaptbx { namespace python {

namespace bp = boost::python;

/// A stream buffer getting data from and putting data into a Python file object
/** The aims are as follow:

    - Given a C++ function acting on a standard stream, e.g.

      \code
      void read_inputs(std::istream& input) {
        ...
        input >> something >> something_else;
      }
      \endcode

      and given a piece of Python code which creates a file-like object,
      to be able to pass this file object to that C++ function, e.g.

      \code
      import gzip
      gzip_file_obj = gzip.GzipFile(...)
      read_inputs(gzip_file_obj)
      \endcode

      and have the standard stream pull data from and put data into the Python
      file object.

    - When Python \c read_inputs() returns, the Python object is able to
      continue reading or writing where the C++ code left off.

    - Operations in C++ on mere files should be competitively fast compared
      to the direct use of \c std::fstream.


    \b Motivation

      - the standard Python library offer of file-like objects (files,
        compressed files and archives, network, ...) is far superior to the
        offer of streams in the C++ standard library and Boost C++ libraries.

      - i/o code involves a fair amount of text processing which is more
        efficiently prototyped in Python but then one may need to rewrite
        a time-critical part in C++, in as seamless a manner as possible.

    \b Usage

      - the Python side

        \code
          from rosetta.utility import istream, ostream
          ...
          read_inputs(istream(in_file_obj))
          ...
          write_outputs(ostream(out_file_obj))
        \endcode

  Note: references are to the C++ standard (the numbers between parentheses
  at the end of references are margin markers).
*/
class streambuf : public std::basic_streambuf<char>
{
  private:
    typedef std::basic_streambuf<char> base_t;

  public:
    /* The syntax
        using base_t::char_type;
       would be nicer but Visual Studio C++ 8 chokes on it
    */
    typedef base_t::char_type   char_type;
    typedef base_t::int_type    int_type;
    typedef base_t::pos_type    pos_type;
    typedef base_t::off_type    off_type;
    typedef base_t::traits_type traits_type;

    // work around Visual C++ 7.1 problem
    inline static int
    traits_type_eof() { return traits_type::eof(); }

    /// The default size of the read and write buffer.
    /** They are respectively used to buffer data read from and data written to
        the Python file object. It can be modified from Python.
    */
    static platform::Size default_buffer_size;

    /// Construct from a Python file object
    /** if buffer_size is 0 the current default_buffer_size is used.
    */
    streambuf(
      bp::object& python_file_obj,
      platform::Size buffer_size_=0)
    :
      py_read (getattr(python_file_obj, "read",  bp::object())),
      py_write(getattr(python_file_obj, "write", bp::object())),
      py_seek (getattr(python_file_obj, "seek",  bp::object())),
      py_tell (getattr(python_file_obj, "tell",  bp::object())),
      buffer_size(buffer_size_ != 0 ? buffer_size_ : default_buffer_size),
      write_buffer(0),
      pos_of_read_buffer_end_in_py_file(0),
      pos_of_write_buffer_end_in_py_file(buffer_size),
      farthest_pptr(0)
    {
      /* Some Python file objects (e.g. sys.stdout and sys.stdin)
         have non-functional seek and tell. If so, assign None to
         py_tell and py_seek.
       */
      if (py_tell != bp::object()) {
        try {
          py_tell();
        }
        catch (bp::error_already_set&) {
          py_tell = bp::object();
          py_seek = bp::object();
          /* Boost.Python does not do any Python exception handling whatsoever
             So we need to catch it by hand like so.
           */
          PyErr_Clear();
        }
      }

      if (py_write != bp::object()) {
        // C-like string to make debugging easier
        write_buffer = new char[buffer_size + 1];
        write_buffer[buffer_size] = '\0';
        setp(write_buffer, write_buffer + buffer_size);  // 27.5.2.4.5 (5)
        farthest_pptr = pptr();
      }
      else {
        // The first attempt at output will result in a call to overflow
        setp(0, 0);
      }

      if (py_tell != bp::object()) {
        off_type py_pos = bp::extract<off_type>(py_tell());
        pos_of_read_buffer_end_in_py_file = py_pos;
        pos_of_write_buffer_end_in_py_file = py_pos;
      }
    }

    /// Mundane destructor freeing the allocated resources
    virtual ~streambuf() {
      if (write_buffer) delete[] write_buffer;
    }

    /// C.f. C++ standard section 27.5.2.4.3
    /** It is essential to override this virtual function for the stream
        member function readsome to work correctly (c.f. 27.6.1.3, alinea 30)
     */
    virtual std::streamsize showmanyc() {
      int_type const failure = traits_type::eof();
      int_type status = underflow();
      if (status == failure) return -1;
      return egptr() - gptr();
    }

    /// C.f. C++ standard section 27.5.2.4.3
    virtual int_type underflow() {
      int_type const failure = traits_type::eof();
      if (py_read == bp::object()) {
        throw std::invalid_argument(
          "That Python file object has no 'read' attribute");
      }
      read_buffer = py_read(buffer_size);
      char *read_buffer_data;
      bp::ssize_t py_n_read;
      if (PyString_AsStringAndSize(read_buffer.ptr(),
                                   &read_buffer_data, &py_n_read) == -1) {
        setg(0, 0, 0);
        throw std::invalid_argument(
          "The method 'read' of the Python file object "
          "did not return a string.");
      }
      off_type n_read = (off_type)py_n_read;
      pos_of_read_buffer_end_in_py_file += n_read;
      setg(read_buffer_data, read_buffer_data, read_buffer_data + n_read);
      // ^^^27.5.2.3.1 (4)
      if (n_read == 0) return failure;
      return traits_type::to_int_type(read_buffer_data[0]);
    }

    /// C.f. C++ standard section 27.5.2.4.5
    virtual int_type overflow(int_type c=traits_type_eof()) {
      if (py_write == bp::object()) {
        throw std::invalid_argument(
          "That Python file object has no 'write' attribute");
      }
      farthest_pptr = std::max(farthest_pptr, pptr());
      off_type n_written = (off_type)(farthest_pptr - pbase());
      bp::str chunk(pbase(), farthest_pptr);
      py_write(chunk);
      if (!traits_type::eq_int_type(c, traits_type::eof())) {
        py_write(traits_type::to_char_type(c));
        n_written++;
      }
      if (n_written) {
        pos_of_write_buffer_end_in_py_file += n_written;
        setp(pbase(), epptr());
        // ^^^ 27.5.2.4.5 (5)
        farthest_pptr = pptr();
      }
      return traits_type::eq_int_type(
        c, traits_type::eof()) ? traits_type::not_eof(c) : c;
    }

    /// Update the python file to reflect the state of this stream buffer
    /** Empty the write buffer into the Python file object and set the seek
        position of the latter accordingly (C++ standard section 27.5.2.4.2).
        If there is no write buffer or it is empty, but there is a non-empty
        read buffer, set the Python file object seek position to the
        seek position in that read buffer.
    */
    virtual int sync() {
      int result = 0;
      farthest_pptr = std::max(farthest_pptr, pptr());
      if (farthest_pptr && farthest_pptr > pbase()) {
        off_type delta = pptr() - farthest_pptr;
        int_type status = overflow();
        if (traits_type::eq_int_type(status, traits_type::eof())) result = -1;
        if (py_seek != bp::object()) py_seek(delta, 1);
      }
      else if (gptr() && gptr() < egptr()) {
        if (py_seek != bp::object()) py_seek(gptr() - egptr(), 1);
      }
      return result;
    }

    /// C.f. C++ standard section 27.5.2.4.2
    /** This implementation is optimised to look whether the position is within
        the buffers, so as to avoid calling Python seek or tell. It is
        important for many applications that the overhead of calling into Python
        is avoided as much as possible (e.g. parsers which may do a lot of
        backtracking)
    */
    virtual
    pos_type seekoff(off_type off, std::ios_base::seekdir way,
                     std::ios_base::openmode which=  std::ios_base::in
                                                   | std::ios_base::out)
    {
      /* In practice, "which" is either std::ios_base::in or out
         since we end up here because either seekp or seekg was called
         on the stream using this buffer. That simplifies the code
         in a few places.
      */
      int const failure = off_type(-1);

      if (py_seek == bp::object()) {
        throw std::invalid_argument(
          "That Python file object has no 'seek' attribute");
      }

      // we need the read buffer to contain something!
      if (which == std::ios_base::in && !gptr()) {
        if (traits_type::eq_int_type(underflow(), traits_type::eof())) {
          return failure;
        }
      }

      // compute the whence parameter for Python seek
      int whence;
      switch (way) {
        case std::ios_base::beg:
          whence = 0;
          break;
        case std::ios_base::cur:
          whence = 1;
          break;
        case std::ios_base::end:
          whence = 2;
          break;
        default:
          return failure;
      }

      // Let's have a go
      boost::optional<off_type> result = seekoff_without_calling_python(
        off, way, which);
      if (!result) {
        // we need to call Python
        if (which == std::ios_base::out) overflow();
        if (way == std::ios_base::cur) {
          if      (which == std::ios_base::in)  off -= egptr() - gptr();
          else if (which == std::ios_base::out) off += pptr() - pbase();
        }
        py_seek(off, whence);
        result = off_type(bp::extract<off_type>(py_tell()));
        if (which == std::ios_base::in) underflow();
      }
      return *result;
    }

    /// C.f. C++ standard section 27.5.2.4.2
    virtual
    pos_type seekpos(pos_type sp,
                     std::ios_base::openmode which=  std::ios_base::in
                                                   | std::ios_base::out)
    {
      return streambuf::seekoff(sp, std::ios_base::beg, which);
    }

  private:
    bp::object py_read, py_write, py_seek, py_tell;

    platform::Size buffer_size;

    /* This is actually a Python string and the actual read buffer is
       its internal data, i.e. an array of characters. We use a Boost.Python
       object so as to hold on it: as a result, the actual buffer can't
       go away.
    */
    bp::object read_buffer;

    /* A mere array of char's allocated on the heap at construction time and
       de-allocated only at destruction time.
    */
    char *write_buffer;

    off_type pos_of_read_buffer_end_in_py_file,
             pos_of_write_buffer_end_in_py_file;

    // the farthest place the buffer has been written into
    char *farthest_pptr;


    boost::optional<off_type> seekoff_without_calling_python(
      off_type off,
      std::ios_base::seekdir way,
      std::ios_base::openmode which)
    {
      boost::optional<off_type> const failure;

      // Buffer range and current position
      off_type buf_begin, buf_end, buf_cur, upper_bound;
      off_type pos_of_buffer_end_in_py_file;
      if (which == std::ios_base::in) {
        pos_of_buffer_end_in_py_file = pos_of_read_buffer_end_in_py_file;
        buf_begin = reinterpret_cast<std::streamsize>(eback());
        buf_cur = reinterpret_cast<std::streamsize>(gptr());
        buf_end = reinterpret_cast<std::streamsize>(egptr());
        upper_bound = buf_end;
      }
      else if (which == std::ios_base::out) {
        pos_of_buffer_end_in_py_file = pos_of_write_buffer_end_in_py_file;
        buf_begin = reinterpret_cast<std::streamsize>(pbase());
        buf_cur = reinterpret_cast<std::streamsize>(pptr());
        buf_end = reinterpret_cast<std::streamsize>(epptr());
        farthest_pptr = std::max(farthest_pptr, pptr());
        upper_bound = reinterpret_cast<std::streamsize>(farthest_pptr) + 1;
      }

      // Sought position in "buffer coordinate"
      off_type buf_sought;
      if (way == std::ios_base::cur) {
        buf_sought = buf_cur + off;
      }
      else if (way == std::ios_base::beg) {
        buf_sought = buf_end + (off - pos_of_buffer_end_in_py_file);
      }
      else if (way == std::ios_base::end) {
        return failure;
      }

      // if the sought position is not in the buffer, give up
      if (buf_sought < buf_begin || buf_sought >= upper_bound) return failure;

      // we are in wonderland
      if      (which == std::ios_base::in)  gbump(buf_sought - buf_cur);
      else if (which == std::ios_base::out) pbump(buf_sought - buf_cur);
      return pos_of_buffer_end_in_py_file + (buf_sought - buf_end);
    }

  public:

    class istream : public std::istream
    {
      public:
        istream(streambuf& buf) : std::istream(&buf)
        {
          exceptions(std::ios_base::badbit);
        }

        ~istream() { if (this->good()) this->sync(); }
    };

    class ostream : public std::ostream
    {
      public:
        ostream(streambuf& buf) : std::ostream(&buf)
        {
          exceptions(std::ios_base::badbit);
        }

        ~ostream() { if (this->good()) this->flush(); }
    };
};

platform::Size streambuf::default_buffer_size = 1024;

struct streambuf_capsule
{
  streambuf python_streambuf;

  streambuf_capsule(
    bp::object& python_file_obj,
    platform::Size buffer_size=0)
  :
    python_streambuf(python_file_obj, buffer_size)
  {}
};

struct ostream : private streambuf_capsule, streambuf::ostream
{
  ostream(
    bp::object& python_file_obj,
    platform::Size buffer_size = 0)
  :
    streambuf_capsule(python_file_obj, buffer_size),
    streambuf::ostream(python_streambuf)
  {}

  ~ostream()
  {
    try {
      if (this->good()) this->flush();
    }
    catch (bp::error_already_set&) {
      PyErr_Clear();
      throw std::runtime_error(
        "Problem closing python ostream.\n"
        "  Known limitation: the error is unrecoverable. Sorry.\n"
        "  Suggestion for programmer: add ostream.flush() before"
        " returning.");
    }
  }
};

struct istream : private streambuf_capsule, streambuf::istream
{
  istream(
    bp::object& python_file_obj,
    platform::Size buffer_size = 0)
  :
    streambuf_capsule(python_file_obj, buffer_size),
    streambuf::istream(python_streambuf)
  {}

  ~istream()
  {
    try {
      if (this->good()) this->sync();
    }
    catch (bp::error_already_set&) {
      PyErr_Clear();
      throw std::runtime_error(
        "Problem closing python istream.\n"
        "  Known limitation: the error is unrecoverable. Sorry.\n"
        "  Suggestion for programmer: add istream.sync() before"
        " returning.");
    }
  }
};

struct python_ostream_wrapper
{
  typedef ostream wt;

  static void
  wrap()
  {
    using namespace boost::python;
    class_<std::ostream, boost::noncopyable>("std_ostream", no_init);
    class_<wt, boost::noncopyable, bases<std::ostream> >("ostream", "Buffered ostream wrapper for file objects. Buffer is flushed on deletion of the wrapper object.", no_init)
      .def( init<object&>(args("file"), "Initialize an ostream wrapper with the default buffer size."))
      .def( init<object&, platform::Size>(args("file", "buffer_size"), "Initialize an ostream wrapper with the given buffer size in bytes."));
  }
};

struct python_istream_wrapper
{
  typedef istream wt;

  static void
  wrap()
  {
    using namespace boost::python;
    class_<std::istream, boost::noncopyable>("std_istream", no_init);
    class_<wt, boost::noncopyable, bases<std::istream> >("istream", "Buffered istream wrapper for file objects. Buffer is flushed on deletion of the wrapper object.", no_init)
      .def( init<object&>(args("file"), "Initialize an istream wrapper with the default buffer size."))
      .def( init<object&, platform::Size>(args("file", "buffer_size"), "Initialize an istream wrapper with the given buffer size in bytes."));
  }
};



}} // boost_adaptbx::python

#endif // GUARD



// template< class T >
// T * getCAP( utility::pointer::access_ptr<T> rs ) {
//   T & rs_ref( *rs );
//   T * rs_ptr = &rs_ref;
//   return rs_ptr;
// }

// std::pair ---------------------------------------------------------------------------------------------------
//

//vector1 wrapper requires ostream '<<' operator

template< class T1, class T2 >
std::ostream& operator<<(std::ostream& strm, const std::pair< T1, T2>& kvPair)
{
  strm << "(" << kvPair.first << ", " << kvPair.second << ")";
  return strm;
}

template <class T1, class T2>
std::string pair_repr(std::pair<T1, T2> const & v)
{
    std::ostringstream os;

    os << v;

    return os.str();
}

template< class T1, class T2 >
void wrap_std_pair(std::string name)
{
    bp::class_< std::pair< T1, T2 > >(name.c_str())
			.def( bp::init< T1 const &, T2 const & >(( bp::arg("__a"),
bp::arg("__b") )))
			.def_readwrite( "first", &std::pair< T1, T2 >::first)
			.def_readwrite( "second", &std::pair< T1, T2 >::second)
      .def("__str__", &pair_repr<T1, T2> );
}

/*
template <class T>
std::ostream& operator <<(std::ostream &os, utility::vector1<T> const & v)
{
    os << "[";
    for(unsigned int i=1; i<=v.size(); i++) {
        os << v[i] << ", ";
    }
    os << "]";
    return os;
} */

template <class T>
std::string vector1_repr(utility::vector1<T> const & v)
{
    std::ostringstream os;

    os << "[";
    for(unsigned int i=1; i<=v.size(); i++) {
        os << v[i] << ", ";
    }
    os << "]";
    return os.str();
}

template< class TT > inline void vector1_set( utility::vector1<TT> & v, platform::Size const & i, TT const & val ) { v[i] = val; }
template< class TT > inline platform::Size vector1_len( utility::vector1<TT> & v ) { return v.size(); }

template< class TT > inline std::string vector1_str( utility::vector1<TT> & v ) { std::ostringstream s; s<<v; return s.str(); }

template< class TT > inline typename utility::vector1<TT>::iterator vector1_begin( utility::vector1<TT> & v ) { return v.begin(); }
template< class TT > inline typename utility::vector1<TT>::iterator vector1_end  ( utility::vector1<TT> & v ) { return v.end(); }

template< class TT > inline void vector1_reserve( utility::vector1<TT> & v, platform::Size n) { v.reserve(n); }
template< class TT > inline void vector1_resize( utility::vector1<TT> & v, platform::Size n) { v.resize(n); }

template< class Htype, class CP, class CP_const>
void wrap_vector1(std::string name) {
  typedef utility::vector1<Htype> Ttype;
  typedef utility::vectorL<1,Htype, std::allocator<Htype> > Btype;
  typedef std::vector<Htype> Vtype;
  bp::class_<Ttype>(name.c_str())
    .def( bp::init< platform::Size >() )
    .def( bp::init< utility::vector1<Htype> const & >() )
    // .def( bp::init< platform::Size, TT >() )
    .def("__getitem__"
        , (Htype const & (Ttype::*)(platform::Size const) const)( &Ttype::at )
        , CP_const()    )
    .def("__getitem__"
        , (Htype & (Ttype::*)(platform::Size const))( &Ttype::at )
        , CP()        )
    .def("__setitem__"
        , &vector1_set<Htype> )
    .def("append"
        , (Btype & (Btype::*)(Htype const &))( &Btype::add_back )
        , bp::return_value_policy< bp::reference_existing_object >()        )
    .def("__len__", & vector1_len<Htype> )
    .def("__iter__", bp::range(&vector1_begin<Htype>,&vector1_end<Htype>) )

    //.def("__str__", & vector1_str<Htype> )
    .def("__str__", & vector1_repr<Htype> )
    //.def( bp::self_ns::str( bp::self ) )

    .def("reserve", &vector1_reserve<Htype> )
    .def("resize", &vector1_resize<Htype> )

  ;
}

template< class Htype, class CP, class CP_const>
void wrap_vector1_part(const char * name) {
  typedef utility::vector1<Htype> Ttype;
  typedef utility::vectorL<1,Htype, std::allocator<Htype> > Btype;
  typedef std::vector<Htype> Vtype;
  bp::class_<Ttype>(name)
    .def("__getitem__"
        , (Htype const & (Ttype::*)(platform::Size const) const)( &Ttype::at )
        , CP_const()    )
    .def("__getitem__"
        , (Htype & (Ttype::*)(platform::Size const))( &Ttype::at )
        , CP()        )
    .def("__setitem__"
        , &vector1_set<Htype> )
    .def("append"
        , (Btype & (Btype::*)(Htype const &))( &Btype::add_back )
        , bp::return_value_policy< bp::reference_existing_object >()        )
    .def("__len__", & vector1_len<Htype> )
    .def("__iter__", bp::range(&vector1_begin<Htype>,&vector1_end<Htype>) )
  ;
}

// std::vector -----------------------------------------------------------------------------------------------------------------
template <class T>
std::string std_vector_repr(std::vector<T> const & v)
{
    std::ostringstream os;

    os << "[";
    for(unsigned int i=1; i<=v.size(); i++) {
        os << v[i] << ", ";
    }
    os << "]";
    return os.str();
}

template< class TT > inline void std_vector_set( std::vector<TT> & v, platform::Size const & i, TT const & val ) { v[i] = val; }
template< class TT > inline platform::Size std_vector_len( std::vector<TT> & v ) { return v.size(); }
template< class TT > inline TT& std_vector_at(std::vector<TT> & v, platform::Size i ) { return v.at(i); }

template< class TT > inline std::string std_vector_str( std::vector<TT> & v ) { std::ostringstream s; s<<v; return s.str(); }

template< class TT > inline typename std::vector<TT>::iterator std_vector_begin( std::vector<TT> & v ) { return v.begin(); }
template< class TT > inline typename std::vector<TT>::iterator std_vector_end  ( std::vector<TT> & v ) { return v.end(); }

template< class TT > inline void std_vector_reserve( std::vector<TT> & v, platform::Size n) { v.reserve(n); }
template< class TT > inline void std_vector_resize( std::vector<TT> & v, platform::Size n) { v.resize(n); }

template< class Htype, class CP, class CP_const>
void wrap_std_vector(std::string name) {
  typedef std::vector<Htype> Ttype;
  typedef utility::vectorL<1,Htype, std::allocator<Htype> > Btype;
  bp::class_<Ttype>(name.c_str())
    .def( bp::init< platform::Size >() )
	  .def( bp::init< std::vector<Htype> const & >() )
    // .def( bp::init< platform::Size, TT >() )


    // .def("__getitem__"
    //     , (Htype const & (Ttype::*)(platform::Size const) const)( &Ttype::at )
    //     , CP_const()    )
    // .def("__getitem__"
    //     ,  &std_vector_at<Htype>
    //     , CP_const()        )
    .def("__getitem__"
        ,  &std_vector_at<Htype>
        , CP()        )
    .def("__setitem__"
        , &std_vector_set<Htype> )
    .def("append"
        , (Btype & (Btype::*)(Htype const &))( &Btype::add_back )
        , bp::return_value_policy< bp::reference_existing_object >()        )


    .def("__len__", & std_vector_len<Htype> )
    .def("__iter__", bp::range(&std_vector_begin<Htype>,&std_vector_end<Htype>) )

    //.def("__str__", & vector1_str<Htype> )
    .def("__str__", & std_vector_repr<Htype> )
    //.def( bp::self_ns::str( bp::self ) )

    .def("reserve", &std_vector_reserve<Htype> )
    .def("resize", &std_vector_resize<Htype> )

  ;
}


// std::map --------------------------------------------------------------------------------------------------------------------

template< class Key, class Val >
struct map_wrapper
{
  typedef std::map<Key,Val> Map;

  static boost::python::list keys(Map const& self)
  {
    boost::python::list t;

    for(typename Map::const_iterator it = self.begin(); it!=self.end(); ++it)
    {
      t.append(it->first);
    }

    return t;
  }

  static boost::python::list values(Map const& self)
  {
    boost::python::list t;

    for(typename Map::const_iterator it = self.begin(); it!=self.end(); ++it)
    {
      t.append(it->second);
    }

    return t;
  }

  static boost::python::list items(Map const& self)
  {
    boost::python::list t;

    for(typename Map::const_iterator it = self.begin(); it!=self.end(); ++it)
    {
        t.append( boost::python::make_tuple(it->first, it->second) );
    }

    return t;
  }
};

template< class Key, class Val >
void wrap_std_map(std::string name)
{
	bp::class_< std::map< Key, Val > >(name.c_str())
		.def( bp::map_indexing_suite< std::map< Key, Val > >())
    .def( "keys", &map_wrapper<Key, Val>::keys )
    .def( "values", &map_wrapper<Key, Val>::values )
    .def( "items", &map_wrapper<Key, Val>::items )
    ;
}

// std::set --------------------------------------------------------------------------------------------------------------------
template< class T > void add_to_set(std::set<T> & s, const T & v) { s.insert(v); }
template< class T > void erase_from_set(std::set<T> & s, const T & v) { s.erase(v); }

template< class TT > inline typename std::set<TT>::iterator set_begin( std::set<TT> & v ) { return v.begin(); }
template< class TT > inline typename std::set<TT>::iterator set_end  ( std::set<TT> & v ) { return v.end(); }

template <class T>
std::string set_repr(std::set<T> const & s)
{
	typedef std::set<T> Stype;
	typedef typename std::set<T>::iterator Stype_iterator;

    std::ostringstream os;
    os << "<set>[";
    for(Stype_iterator p=s.begin(); p!=s.end(); ++p) os << *p << ", ";
    os << "]";
    return os.str();
}


template< class Htype, class CP, class CP_const>
void wrap_std_set(std::string name)
{
	typedef std::set<Htype> Ttype;
	bp::class_<Ttype>(name.c_str())
	.def( bp::init< >() )
	.def( bp::init< std::set<Htype> const & >() )

	.def("__contains__", &std::set<Htype>::count )
	.def("add", &add_to_set<Htype> )
	.def("erase",  &erase_from_set<Htype> )
	.def("__len__",  &std::set<Htype>::size )
	.def("__iter__", bp::range(&set_begin<Htype>, &set_end<Htype> ) )
	.def("__str__", &set_repr<Htype> )
	;
}

template< class Type >
void wrap_owning_pointer(char * name)
{
    bp::class_<Type>(name)
        //.def("get", &Type::get)
        .def("reset_to_null", &Type::reset_to_null)
    ;
}


/*
template< class T >
void wrap_access_pointer(std::string class_name)
{
    boost::python::implicitly_convertible< utility::pointer::access_ptr< T >
                                         , utility::pointer::access_ptr< T const > >();

    bp::class_< utility::pointer::access_ptr< T > >( std::string(class_name+"AP").c_str() )
        .def("get", (  T * (*)( utility::pointer::access_ptr<T> )  )( & wrap_access_pointer_get_function<T> )
             , bp::return_value_policy< bp::reference_existing_object >() );

    bp::class_< utility::pointer::access_ptr< T const > >( std::string(class_name+"CAP").c_str() )
        .def("get", (  T const * (*)( utility::pointer::access_ptr<T const > )  )( & wrap_access_pointer_get_function<T const> )
             , bp::return_value_policy< bp::reference_existing_object >() );
}
*/

// .def("__iter__", bp::range( &core::pose::Pose::res_begin, &core::pose::Pose::res_end));

inline bool vector1_bool_get ( utility::vector1<bool> & v, int i ) { if(v[i]) return true; else return false; }
inline void vector1_bool_push( utility::vector1<bool> & v, bool h ) { return v.push_back(h); }

template <class T>
void expose_basic_type(std::string name)
{
  typedef bp::return_value_policy< bp::reference_existing_object > CP_REF;
  typedef bp::return_value_policy< bp::copy_const_reference >      CP_CCR;
  typedef bp::return_value_policy< bp::copy_non_const_reference >  CP_CNCR;

  wrap_std_vector< T, CP_CNCR, CP_CCR >("std_vector_" + name);

  wrap_vector1< T, CP_CNCR, CP_CCR >("vector1_" + name);
  wrap_vector1< utility::vector1<T>, CP_REF, CP_REF >("vec1_vec1_" + name);
  wrap_std_set< T, CP_CNCR, CP_CCR >("set_" + name);
}





void pyexit_callback(void)
{
    //throw "RosettaException";
	throw utility::excn::EXCN_Msg_Exception("PyExitCallbackException");
}

void set_pyexit_callback(void)
{
	utility::set_main_exit_callback(pyexit_callback);
}

void py_xinc_ref(PyObject *o) { Py_XINCREF(o); };
void py_xdec_ref(PyObject *o) { Py_XDECREF(o); };

// Static int test for Windows DLL's
static int __static_int_ = 42;
static std::string __static_string_("static string...");
int __static_int_test() { 	return __static_int_; }
std::string __static_string_test() { 	return __static_string_; }

void __utility_by_hand_beginning__()
{
	// #ifdef DEBUG
	// 	std::cerr << "__utility_by_hand_beginning__..." << std::endl;
	// #endif
    bp::def("__static_int_test", __static_int_test);
    bp::def("__static_string_test", __static_string_test);

    bp::def("set_pyexit_callback", set_pyexit_callback);

    bp::def("py_xinc_ref", py_xinc_ref);
    bp::def("py_xdec_ref", py_xdec_ref);


    // //bp::def("set_main_exit_callback", set_main_exit_callback);
    // bp::def("pyexit_callback", pyexit_callback);

    //wrap_owning_pointer<core::pack::task::PackerTaskOP>("PackerTaskOP");

    // istream and ostream converters for python file objects
    // See documentation in class definition for usage.
    boost_adaptbx::python::python_ostream_wrapper::wrap();
    boost_adaptbx::python::python_istream_wrapper::wrap();

    // std::stringstream wrappers ---------------------------------------------------------------------
    bp::class_< std::ostream, boost::noncopyable >("OStream", bp::no_init);
    bp::class_< std::iostream, boost::noncopyable >("IOStream", bp::no_init);

    typedef void ( std::ostringstream::*ostringstream_str_set_function_type )( std::string const & );
    typedef std::string ( std::ostringstream::*ostringstream_str_get_function_type )( ) const;

    bp::class_< std::ostringstream, bp::bases<std::ostream>, boost::noncopyable >("OStringStream")
        .def("str", ostringstream_str_set_function_type( &::std::ostringstream::str ) )
        .def("str", ostringstream_str_get_function_type( &::std::ostringstream::str ) )
    ;

    typedef void ( std::stringstream::*stringstream_str_set_function_type )( std::string const & );
    typedef std::string ( std::stringstream::*stringstream_str_get_function_type )( ) const;

    bp::class_< std::stringstream, bp::bases<std::iostream>, boost::noncopyable >("StringStream")
        .def( bp::init<std::string>(bp::args("str"), "Initialize stringstream with the given string."))
        .def("str", stringstream_str_set_function_type( &::std::stringstream::str ) )
        .def("str", stringstream_str_get_function_type( &::std::stringstream::str ) );


	// some holders, todo: maybe add some values/functions later
	// bp::class_< std::ios_base::openmode, boost::noncopyable >("std_ios_base_openmode");
	// bp::class_< std::_Ios_Openmode, boost::noncopyable >("std__Ios_Openmode");


    utility::wrap_access_pointer< utility::vector1< bool > >("utility_vector1_bool");


	// Dummy bindings to simplify import orders
	boost::python::class_< utility::SingletonBase<utility::Inline_File_Provider>, boost::noncopyable >( "__utility_SingletonBase_utility_Inline_File_Provider__");
}
