#ifndef TBXX_PRETTY_TYPE_NAME_HPP
#define TBXX_PRETTY_TYPE_NAME_HPP

#include <typeinfo>
#include <string>
#include <stdexcept>

#if defined(_MS_VER)
// Need do nothing as typeid(..).name() is human-readable with MSVC++

#elif defined(__clang__)
#if __has_include(<cxxabi.h>)
#define TBXX_PRETTY_TYPE_NAME_USE_CXXABI_H
#endif

#elif defined(__GNUC__)
// test copied from libc_backtrace.hpp
#if ((__GNUC__ > 3) || ((__GNUC__ == 3) && (__GNUC_MINOR__ >= 1))) && !defined(__EDG_VERSION__)
#define TBXX_PRETTY_TYPE_NAME_USE_CXXABI_H
#endif

#endif

#ifdef TBXX_PRETTY_TYPE_NAME_USE_CXXABI_H
#include <cxxabi.h>
#endif


namespace tbxx {


  /// A string representing the given type T that is like its declaration in C++ code.
  template <typename T>
  class pretty_type_name : public std::string
  {
  private:
    std::string value() {
      char const *maybe_mangled = typeid(T).name();
#ifdef TBXX_PRETTY_TYPE_NAME_USE_CXXABI_H
      int status;
      char *demangled = abi::__cxa_demangle(maybe_mangled, 0, 0, &status);
      std::string result(status == -2 ? maybe_mangled : demangled);
      std::free(demangled);
#else
      std::string result(maybe_mangled);
#endif
      return result;
    }

  public:
    pretty_type_name()
    : std::string(value())
    {}
  };

  /// A string representing the type of o that is like the declaration of that type in C++ code.
  template <class T>
  pretty_type_name<T> pretty_type_name_of(T const &o) {
    return pretty_type_name<T>();
  }

}
#endif
