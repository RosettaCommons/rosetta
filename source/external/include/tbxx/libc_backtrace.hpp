#ifndef TBXX_LIBC_BACKTRACE_HPP
#define TBXX_LIBC_BACKTRACE_HPP

#include <ostream>
#include <cstdlib>
#include <cstring>

#if defined(__GNUC__)
#if defined(__linux) \
 || (defined(__APPLE_CC__) && __APPLE_CC__ >= 5465)
#include <execinfo.h>
#define TBXX_LIBC_BACKTRACE_HAVE_EXECINFO_H
#if ((__GNUC__ > 3) || ((__GNUC__ == 3) && (__GNUC_MINOR__ >= 1))) \
 && !defined(__EDG_VERSION__)
#include <cxxabi.h>
#define TBXX_LIBC_BACKTRACE_HAVE_CXXABI_H
#endif
#endif
#endif

namespace tbxx { namespace libc_backtrace {

  inline
  bool
  show_if_possible(
    std::ostream& ostream,
    int n_frames_skip=0)
  {
    static bool active = false;
    if (active) return false;
    active = true;
    bool result = false;
#if defined(TBXX_LIBC_BACKTRACE_HAVE_EXECINFO_H)
    static const int max_frames = 1024;
    void *array[max_frames];
    int size = backtrace(array, max_frames);
    ostream
      << "libc backtrace ("
      << size - n_frames_skip
      << " frames, most recent call last):"
      << std::endl;
    char **strings = backtrace_symbols(array, size);
    for(int i=size-1;i>=n_frames_skip;i--) {
      char* s = strings[i];
#if defined(TBXX_LIBC_BACKTRACE_HAVE_CXXABI_H)
      const char* m_bgn = 0;
#if defined(__APPLE_CC__)
      if (std::strlen(s) >= 52 && std::strncmp(s+40, "0x", 2) == 0) {
        m_bgn = std::strchr(s+40, ' ');
      }
#else // __linux
      m_bgn = std::strchr(s, '(');
#endif
      if (m_bgn != 0) {
        m_bgn++;
        const char* m_end = std::strchr(m_bgn,
#if defined(__APPLE_CC__)
        ' '
#else // __linux
        '+'
#endif
        );
        long n = m_end - m_bgn;
        if (n > 0) {
          char* mangled = static_cast<char*>(std::malloc(n+1));
          if (mangled != 0) {
            std::strncpy(mangled, m_bgn, n);
            mangled[n] = '\0';
            char* demangled = abi::__cxa_demangle(mangled, 0, 0, 0);
            std::free(mangled);
            if (demangled != 0) {
              long n1 = m_bgn - s;
              long n2 = std::strlen(demangled);
              long n3 = std::strlen(m_end);
              char* b = static_cast<char*>(
                std::malloc(static_cast<size_t>(n1+n2+n3+1)));
              if (b != 0) {
                std::strncpy(b, s, n1);
                std::strncpy(b+n1, demangled, n2);
                std::strncpy(b+n1+n2, m_end, n3);
                b[n1+n2+n3] = '\0';
                s = b;
              }
              std::free(demangled);
            }
          }
        }
      }
#endif // TBXX_LIBC_BACKTRACE_HAVE_CXXABI_H
      ostream << "  " << s << std::endl;
      if (s != strings[i]) std::free(s);
      result = true;
    }
    std::free(strings);
#endif // TBXX_LIBC_BACKTRACE_HAVE_EXECINFO_H
    active = false;
    return result;
  }

}} // tbxx::error_utils

#endif // GUARD
