#ifndef TBXX_ALIGNMENT_H
#define TBXX_ALIGNMENT_H

#include <cstddef>

namespace tbxx {

  /// Syntactic sugar to test pointer alignment
  /** Example:
   \code
   double *a = ...;
   SCITBX_ASSERT(tbxx::pointer_alignment(p) == 16);
   \endcode
   */
  class pointer_alignment
  {
    void const *p;

  public:
    pointer_alignment(void const *p)
    : p(p)
    {}

    bool operator==(unsigned short n) {
      return reinterpret_cast<std::size_t>(p) % n == 0;
    }


  };

}

#endif
