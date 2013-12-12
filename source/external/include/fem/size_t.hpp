#ifndef FEM_SIZE_T_HPP
#define FEM_SIZE_T_HPP

#include <limits>
#include <cstddef>

namespace fem {

  typedef std::size_t size_t;
  typedef std::ptrdiff_t ssize_t;

  static const size_t size_t_max = std::numeric_limits<size_t>::max();
  static const ssize_t ssize_t_max = std::numeric_limits<ssize_t>::max();

  template <typename T>
  struct array_of_2
  {
    T elems[2];

    array_of_2() {}

    array_of_2(T const& i, T const& j) { elems[0] = i; elems[1] = j; }
  };

  typedef array_of_2<size_t> size_t_2;
  typedef array_of_2<ssize_t> ssize_t_2;

} // namespace fem

#endif // GUARD
