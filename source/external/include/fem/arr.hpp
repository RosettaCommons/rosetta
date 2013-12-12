#ifndef FEM_ARR_HPP
#define FEM_ARR_HPP

#include <fem/arr_ref.hpp>

namespace fem {

  template <typename T, size_t Ndims=1>
  struct arr : arr_ref<T, Ndims>
  {
    private:
      arr(
        arr const&);

      arr const&
      operator=(
        arr const&);

      public:

    template <size_t BufferNdims>
    explicit
    arr(
      dims<BufferNdims> const& dims)
    :
      arr_ref<T, Ndims>(*(new T[dims.size_1d(Ndims)]), dims)
    {}

    template <size_t BufferNdims>
    arr(
      dims<BufferNdims> const& dims,
      no_fill0_type const&)
    :
      arr_ref<T, Ndims>(*(new T[dims.size_1d(Ndims)]), dims, no_fill0)
    {}

    template <size_t BufferNdims>
    arr(
      dims<BufferNdims> const& dims,
      fill0_type const&)
    :
      arr_ref<T, Ndims>(*(new T[dims.size_1d(Ndims)]), dims, fill0)
    {}

    ~arr()
    {
      delete[] this->elems_;
    }

    operator
    T&() { return *this->begin(); }
  };

} // namespace fem

#endif // GUARD
