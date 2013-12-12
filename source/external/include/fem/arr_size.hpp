#ifndef FEM_ARR_SIZE_HPP
#define FEM_ARR_SIZE_HPP

#include <fem/arr_ref.hpp>

namespace fem {

  template <size_t Size, typename T, size_t Ndims=1>
  struct arr_size : arr_ref<T, Ndims>
  {
    T elems_memory[Size];

    private:
      arr_size(
        arr_size const&);

      arr_size const&
      operator=(
        arr_size const&);

      public:

    arr_size()
    :
      arr_ref<T, Ndims>(*elems_memory)
    {}

    template <size_t BufferNdims>
    explicit
    arr_size(
      dims<BufferNdims> const& dims)
    :
      arr_ref<T, Ndims>(*elems_memory, dims)
    {}

    template <size_t BufferNdims>
    arr_size(
      dims<BufferNdims> const& dims,
      no_fill0_type const&)
    :
      arr_ref<T, Ndims>(*elems_memory, dims, no_fill0)
    {}

    template <size_t BufferNdims>
    arr_size(
      dims<BufferNdims> const& dims,
      fill0_type const&)
    :
      arr_ref<T, Ndims>(*elems_memory, dims, fill0)
    {}

    operator
    T&() { return *elems_memory; }
  };

  template <size_t Size, typename T>
  struct arr_1d : arr_size<Size, T, 1>
  {
    typedef arr_size<Size, T, 1> base_t;

    arr_1d()
    :
      base_t(dimension(Size))
    {}

    explicit
    arr_1d(
      no_fill0_type const&)
    :
      base_t(dimension(Size), no_fill0)
    {}

    explicit
    arr_1d(
      fill0_type const&)
    :
      base_t(dimension(Size), fill0)
    {}

    template <size_t BufferNdims>
    explicit
    arr_1d(
      dims<BufferNdims> const& dims)
    :
      base_t(dims)
    {}

    template <size_t BufferNdims>
    arr_1d(
      dims<BufferNdims> const& dims,
      no_fill0_type const&)
    :
      base_t(dims, no_fill0)
    {}

    template <size_t BufferNdims>
    arr_1d(
      dims<BufferNdims> const& dims,
      fill0_type const&)
    :
      base_t(dims, fill0)
    {}
  };

  template <size_t Size1, size_t Size2, typename T>
  struct arr_2d : arr_size<Size1*Size2, T, 2>
  {
    typedef arr_size<Size1*Size2, T, 2> base_t;

    arr_2d()
    :
      base_t(dimension(Size1, Size2))
    {}

    explicit
    arr_2d(
      no_fill0_type const&)
    :
      base_t(dimension(Size1, Size2), no_fill0)
    {}

    explicit
    arr_2d(
      fill0_type const&)
    :
      base_t(dimension(Size1, Size2), fill0)
    {}

    template <size_t BufferNdims>
    explicit
    arr_2d(
      dims<BufferNdims> const& dims)
    :
      base_t(dims)
    {}

    template <size_t BufferNdims>
    arr_2d(
      dims<BufferNdims> const& dims,
      no_fill0_type const&)
    :
      base_t(dims, no_fill0)
    {}

    template <size_t BufferNdims>
    arr_2d(
      dims<BufferNdims> const& dims,
      fill0_type const&)
    :
      base_t(dims, fill0)
    {}
  };

  template <size_t Size1, size_t Size2, size_t Size3, typename T>
  struct arr_3d : arr_size<Size1*Size2*Size3, T, 3>
  {
    typedef arr_size<Size1*Size2*Size3, T, 3> base_t;

    arr_3d()
    :
      base_t(dimension(Size1, Size2, Size3))
    {}

    explicit
    arr_3d(
      no_fill0_type const&)
    :
      base_t(dimension(Size1, Size2, Size3), no_fill0)
    {}

    explicit
    arr_3d(
      fill0_type const&)
    :
      base_t(dimension(Size1, Size2, Size3), fill0)
    {}

    template <size_t BufferNdims>
    explicit
    arr_3d(
      dims<BufferNdims> const& dims)
    :
      base_t(dims)
    {}

    template <size_t BufferNdims>
    arr_3d(
      dims<BufferNdims> const& dims,
      no_fill0_type const&)
    :
      base_t(dims, no_fill0)
    {}

    template <size_t BufferNdims>
    arr_3d(
      dims<BufferNdims> const& dims,
      fill0_type const&)
    :
      base_t(dims, fill0)
    {}
  };

} // namespace fem

#endif // GUARD
