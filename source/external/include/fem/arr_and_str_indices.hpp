#ifndef FEM_ARR_AND_STR_INDICES_HPP
#define FEM_ARR_AND_STR_INDICES_HPP

#include <fem/size_t.hpp>

namespace fem {

  static const size_t arr_dim_max = 6;

  template <size_t Ndims>
  struct arr_and_str_indices;

  template <size_t Ndims>
  struct arr_index_data
  {
    size_t actual_number_of_dimensions;
    ssize_t elems[Ndims];

    arr_index_data() : actual_number_of_dimensions(0) {}

    explicit
    arr_index_data(
      ssize_t i1) : actual_number_of_dimensions(1)
    {
      elems[0] = i1;
    }

    arr_index_data(
      ssize_t i1,
      ssize_t i2) : actual_number_of_dimensions(2)
    {
      elems[0] = i1;
      elems[1] = i2;
    }

    arr_index_data(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3) : actual_number_of_dimensions(3)
    {
      elems[0] = i1;
      elems[1] = i2;
      elems[2] = i3;
    }

    arr_index_data(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3,
      ssize_t i4) : actual_number_of_dimensions(4)
    {
      elems[0] = i1;
      elems[1] = i2;
      elems[2] = i3;
      elems[3] = i4;
    }

    arr_index_data(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3,
      ssize_t i4,
      ssize_t i5) : actual_number_of_dimensions(5)
    {
      elems[0] = i1;
      elems[1] = i2;
      elems[2] = i3;
      elems[3] = i4;
      elems[4] = i5;
    }

    arr_index_data(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3,
      ssize_t i4,
      ssize_t i5,
      ssize_t i6) : actual_number_of_dimensions(6)
    {
      elems[0] = i1;
      elems[1] = i2;
      elems[2] = i3;
      elems[3] = i4;
      elems[4] = i5;
      elems[5] = i6;
    }

    inline
    arr_and_str_indices<Ndims>
    operator()(
      int first,
      int last) const;
  };

  typedef arr_index_data<arr_dim_max> arr_index;

  struct str_index
  {
    int first;
    int last;

    str_index(int first_, int last_) : first(first_), last(last_) {}
  };

  template <size_t Ndims>
  struct arr_and_str_indices
  {
    arr_index_data<Ndims> arr_ix;
    str_index str_ix;

    arr_and_str_indices(
      arr_index_data<Ndims> const& arr_ix_,
      str_index const& str_ix_)
    :
      arr_ix(arr_ix_),
      str_ix(str_ix_)
    {}
  };

  template <size_t Ndims>
  inline
  arr_and_str_indices<Ndims>
  arr_index_data<Ndims>::operator()(
    int first,
    int last) const
  {
    return arr_and_str_indices<Ndims>(*this, str_index(first, last));
  }

} // namespace fem

#endif // GUARD
