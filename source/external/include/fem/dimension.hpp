#ifndef FEM_DIMENSION_HPP
#define FEM_DIMENSION_HPP

#include <fem/arr_and_str_indices.hpp>
#include <fem/error_utils.hpp>
#include <fem/star.hpp>
#include <algorithm>

namespace fem {

  template <size_t Ndims>
  struct dims_base_1
  {
    size_t all[Ndims];

    size_t
    size_1d(
      size_t n_dims=Ndims) const
    {
      size_t result = 1;
      for(size_t i=0;i<n_dims;i++) result *= all[i];
      return result;
    }

    template <size_t BufferNdims>
    void
    set_dims(
      dims_base_1<BufferNdims> const& source)
    {
      for(size_t i=0;i<Ndims;i++) all[i] = source.all[i];
    }
  };

  template <size_t Ndims>
  struct dim_data : dims_base_1<Ndims>
  {
    ssize_t origin[Ndims];

    template <size_t BufferNdims>
    void
    set_dims(
      dim_data<BufferNdims> const& source)
    {
      dims_base_1<Ndims>::set_dims(source);
      for(size_t i=0;i<Ndims;i++) origin[i] = source.origin[i];
    }
  };

  template <size_t Ndims>
  struct dims;

  template <>
  struct dims<1> : dim_data<1>
  {
    size_t
    index_1d(
      ssize_t i1) const
    {
      return i1 - this->origin[0];
    }
  };

  template <>
  struct dims<2> : dim_data<2>
  {
    size_t
    index_1d(
      ssize_t i1,
      ssize_t i2) const
    {
      return
          (i2 - this->origin[1]) * this->all[0]
        + (i1 - this->origin[0]);
    }

  };

  template <>
  struct dims<3> : dim_data<3>
  {
    size_t
    index_1d(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3) const
    {
      return
          ((i3 - this->origin[2])  * this->all[1]
        +  (i2 - this->origin[1])) * this->all[0]
        +  (i1 - this->origin[0]);
    }
  };

  template <>
  struct dims<4> : dim_data<4>
  {
    size_t
    index_1d(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3,
      ssize_t i4) const
    {
      return
          (((i4 - this->origin[3])  * this->all[2]
        +   (i3 - this->origin[2])) * this->all[1]
        +   (i2 - this->origin[1])) * this->all[0]
        +   (i1 - this->origin[0]);
    }
  };

  template <>
  struct dims<5> : dim_data<5>
  {
    size_t
    index_1d(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3,
      ssize_t i4,
      ssize_t i5) const
    {
      return
          ((((i5 - this->origin[4])  * this->all[3]
        +    (i4 - this->origin[3])) * this->all[2]
        +    (i3 - this->origin[2])) * this->all[1]
        +    (i2 - this->origin[1])) * this->all[0]
        +    (i1 - this->origin[0]);
    }
  };

  template <>
  struct dims<6> : dim_data<6>
  {
    size_t
    index_1d(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3,
      ssize_t i4,
      ssize_t i5,
      ssize_t i6) const
    {
      return
          (((((i6 - this->origin[5])  * this->all[4]
        +     (i5 - this->origin[4])) * this->all[3]
        +     (i4 - this->origin[3])) * this->all[2]
        +     (i3 - this->origin[2])) * this->all[1]
        +     (i2 - this->origin[1])) * this->all[0]
        +     (i1 - this->origin[0]);
    }
  };

  struct dim_buffer : dims<arr_dim_max>
  {
    size_t actual_number_of_dimensions;

    dim_buffer() : actual_number_of_dimensions(0) {}

    dim_buffer(
      dim_buffer const& other)
    :
      actual_number_of_dimensions(other.actual_number_of_dimensions)
    {
      std::copy(other.all, other.all+actual_number_of_dimensions, all);
      std::copy(other.origin, other.origin+actual_number_of_dimensions,origin);
    }

    template <size_t BufferNdims>
    dim_buffer(
      dims<BufferNdims> const& source)
    :
      actual_number_of_dimensions(BufferNdims)
    {
      this->set_dims(source);
    }

    template <unsigned I>
    void
    set_origin_all_star()
    {
      actual_number_of_dimensions = (I+1);
      all[I] = size_t_max;
      origin[I] = 1;
    }

    template <unsigned I>
    void
    set_origin_all(
      ssize_t first,
      star_type const&)
    {
      actual_number_of_dimensions = (I+1);
      all[I] = size_t_max - 1;
      origin[I] = first;
    }

    template <unsigned I>
    void
    set_origin_all(
      ssize_t first,
      ssize_t last)
    {
      actual_number_of_dimensions = (I+1);
      all[I] = last - first + 1;
      origin[I] = first;
    }

    size_t
    actual_size_1d() const
    {
      return size_1d(actual_number_of_dimensions);
    }

    template <size_t MaxNdims>
    size_t
    actual_index_1d(
      arr_index_data<MaxNdims> const& arr_index) const
    {
      TBXX_ASSERT(
        arr_index.actual_number_of_dimensions == actual_number_of_dimensions);
      size_t i = actual_number_of_dimensions;
      if (i == 0) return 0;
      i--;
      size_t result = arr_index.elems[i] - origin[i];
      while (i != 0) {
        i--;
        result *= all[i];
        result += arr_index.elems[i] - origin[i];
      }
      return result;
    }
  };

  struct dim_chain : dim_buffer
  {
    dim_chain(
      star_type const&)
    {
      this->set_origin_all_star<0>();
    }

    dim_chain(
      size_t n)
    {
      this->set_origin_all<0>(1, n);
    }

    dim_chain(
      ssize_t first,
      star_type const&)
    {
      this->set_origin_all<0>(first, star_type());
    }

    dim_chain(
      ssize_t first,
      ssize_t last)
    {
      this->set_origin_all<0>(first, last);
    }

    dim_chain&
    dim2(
      star_type const&)
    {
      this->set_origin_all_star<1>();
      return *this;
    }

    dim_chain&
    dim2(
      size_t n)
    {
      this->set_origin_all<1>(1, n);
      return *this;
    }

    dim_chain&
    dim2(
      ssize_t first,
      star_type const&)
    {
      this->set_origin_all<1>(first, star_type());
      return *this;
    }

    dim_chain&
    dim2(
      ssize_t first,
      ssize_t last)
    {
      this->set_origin_all<1>(first, last);
      return *this;
    }

    dim_chain&
    dim3(
      star_type const&)
    {
      this->set_origin_all_star<2>();
      return *this;
    }

    dim_chain&
    dim3(
      size_t n)
    {
      this->set_origin_all<2>(1, n);
      return *this;
    }

    dim_chain&
    dim3(
      ssize_t first,
      star_type const&)
    {
      this->set_origin_all<2>(first, star_type());
      return *this;
    }

    dim_chain&
    dim3(
      ssize_t first,
      ssize_t last)
    {
      this->set_origin_all<2>(first, last);
      return *this;
    }

    dim_chain&
    dim4(
      star_type const&)
    {
      this->set_origin_all_star<3>();
      return *this;
    }

    dim_chain&
    dim4(
      size_t n)
    {
      this->set_origin_all<3>(1, n);
      return *this;
    }

    dim_chain&
    dim4(
      ssize_t first,
      star_type const&)
    {
      this->set_origin_all<3>(first, star_type());
      return *this;
    }

    dim_chain&
    dim4(
      ssize_t first,
      ssize_t last)
    {
      this->set_origin_all<3>(first, last);
      return *this;
    }

    dim_chain&
    dim5(
      star_type const&)
    {
      this->set_origin_all_star<4>();
      return *this;
    }

    dim_chain&
    dim5(
      size_t n)
    {
      this->set_origin_all<4>(1, n);
      return *this;
    }

    dim_chain&
    dim5(
      ssize_t first,
      star_type const&)
    {
      this->set_origin_all<4>(first, star_type());
      return *this;
    }

    dim_chain&
    dim5(
      ssize_t first,
      ssize_t last)
    {
      this->set_origin_all<4>(first, last);
      return *this;
    }

    dim_chain&
    dim6(
      star_type const&)
    {
      this->set_origin_all_star<5>();
      return *this;
    }

    dim_chain&
    dim6(
      size_t n)
    {
      this->set_origin_all<5>(1, n);
      return *this;
    }

    dim_chain&
    dim6(
      ssize_t first,
      star_type const&)
    {
      this->set_origin_all<5>(first, star_type());
      return *this;
    }

    dim_chain&
    dim6(
      ssize_t first,
      ssize_t last)
    {
      this->set_origin_all<5>(first, last);
      return *this;
    }
  };

  inline
  dim_chain
  dim1(
    star_type const&)
  {
    return dim_chain(star_type());
  }

  inline
  dim_chain
  dim1(
    size_t n)
  {
    return dim_chain(n);
  }

  inline
  dim_chain
  dim1(
    ssize_t first,
    star_type const&)
  {
    return dim_chain(first, star_type());
  }

  inline
  dim_chain
  dim1(
    ssize_t first,
    ssize_t last)
  {
    return dim_chain(first, last);
  }

  inline
  dims<1>
  dimension(
    star_type const&)
  {
    dims<1> result;
    result.all[0] = size_t_max;
    result.origin[0] = 1;
    return result;
  }

  inline
  dims<1>
  dimension(
    size_t n)
  {
    dims<1> result;
    result.all[0] = n;
    result.origin[0] = 1;
    return result;
  }

  inline
  dims<2>
  dimension(
    size_t n1,
    star_type const&)
  {
    dims<2> result;
    result.all[0] = n1;
    result.all[1] = size_t_max;
    result.origin[0] = 1;
    result.origin[1] = 1;
    return result;
  }

  inline
  dims<2>
  dimension(
    size_t n1,
    size_t n2)
  {
    dims<2> result;
    result.all[0] = n1;
    result.all[1] = n2;
    result.origin[0] = 1;
    result.origin[1] = 1;
    return result;
  }

  inline
  dims<3>
  dimension(
    size_t n1,
    size_t n2,
    star_type const&)
  {
    dims<3> result;
    result.all[0] = n1;
    result.all[1] = n2;
    result.all[2] = size_t_max;
    result.origin[0] = 1;
    result.origin[1] = 1;
    result.origin[2] = 1;
    return result;
  }

  inline
  dims<3>
  dimension(
    size_t n1,
    size_t n2,
    size_t n3)
  {
    dims<3> result;
    result.all[0] = n1;
    result.all[1] = n2;
    result.all[2] = n3;
    result.origin[0] = 1;
    result.origin[1] = 1;
    result.origin[2] = 1;
    return result;
  }

  inline
  dims<4>
  dimension(
    size_t n1,
    size_t n2,
    size_t n3,
    star_type const&)
  {
    dims<4> result;
    result.all[0] = n1;
    result.all[1] = n2;
    result.all[2] = n3;
    result.all[3] = size_t_max;
    result.origin[0] = 1;
    result.origin[1] = 1;
    result.origin[2] = 1;
    result.origin[3] = 1;
    return result;
  }

  inline
  dims<4>
  dimension(
    size_t n1,
    size_t n2,
    size_t n3,
    size_t n4)
  {
    dims<4> result;
    result.all[0] = n1;
    result.all[1] = n2;
    result.all[2] = n3;
    result.all[3] = n4;
    result.origin[0] = 1;
    result.origin[1] = 1;
    result.origin[2] = 1;
    result.origin[3] = 1;
    return result;
  }

  inline
  dims<5>
  dimension(
    size_t n1,
    size_t n2,
    size_t n3,
    size_t n4,
    star_type const&)
  {
    dims<5> result;
    result.all[0] = n1;
    result.all[1] = n2;
    result.all[2] = n3;
    result.all[3] = n4;
    result.all[4] = size_t_max;
    result.origin[0] = 1;
    result.origin[1] = 1;
    result.origin[2] = 1;
    result.origin[3] = 1;
    result.origin[4] = 1;
    return result;
  }

  inline
  dims<5>
  dimension(
    size_t n1,
    size_t n2,
    size_t n3,
    size_t n4,
    size_t n5)
  {
    dims<5> result;
    result.all[0] = n1;
    result.all[1] = n2;
    result.all[2] = n3;
    result.all[3] = n4;
    result.all[4] = n5;
    result.origin[0] = 1;
    result.origin[1] = 1;
    result.origin[2] = 1;
    result.origin[3] = 1;
    result.origin[4] = 1;
    return result;
  }

  inline
  dims<6>
  dimension(
    size_t n1,
    size_t n2,
    size_t n3,
    size_t n4,
    size_t n5,
    star_type const&)
  {
    dims<6> result;
    result.all[0] = n1;
    result.all[1] = n2;
    result.all[2] = n3;
    result.all[3] = n4;
    result.all[4] = n5;
    result.all[5] = size_t_max;
    result.origin[0] = 1;
    result.origin[1] = 1;
    result.origin[2] = 1;
    result.origin[3] = 1;
    result.origin[4] = 1;
    result.origin[5] = 1;
    return result;
  }

  inline
  dims<6>
  dimension(
    size_t n1,
    size_t n2,
    size_t n3,
    size_t n4,
    size_t n5,
    size_t n6)
  {
    dims<6> result;
    result.all[0] = n1;
    result.all[1] = n2;
    result.all[2] = n3;
    result.all[3] = n4;
    result.all[4] = n5;
    result.all[5] = n6;
    result.origin[0] = 1;
    result.origin[1] = 1;
    result.origin[2] = 1;
    result.origin[3] = 1;
    result.origin[4] = 1;
    result.origin[5] = 1;
    return result;
  }

} // namespace fem

#endif // GUARD
