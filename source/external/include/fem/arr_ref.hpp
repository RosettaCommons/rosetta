#ifndef FEM_ARR_REF_HPP
#define FEM_ARR_REF_HPP

#include <fem/dimension.hpp>
#include <fem/zero.hpp>
#include <cstring>

namespace fem {

  template <typename T, size_t Ndims=1>
  struct arr_cref
  {
    protected:
      T const* elems_;
      dims<Ndims> dims_;

      arr_cref() {}

		public:

    arr_cref(
      T const& val)
    :
      elems_(&val)
    {}

    template <size_t OtherNdims>
    arr_cref(
      arr_cref<T, OtherNdims> const& other)
    :
      elems_(other.begin())
    {}

    template <size_t OtherNdims, size_t BufferNdims>
    arr_cref(
      arr_cref<T, OtherNdims> const& other,
      dims<BufferNdims> const& dims)
    :
      elems_(other.begin())
    {
      (*this)(dims);
    }

    template <size_t BufferNdims>
    arr_cref(
      T const& val,
      dims<BufferNdims> const& dims)
    :
      elems_(&val)
    {
      (*this)(dims);
    }

    T const*
    begin() const { return elems_; }

    T const*
    end() const { return elems_ + size_1d(); }

#if defined(_MSC_VER)
    // workaround, breaks const-correctness, but only on this platform
    operator
    T&() const { return const_cast<T*>(elems_)[0]; }
#else
    operator
    T const&() const { return elems_[0]; }
#endif

    template <size_t BufferNdims>
    void
    operator()(
      dims<BufferNdims> const& dims)
    {
      dims_.set_dims(dims);
    }

    size_t
    size_1d() const { return dims_.size_1d(); }

    T const&
    operator()(
      ssize_t i1) const
    {
      return elems_[dims_.index_1d(i1)];
    }

    T const&
    operator()(
      ssize_t i1,
      ssize_t i2) const
    {
      return elems_[dims_.index_1d(i1, i2)];
    }

    T const&
    operator()(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3) const
    {
      return elems_[dims_.index_1d(i1, i2, i3)];
    }

    T const&
    operator()(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3,
      ssize_t i4) const
    {
      return elems_[dims_.index_1d(i1, i2, i3, i4)];
    }

    T const&
    operator()(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3,
      ssize_t i4,
      ssize_t i5) const
    {
      return elems_[dims_.index_1d(i1, i2, i3, i4, i5)];
    }

    T const&
    operator()(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3,
      ssize_t i4,
      ssize_t i5,
      ssize_t i6) const
    {
      return elems_[dims_.index_1d(i1, i2, i3, i4, i5, i6)];
    }
  };

  template <typename T, size_t Ndims=1>
  struct arr_ref : arr_cref<T, Ndims>
  {
    protected:
      arr_ref() {}

		public:

    arr_ref(
      T& val)
    :
      arr_cref<T, Ndims>(val)
    {}

    template <size_t OtherNdims>
    arr_ref(
      arr_ref<T, OtherNdims> const& other)
    :
      arr_cref<T, Ndims>(other)
    {}

    template <size_t BufferNdims>
    arr_ref(
      T& val,
      dims<BufferNdims> const& dims)
    :
      arr_cref<T, Ndims>(val, dims)
    {}

    template <size_t BufferNdims>
    arr_ref(
      T& val,
      dims<BufferNdims> const& dims,
      no_fill0_type const&)
    :
      arr_cref<T, Ndims>(val, dims)
    {}

    template <size_t BufferNdims>
    arr_ref(
      T& val,
      dims<BufferNdims> const& dims,
      fill0_type const&)
    :
      arr_cref<T, Ndims>(val, dims)
    {
      std::memset(this->begin(), 0, this->dims_.size_1d() * sizeof(T));
    }

    T*
    begin() const { return const_cast<T*>(this->elems_); }

#if !defined(_MSC_VER)
    operator
    T&() const { return *(this->begin()); }
#endif

    template <size_t BufferNdims>
    void
    operator()(
      dims<BufferNdims> const& dims)
    {
      arr_cref<T, Ndims>::operator()(dims);
    }

    T&
    operator()(
      ssize_t i1) const
    {
      return this->begin()[this->dims_.index_1d(i1)];
    }

    T&
    operator()(
      ssize_t i1,
      ssize_t i2) const
    {
      return this->begin()[this->dims_.index_1d(i1, i2)];
    }

    T&
    operator()(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3) const
    {
      return this->begin()[this->dims_.index_1d(i1, i2, i3)];
    }

    T&
    operator()(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3,
      ssize_t i4) const
    {
      return this->begin()[this->dims_.index_1d(i1, i2, i3, i4)];
    }

    T&
    operator()(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3,
      ssize_t i4,
      ssize_t i5) const
    {
      return this->begin()[this->dims_.index_1d(i1, i2, i3, i4, i5)];
    }

    T&
    operator()(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3,
      ssize_t i4,
      ssize_t i5,
      ssize_t i6) const
    {
      return this->begin()[this->dims_.index_1d(i1, i2, i3, i4, i5, i6)];
    }
  };

} // namespace fem

#endif // GUARD
