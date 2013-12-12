#ifndef FEM_STR_ARR_REF_HPP
#define FEM_STR_ARR_REF_HPP

#include <fem/arr_ref.hpp>
#include <fem/str.hpp>

namespace fem {

  template <size_t Ndims=1>
  struct str_arr_cref
  {
    protected:
      char const* elems_;
      int len_;
      dims<Ndims> dims_;

      str_arr_cref() {}

      public:

    str_arr_cref(
      char const* elems,
      int len)
    :
      elems_(elems),
      len_(len)
    {}

    template <size_t BufferNdims>
    str_arr_cref(
      char const* elems,
      int len,
      dims<BufferNdims> const& dims)
    :
      elems_(elems),
      len_(len)
    {
      (*this)(dims);
    }

    str_arr_cref(
      str_cref const& other)
    :
      elems_(other.elems()),
      len_(other.len())
    {}

    template <size_t BufferNdims>
    str_arr_cref(
      str_cref const& other,
      dims<BufferNdims> const& dims)
    :
      elems_(other.elems()),
      len_(other.len())
    {
      (*this)(dims);
    }

    str_arr_cref(
      str_cref const& other,
      int len)
    :
      elems_(other.elems()),
      len_(len)
    {}

    template <size_t BufferNdims>
    str_arr_cref(
      str_cref const& other,
      int len,
      dims<BufferNdims> const& dims)
    :
      elems_(other.elems()),
      len_(len)
    {
      (*this)(dims);
    }

    template <int StrLen>
    str_arr_cref(
      str<StrLen> const& other)
    :
      elems_(other.elems),
      len_(StrLen)
    {}

    template <int StrLen, size_t BufferNdims>
    str_arr_cref(
      str<StrLen> const& other,
      dims<BufferNdims> const& dims)
    :
      elems_(other.elems),
      len_(StrLen)
    {
      (*this)(dims);
    }

    template <int StrLen, size_t OtherNdims>
    str_arr_cref(
      arr_cref<str<StrLen>, OtherNdims> const& other)
    :
      elems_(other.begin()->elems),
      len_(StrLen)
    {}

    template <int StrLen, size_t OtherNdims, size_t BufferNdims>
    str_arr_cref(
      arr_cref<str<StrLen>, OtherNdims> const& other,
      dims<BufferNdims> const& dims)
    :
      elems_(other.begin()->elems),
      len_(StrLen)
    {
      (*this)(dims);
    }

    template <size_t OtherNdims, size_t BufferNdims>
    str_arr_cref(
      str_arr_cref<OtherNdims> const& other,
      dims<BufferNdims> const& dims)
    :
      elems_(other.begin()),
      len_(other.len())
    {
      (*this)(dims);
    }

    char const*
    begin() const { return elems_; }

    int
    len() const { return len_; }

    template <size_t BufferNdims>
    void
    operator()(
      dims<BufferNdims> const& dims)
    {
      dims_.set_dims(dims);
    }

    size_t
    size_1d() const { return dims_.size_1d(); }

    str_cref
    operator[](
      size_t i) const
    {
      return str_cref(&elems_[len_ * i], len_);
    }

    str_cref
    operator()(
      ssize_t i1) const
    {
      return str_cref(&elems_[len_ * dims_.index_1d(i1)], len_);
    }

    str_cref
    operator()(
      ssize_t i1,
      ssize_t i2) const
    {
      return str_cref(&elems_[len_ * dims_.index_1d(i1, i2)], len_);
    }

    str_cref
    operator()(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3) const
    {
      return str_cref(&elems_[len_ * dims_.index_1d(i1, i2, i3)], len_);
    }

    str_cref
    operator()(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3,
      ssize_t i4) const
    {
      return str_cref(&elems_[len_ * dims_.index_1d(i1, i2, i3, i4)], len_);
    }

    str_cref
    operator()(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3,
      ssize_t i4,
      ssize_t i5) const
    {
      return str_cref(
        &elems_[len_ * dims_.index_1d(i1, i2, i3, i4, i5)], len_);
    }

    str_cref
    operator()(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3,
      ssize_t i4,
      ssize_t i5,
      ssize_t i6) const
    {
      return str_cref(
        &elems_[len_ * dims_.index_1d(i1, i2, i3, i4, i5, i6)], len_);
    }
  };

  template <size_t Ndims=1>
  struct str_arr_ref : str_arr_cref<Ndims>
  {
    protected:
      str_arr_ref() {}

      public:

    str_arr_ref(
      char* elems,
      int len)
    :
      str_arr_cref<Ndims>(elems, len)
    {}

    template <size_t BufferNdims>
    str_arr_ref(
      char const* elems,
      int len,
      dims<BufferNdims> const& dims)
    :
      str_arr_cref<Ndims>(elems, len)
    {
      (*this)(dims);
    }

    str_arr_ref(
      str_ref const& other)
    :
      str_arr_cref<Ndims>(other)
    {}

    template <size_t BufferNdims>
    str_arr_ref(
      str_ref const& other,
      dims<BufferNdims> const& dims)
    :
      str_arr_cref<Ndims>(other, dims)
    {}

    str_arr_ref(
      str_ref const& other,
      int len)
    :
      str_arr_cref<Ndims>(other, len)
    {}

    template <size_t BufferNdims>
    str_arr_ref(
      str_ref const& other,
      int len,
      dims<BufferNdims> const& dims)
    :
      str_arr_cref<Ndims>(other, len, dims)
    {}

    template <int StrLen>
    str_arr_ref(
      str<StrLen> const& other)
    :
      str_arr_cref<Ndims>(other)
    {}

    template <int StrLen, size_t OtherNdims>
    str_arr_ref(
      arr_ref<str<StrLen>, OtherNdims> const& other)
    :
      str_arr_cref<Ndims>(other)
    {}

    template <int StrLen, size_t BufferNdims>
    str_arr_ref(
      str<StrLen> const& other,
      int len,
      dims<BufferNdims> const& dims)
    :
      str_arr_cref<Ndims>(other.elems, len, dims)
    {}

    template <int StrLen, size_t BufferNdims>
    str_arr_ref(
      str<StrLen> const& other,
      dims<BufferNdims> const& dims)
    :
      str_arr_cref<Ndims>(other, dims)
    {}

    template <int StrLen, size_t BufferNdims>
    str_arr_ref(
      str<StrLen> const& other,
      dims<BufferNdims> const& dims,
      no_fill0_type const&)
    :
      str_arr_cref<Ndims>(other, dims)
    {}

    template <int StrLen, size_t BufferNdims>
    str_arr_ref(
      str<StrLen> const& other,
      dims<BufferNdims> const& dims,
      fill0_type const&)
    :
      str_arr_cref<Ndims>(other, dims)
    {
      std::memset(this->begin(), char0, StrLen * this->size_1d());
    }

    template <int StrLen, size_t OtherNdims, size_t BufferNdims>
    str_arr_ref(
      arr_ref<str<StrLen>, OtherNdims> const& other,
      dims<BufferNdims> const& dims)
    :
      str_arr_cref<Ndims>(other, dims)
    {}

    template <size_t OtherNdims, size_t BufferNdims>
    str_arr_ref(
      str_arr_ref<OtherNdims> const& other,
      dims<BufferNdims> const& dims)
    :
      str_arr_cref<Ndims>(other, dims)
    {}

    char*
    begin() const { return const_cast<char*>(this->elems_); }

    operator
    char*() const { return this->begin(); }

    template <size_t BufferNdims>
    void
    operator()(
      dims<BufferNdims> const& dims)
    {
      str_arr_cref<Ndims>::operator()(dims);
    }

    str_ref
    operator[](
      size_t i) const
    {
      return str_ref(
        &this->begin()[this->len() * i],
        this->len());
    }

    str_ref
    operator()(
      ssize_t i1) const
    {
      return str_ref(
        &this->begin()[this->len() * this->dims_.index_1d(i1)],
        this->len());
    }

    str_ref
    operator()(
      ssize_t i1,
      ssize_t i2) const
    {
      return str_ref(
        &this->begin()[this->len() * this->dims_.index_1d(i1, i2)],
        this->len());
    }

    str_ref
    operator()(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3) const
    {
      return str_ref(
        &this->begin()[this->len() * this->dims_.index_1d(i1, i2, i3)],
        this->len());
    }

    str_ref
    operator()(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3,
      ssize_t i4) const
    {
      return str_ref(
        &this->begin()[this->len() * this->dims_.index_1d(i1, i2, i3, i4)],
        this->len());
    }

    str_ref
    operator()(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3,
      ssize_t i4,
      ssize_t i5) const
    {
      return str_ref(
        &this->begin()[this->len() * this->dims_.index_1d(i1, i2, i3, i4, i5)],
        this->len());
    }

    str_ref
    operator()(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3,
      ssize_t i4,
      ssize_t i5,
      ssize_t i6) const
    {
      return str_ref(
        &this->begin()[this->len() * this->dims_.index_1d(i1,i2,i3,i4,i5,i6)],
        this->len());
    }
  };

} // namespace fem

#endif // GUARD
