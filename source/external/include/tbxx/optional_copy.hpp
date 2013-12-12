#ifndef TBXX_OPTIONAL_COPY_H
#define TBXX_OPTIONAL_COPY_H

namespace tbxx {

  //! Optional allocation of value with new, with value-copy semantics.
  template <typename ValueType>
  class optional_copy
  {
    public:
      typedef ValueType value_type;

    protected:
      ValueType* ptr_;

    public:
      optional_copy() : ptr_(0) {}

      ~optional_copy()
      {
        delete ptr_;
      }

      optional_copy(
        optional_copy const& other)
      :
        ptr_(other.ptr_ == 0 ? 0 : new ValueType(*(other.ptr_)))
      {}

      optional_copy&
      operator=(
        optional_copy const& other)
      {
        delete ptr_;
        ptr_ = (other.ptr_ == 0 ? 0 : new ValueType(*(other.ptr_)));
        return *this;
      }

      explicit
      optional_copy(
        ValueType const& value)
      :
        ptr_(new ValueType(value))
      {}

      optional_copy&
      operator=(
        ValueType const& value)
      {
        delete ptr_;
        ptr_ = new ValueType(value);
        return *this;
      }

      void
      release()
      {
        delete ptr_;
        ptr_ = 0;
      }

      operator bool() const { return (ptr_ != 0); }

      ValueType*
      get()       { return ptr_; }

      ValueType const*
      get() const { return ptr_; }

      ValueType*
      operator->()       { return ptr_; }

      ValueType const*
      operator->() const { return ptr_; }

      ValueType&
      operator*()       { return *ptr_; }

      ValueType const&
      operator*() const { return *ptr_; }
  };

  /*! \brief Same behavior as optional_copy, but with operator[] as
      shortcut for (*get())[]
   */
  template <typename ContainerType>
  class optional_container : public optional_copy<ContainerType>
  {
    public:
      optional_container() {}

      optional_container(
        optional_container const& other)
      :
        optional_copy<ContainerType>(other)
      {}

      optional_container&
      operator=(
        optional_container const& other)
      {
        delete this->ptr_;
        this->ptr_ = (other.ptr_ == 0 ? 0 : new ContainerType(*(other.ptr_)));
        return *this;
      }

      explicit
      optional_container(
        ContainerType const& container)
      :
        optional_copy<ContainerType>(container)
      {}

      optional_container&
      operator=(
        ContainerType const& container)
      {
        delete this->ptr_;
        this->ptr_ = new ContainerType(container);
        return *this;
      }

      typename ContainerType::value_type&
      operator[](typename ContainerType::size_type const& i)
      {
        return (*this->ptr_)[i];
      }

      typename ContainerType::value_type const&
      operator[](typename ContainerType::size_type const& i) const
      {
        return (*this->ptr_)[i];
      }
  };

} // namespace tbxx

#endif // GUARD
