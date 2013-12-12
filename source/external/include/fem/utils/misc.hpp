#ifndef FEM_UTILS_MISC_HPP
#define FEM_UTILS_MISC_HPP

#include <fem/size_t.hpp>

namespace fem { namespace utils {

  template <typename T>
  struct hide
  {
    T hidden;
  };

  struct noncopyable // clone of boost::noncopyable
  {
    protected:
      noncopyable() {}

    private:
      noncopyable(noncopyable const&);
      noncopyable const& operator=(noncopyable const&);
  };

  template <typename T, size_t SmallSize=256>
  struct simple_buffer : noncopyable
  {
    T small_space[SmallSize];
    T* space;

    explicit
    simple_buffer(
      size_t size)
    :
      space(size <= SmallSize ? small_space : new T[size])
    {}

    ~simple_buffer()
    {
      if (space != small_space) delete[] space;
    }
  };

  // simlar to std::auto_ptr, but all member functions const
  template<typename T>
  struct slick_ptr
  {
    protected:
      mutable T* ptr;

      public:

    explicit
    slick_ptr(T* p=0)
    :
      ptr(p)
    {}

    slick_ptr(
      slick_ptr const& other)
    :
      ptr(other.release())
    {}

    slick_ptr&
    operator=(
      slick_ptr const& other)
    {
      reset(other.release());
      return *this;
    }

    ~slick_ptr() { delete ptr; }

    void
    reset(T* p=0) const
    {
      if (p != ptr) slick_ptr<T>(p).swap(*this);
    }

    T*
    release() const
    {
      T* result = ptr;
      ptr = 0;
      return result;
    }

    T*
    get() const { return ptr; }

    T*
    operator->() const { return ptr; }

    T&
    operator*() const { return *ptr; }

    operator bool() const { return (ptr != 0); }

    void
    swap(
      slick_ptr const& other) const
    {
      T* tmp = other.ptr;
      other.ptr = ptr;
      ptr = tmp;
    }
  };

}} // namespace fem::utils

#endif // GUARD
