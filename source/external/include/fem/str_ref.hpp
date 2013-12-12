#ifndef FEM_STR_REF_HPP
#define FEM_STR_REF_HPP

#include <fem/utils/misc.hpp>
#include <fem/utils/string.hpp>

namespace fem {

  struct str_cref
  {
    protected:
      char const* elems_;
      int len_;

      public:

    str_cref(
      char const* elems,
      int len)
    :
      elems_(elems),
      len_(len)
    {}

    str_cref(
      char const* elems)
    :
      elems_(elems),
      len_(std::strlen(elems))
    {}

    str_cref(
      str_cref const& other,
      int len)
    :
      elems_(other.elems()),
      len_(len)
    {}

    char const*
    elems() const { return elems_; }

    int
    len() const { return len_; }

    operator
    std::string() const
    {
      return std::string(elems_, len_);
    }

    char const&
    operator[](
      int i) const
    {
      return elems_[i];
    }

    str_cref
    operator()(
      int first,
      int last) const
    {
      return str_cref(elems_+first-1, last-first+1);
    }

    bool
    operator==(
      char const* rhs) const
    {
      return utils::string_eq(elems_, len_, rhs);
    }

    bool
    operator!=(
      char const* rhs) const
    {
      return !((*this) == rhs);
    }

    bool
    operator==(
      str_cref const& rhs) const
    {
      return utils::string_eq(elems_, len_, rhs.elems_, rhs.len_);
    }

    bool
    operator!=(
      str_cref const& rhs) const
    {
      return !((*this) == rhs);
    }

    bool
    operator<(
      str_cref const& rhs) const
    {
      return (utils::string_compare_lexical(
        elems(), len(), rhs.elems(), rhs.len()) < 0);
    }

    bool
    operator<=(
      str_cref const& rhs) const
    {
      return (utils::string_compare_lexical(
        elems(), len(), rhs.elems(), rhs.len()) <= 0);
    }

    bool
    operator>(
      str_cref const& rhs) const
    {
      return (utils::string_compare_lexical(
        elems(), len(), rhs.elems(), rhs.len()) > 0);
    }

    bool
    operator>=(
      str_cref const& rhs) const
    {
      return (utils::string_compare_lexical(
        elems(), len(), rhs.elems(), rhs.len()) >= 0);
    }
  };

  int
  inline
  len(
    str_cref const& s) { return s.len(); }

  bool
  inline
  operator==(
    char const* lhs,
    str_cref const& rhs) { return (rhs == lhs); }

  bool
  inline
  operator!=(
    char const* lhs,
    str_cref const& rhs) { return (rhs != lhs); }

  struct str_addends
  {
    str_cref lhs, rhs;
    std::string sum;

    str_addends(
      str_cref const& lhs_,
      str_cref const& rhs_)
    :
      lhs(lhs_),
      rhs(rhs_)
    {}

    operator
    str_cref()
    {
      int sum_len = lhs.len() + rhs.len();
      if (sum_len != 0 && sum.size() == 0) {
        sum = std::string(lhs) + std::string(rhs);
      }
      return str_cref(sum.data(), sum_len);
    }
  };

  inline
  str_addends
  operator+(
    str_cref const& lhs,
    str_cref const& rhs)
  {
    return str_addends(lhs, rhs);
  }

  template <int StrLen>
  struct str;

  struct str_ref : str_cref
  {
    str_ref(
      char* elems,
      int len)
    :
      str_cref(elems, len)
    {}

    str_ref(
      str_ref const& other,
      int len)
    :
      str_cref(other, len)
    {}

    char*
    elems() const { return const_cast<char*>(elems_); }

    char&
    operator[](
      int i) const
    {
      return elems()[i];
    }

    void
    operator=(
      char const* rhs)
    {
      utils::copy_with_blank_padding(rhs, elems(), len());
    }

    void
    operator=(
      std::string const& rhs)
    {
      utils::copy_with_blank_padding(rhs.data(), rhs.size(), elems(), len());
    }

    void
    operator=(
      str_cref const& rhs)
    {
      utils::copy_with_blank_padding(rhs.elems(), rhs.len(), elems(), len());
    }

    void
    operator=(
      str_ref const& rhs) { (*this) = static_cast<str_cref>(rhs); }

    void
    operator=(
      str_addends const& addends)
    {
      int ll = addends.lhs.len();
      int n_from_rhs = len() - ll;
      if (n_from_rhs <= 0) {
        std::memmove(elems(), addends.lhs.elems(), len());
      }
      else {
        utils::simple_buffer<char> buffer((len()));
        std::memcpy(buffer.space, addends.lhs.elems(), ll);
        utils::copy_with_blank_padding(
          addends.rhs.elems(), addends.rhs.len(), buffer.space+ll, n_from_rhs);
        std::memcpy(elems(), buffer.space, len());
      }
    }

    template <int StrLen>
    void
    operator=(
      str<StrLen> const& rhs);

    str_ref
    operator()(
      int first,
      int last) const
    {
      return str_ref(elems()+first-1, last-first+1);
    }
  };

} // namespace fem

#endif // GUARD
