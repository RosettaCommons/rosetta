#ifndef FEM_DATA_HPP
#define FEM_DATA_HPP

#include <fem/arr.hpp>
#include <fem/str_arr_ref.hpp>
#include <fem/utils/misc.hpp>
#include <typeinfo>
#include <vector>

namespace fem {

  struct datum
  {
    // similar to boost::any
    struct placeholder
    {
      virtual ~placeholder() {}
      virtual std::type_info const& type() const = 0;
      virtual placeholder* clone() const = 0;
    };

    template <typename T>
    struct holder : placeholder
    {
      T held;

      holder(T const& value) : held(value) {}

      virtual std::type_info const&
      type() const { return typeid(T); }

      virtual placeholder*
      clone() const { return new holder(held); }

      private:
        holder const& operator=(holder const&);
    };

    placeholder* content;
    mutable unsigned repeats;

    datum() : content(0), repeats(0) {}

    datum(
      char const& value)
    :
      content(new holder<char>(value)),
      repeats(1)
    {}

    datum(
      bool const& value)
    :
      content(new holder<bool>(value)),
      repeats(1)
    {}

    datum(
      int const& value)
    :
      content(new holder<int>(value)),
      repeats(1)
    {}

    datum(
      float const& value)
    :
      content(new holder<float>(value)),
      repeats(1)
    {}

    datum(
      double const& value)
    :
      content(new holder<double>(value)),
      repeats(1)
    {}

    template <typename T>
    datum(
      T const& value)
    :
      content(new holder<std::string>(std::string(value))),
      repeats(1)
    {}

    datum(
      datum const& other)
    :
       content(other.content ? other.content->clone() : 0),
       repeats(other.repeats)
    {}

    ~datum() { delete content; }

    // required for std::vector<datum> (Apple g++ 4.0.1, Visual C++ 9.0)
    // but not actually used
    datum const& operator=(datum const& rhs)
    {
      throw TBXX_UNREACHABLE_ERROR();
    }

    void
    set_repeats(
      unsigned value) const { repeats = value; }

    void
    throw_type_mismatch(
      char const* target_type) const
    {
      throw std::runtime_error(
        std::string("DATA type mismatch: target type ")
        + target_type
        + " vs. "
        + content->type().name()
        + " value");
    }
  };

  inline
  datum const&
  operator*(
    unsigned repeats,
    datum const& value)
  {
    value.set_repeats(repeats);
    return value;
  }

  struct data_buffer
  {
    utils::slick_ptr<std::vector<datum> > objects;

    data_buffer()
    :
      objects(new std::vector<datum>)
    {}

    template <typename T>
    data_buffer&
    operator,(
      T const& val)
    {
      objects->push_back(datum(val));
      return *this;
    }
  };

  struct data_values
  {
    data_buffer values;
    size_t value_index;
    size_t repeat_index;

    data_values()
    :
      value_index(0),
      repeat_index(0)
    {}

    data_values(
      data_buffer const& values_)
    :
      values(values_),
      value_index(0),
      repeat_index(0)
    {}

    // XXX
    void
    report_types() const
    {
      size_t n = values.objects->size();
      for(size_t i=0;i!=n;i++) {
        datum const& obj = (*(values.objects))[i];
        std::cout <<
          "data_buffer[" << i << "]: " <<
            obj.content->type().name() << std::endl;
      }
    }

    datum const&
    next_datum()
    {
      datum const& result = (*(values.objects))[value_index];
      TBXX_ASSERT(result.content != 0);
      TBXX_ASSERT(result.repeats > 0);
      repeat_index++;
      if (repeat_index == result.repeats) {
        value_index++;
        repeat_index = 0;
      }
      return result;
    }

    data_values&
    operator,(
      bool& val)
    {
      datum const& tab_val = next_datum();
      if (tab_val.content->type() == typeid(bool)) {
        val = static_cast<datum::holder<bool>*>(tab_val.content)->held;
      }
      else {
        tab_val.throw_type_mismatch("bool");
      }
      return *this;
    }

    data_values&
    operator,(
      char& val)
    {
      datum const& tab_val = next_datum();
      if (tab_val.content->type() == typeid(int)) {
        val = static_cast<char>(
          static_cast<datum::holder<int>*>(tab_val.content)->held);
      }
      else if (tab_val.content->type() == typeid(char)) {
        val = static_cast<datum::holder<char>*>(tab_val.content)->held;
      }
      else {
        tab_val.throw_type_mismatch("char");
      }
      return *this;
    }

    data_values&
    operator,(
      integer_star_2& val)
    {
      datum const& tab_val = next_datum();
      if (tab_val.content->type() == typeid(int)) {
        val = static_cast<integer_star_2>(
          static_cast<datum::holder<int>*>(tab_val.content)->held);
      }
      else {
        tab_val.throw_type_mismatch("integer_star_2");
      }
      return *this;
    }

    data_values&
    operator,(
      int& val)
    {
      datum const& tab_val = next_datum();
      if (tab_val.content->type() == typeid(int)) {
        val = static_cast<datum::holder<int>*>(tab_val.content)->held;
      }
      else {
        tab_val.throw_type_mismatch("int");
      }
      return *this;
    }

    data_values&
    operator,(
      float& val)
    {
      datum const& tab_val = next_datum();
      if (tab_val.content->type() == typeid(int)) {
        val = static_cast<float>(
          static_cast<datum::holder<int>*>(tab_val.content)->held);
      }
      else if (tab_val.content->type() == typeid(float)) {
        val = static_cast<datum::holder<float>*>(tab_val.content)->held;
      }
      else {
        tab_val.throw_type_mismatch("float");
      }
      return *this;
    }

    data_values&
    operator,(
      double& val)
    {
      datum const& tab_val = next_datum();
      if (tab_val.content->type() == typeid(int)) {
        val = static_cast<double>(
          static_cast<datum::holder<int>*>(tab_val.content)->held);
      }
      else if (tab_val.content->type() == typeid(float)) {
        val = static_cast<double>(
          static_cast<datum::holder<float>*>(tab_val.content)->held);
      }
      else if (tab_val.content->type() == typeid(double)) {
        val = static_cast<datum::holder<double>*>(tab_val.content)->held;
      }
      else {
        tab_val.throw_type_mismatch("double");
      }
      return *this;
    }

    data_values&
    operator,(
      str_ref val)
    {
      datum const& tab_val = next_datum();
      if (tab_val.content->type() == typeid(std::string)) {
        std::string const& tab_str = static_cast<datum::holder<std::string>*>(
          tab_val.content)->held;
        val = tab_str.c_str();
      }
      else {
        tab_val.throw_type_mismatch("string");
      }
      return *this;
    }

    template <typename T, size_t Ndims>
    data_values&
    operator,(
      arr_ref<T, Ndims>& val)
    {
      size_t n = val.size_1d();
      T* val_begin = val.begin();
      for(size_t i=0;i<n;i++) (*this), val_begin[i];
      return *this;
    }

    template <size_t Ndims>
    data_values&
    operator,(
      str_arr_ref<Ndims>& val)
    {
      size_t n = val.size_1d();
      for(size_t i=0;i<n;i++) (*this), val[i];
      return *this;
    }
  };

  struct values_type
  {
    template <typename T>
    data_buffer
    operator,(
      T const& val) const
    {
      return data_buffer(), val;
    }
  };

  static const values_type values = values_type();

  inline
  data_values
  data(
    data_buffer const& values)
  {
    return data_values(values);
  }
}

#endif // GUARD
