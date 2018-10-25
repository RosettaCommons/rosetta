#pragma once

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

#include "Python.h"
#include "ndarray.h"
#include "ndarray/converter/pointer.h"

NAMESPACE_BEGIN(pybind11)
NAMESPACE_BEGIN(detail)

/* @brief A pybind11 type_caster for ndarray::Array
 */
template <typename Element, int N, int C>
class type_caster< ndarray::Array<Element,N,C> > {
public:
    typedef ndarray::Array<Element,N,C> array_type;
    static constexpr auto array_name = _("ndarray::Array[") +
      make_caster<Element>::name + _(", ") +
      _("ndim=") + _<N>() + _(", ") +
      _("contiguous=") + _<C>() + _("]");

    PYBIND11_TYPE_CASTER(array_type, array_name);

    bool load(handle src, bool) {

        auto buf = array_t<Element>::ensure(src);
        if (!buf) {
            return false;
        }

        auto dims = buf.ndim();
        if (dims != N) {
            return false;
        }

        for(std::size_t s = 0; s < dims; ++s) {
          if(buf.strides(s) % sizeof(Element) != 0) {
            return false;
          }
        }

        if (!boost::is_const<Element>::value && !buf.writeable()){
          return false;
        }

        if (C > 0) {
          ndarray::Offset requiredStride = sizeof(Element);
            for (int i = 0; i < C; ++i) {
              ndarray::Offset actualStride = buf.strides(N-i-1);
                if (actualStride != requiredStride) {
                    return false;
                }
                requiredStride = buf.shape(N-i-1) * buf.strides(N-i-1);
            }
        } else if (C < 0) {
          ndarray::Offset requiredStride = sizeof(Element);
            for (int i = 0; i < -C; ++i) {
              ndarray::Offset actualStride = buf.strides(i);
                if (actualStride != requiredStride) {
                    return false;
                }
                requiredStride = buf.shape(i) * buf.strides(i);
            }
        }

        ndarray::Vector<ndarray::Size, N> shape;
        ndarray::Vector<ndarray::Size, N> strides;
        for(int i = 0; i < N; ++i) {
          shape[i] = buf.shape(i);
          strides[i] = buf.strides(i) / sizeof(Element);
        }

        ndarray::PyPtr buf_ptr;
        buf_ptr.reset(buf.ptr());

        value = external(
          const_cast<Element*>(buf.data()),
          shape, strides, buf_ptr
        );

        return true;
    }
    static handle cast(const ndarray::Array<Element,N,C> &src, return_value_policy /* policy */, handle /* parent */) {
        std::vector<size_t> shape;
        std::vector<size_t> strides;
        for(int i = 0; i < N; ++i) {
          shape.push_back(src.getShape()[i]);
          strides.push_back(src.getStrides()[i] * sizeof(Element));
        }

        handle owner;
        if (src.getManager()) {
          owner = PyCapsule_New(
              new ndarray::Manager::Ptr(src.getManager()),
              "ndarray.Manager",
              ndarray::detail::destroyCapsule
          );
        } else {
          owner = handle();
        }

        array_t<Element> result(
          std::move(shape),
          std::move(strides),
          src.getData(),
          owner);

        if (boost::is_const<Element>::value) {
            array_proxy(result.ptr())->flags &= ~detail::npy_api::NPY_ARRAY_WRITEABLE_;
        }

        return result.release();
    }

};

NAMESPACE_END(detail)
NAMESPACE_END(pybind11)
