// -*- c++ -*-
/*
 * Copyright (c) 2010-2012, Jim Bosch
 * All rights reserved.
 *
 * ndarray is distributed under a simple BSD-like license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/ndarray/ndarray
 */
#ifndef NDARRAY_EIGEN_BP_ReturnInternal_h_INCLUDED
#define NDARRAY_EIGEN_BP_ReturnInternal_h_INCLUDED

#if defined __GNUC__ && __GNUC__>=6
 #pragma GCC diagnostic ignored "-Wignored-attributes"
 #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

#include "boost/numpy.hpp"
#include <Eigen/Core>

namespace ndarray { namespace detail {

template <typename Input>
struct ShallowToBoostPython {

    typedef boost::numpy::ndarray result_type;

    typedef typename Input::Scalar Scalar;

    static boost::numpy::ndarray apply(Input const & input, bool writeable=true) {
        boost::numpy::dtype dtype = boost::numpy::dtype::get_builtin<Scalar>();
        int itemsize = dtype.get_itemsize();
        std::vector<Py_intptr_t> shape(2);
        shape[0] = input.rows();
        shape[1] = input.cols();
        std::vector<Py_intptr_t> strides(2, itemsize);
        strides[0] *= input.rowStride();
        strides[1] *= input.colStride();
        boost::python::object none;
        boost::numpy::ndarray result = (writeable) ?
            boost::numpy::from_data(const_cast<Scalar const*>(input.data()), dtype, shape, strides, none) :
            boost::numpy::from_data(const_cast<Scalar*>(input.data()), dtype, shape, strides, none);
        if (Input::RowsAtCompileTime == 1 || Input::ColsAtCompileTime == 1) return result.squeeze();
        return result;
    }
};

} // namespace detail

template <std::size_t owner_arg = 1, typename Base = boost::python::default_call_policies>
struct ReturnInternal : public Base {

    /**
     *  \brief A model of ResultConverter that creates a numpy.matrix view into an Eigen object.
     *
     *  Note that the returned view does not own its own data, and doesn't (yet) own a reference
     *  to the object that does.  For this reason, this result converter should ONLY be used by
     *  return_internal_matrix or with another model of CallPolicies with a postcall operations
     *  that installs a reference to the owning object.
     */
    template <typename Input, bool writeable>
    struct MakeShallowMatrix {

        inline bool convertible() const { return true; }

        inline PyObject * operator()(Input const & x) const {
            try {
                boost::python::object r = detail::ShallowToBoostPython< Input >::apply(x, writeable);
                Py_INCREF(r.ptr());
                return r.ptr();
            } catch (boost::python::error_already_set & exc) {
                boost::python::handle_exception();
                return NULL;
            }
        }

        inline PyTypeObject const * get_pytype() const {
            return boost::python::converter::object_manager_traits<
                typename detail::ShallowToBoostPython< Input >::result_type
            >::get_pytype();
        }

    };

    /// \brief Model of ResultConverterGenerator for make_shallow_matrix.
    struct result_converter {
        template <typename T>
        struct apply {
            typedef typename boost::remove_reference<T>::type dereferenced;
            
            typedef MakeShallowMatrix<
                typename boost::remove_const<dereferenced>::type,
                boost::is_const<dereferenced>::value
            > type;
        };
    };

    template <typename ArgumentPackage>
    static PyObject * postcall(ArgumentPackage const & args_, PyObject * result) {
        std::size_t arity_ = boost::python::detail::arity(args_);
        if (owner_arg > arity_ || owner_arg < 1) {
            PyErr_SetString(PyExc_IndexError,
                            "ndarray::ReturnInternal: argument out of range.");
            return NULL;
        }
        PyObject * owner_ref = boost::python::detail::get_prev<owner_arg>::execute(args_, result);
        if (owner_ref == NULL) return 0;
        result = Base::postcall(args_, result);
        if (result == NULL) return 0;
        try {
            boost::python::object owner(
                reinterpret_cast<boost::python::detail::borrowed_reference>(owner_ref)
            );
            boost::numpy::ndarray array(
                reinterpret_cast<boost::python::detail::borrowed_reference>(result)
            );
            array.set_base(owner);
            Py_INCREF(array.ptr());
            result = array.ptr();
        } catch (boost::python::error_already_set & err) {
            boost::python::handle_exception();
            return NULL;
        }
        return result;
    }
};

} // namespace ndarray

#endif // !NDARRAY_EIGEN_BP_ReturnInternal_h_INCLUDED
