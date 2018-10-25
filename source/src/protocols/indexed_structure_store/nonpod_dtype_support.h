// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/indexed_structure_store/nonpod_dtype_support.hh
/// @brief Override-template allowing numpy-access to non-pod data structures.
/// @author Alex Ford (fordas@uw.edu)
//
#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

#pragma once

namespace pybind11 { namespace detail {

template <typename T> struct unsafe_nonpod_npy_format_descriptor {
    static constexpr auto name = make_caster<T>::name;

    static pybind11::dtype dtype() {
        return reinterpret_borrow<pybind11::dtype>(dtype_ptr());
    }

    static std::string format() {
        static auto format_str = get_numpy_internals().get_type_info<T>(true)->format_str;
        return format_str;
    }

    static void register_dtype(any_container<field_descriptor> fields) {
        register_structured_dtype(std::move(fields), typeid(typename std::remove_cv<T>::type),
                                  sizeof(T), &direct_converter);
    }

private:
    static PyObject* dtype_ptr() {
        static PyObject* ptr = get_numpy_internals().get_type_info<T>(true)->dtype_ptr;
        return ptr;
    }

    static bool direct_converter(PyObject *obj, void*& value) {
        auto& api = npy_api::get();
        if (!PyObject_TypeCheck(obj, api.PyVoidArrType_Type_))
            return false;
        if (auto descr = reinterpret_steal<object>(api.PyArray_DescrFromScalar_(obj))) {
            if (api.PyArray_EquivTypes_(dtype_ptr(), descr.ptr())) {
                value = ((PyVoidScalarObject_Proxy *) obj)->obval;
                return true;
            }
        }
        return false;
    }
};

}}

#define DECLARE_NONPOD_DTYPE_SUPPORT(TYPE_NAME) \
namespace pybind11 { namespace detail { \
    template <> struct npy_format_descriptor<TYPE_NAME, void> : \
      unsafe_nonpod_npy_format_descriptor<TYPE_NAME> { \
    }; \
  } \
 \
  template <> struct format_descriptor<TYPE_NAME, void> { \
      static std::string format() { \
          return detail::npy_format_descriptor<typename std::remove_cv<TYPE_NAME>::type>::format(); \
      } \
  }; \
}
