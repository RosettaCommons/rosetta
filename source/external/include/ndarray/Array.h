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
#ifndef NDARRAY_Array_h_INCLUDED
#define NDARRAY_Array_h_INCLUDED

/** 
 *  @file ndarray/Array.h
 *
 *  @brief Definitions for Array.
 */

#include "ndarray_fwd.h"
#include "ndarray/ArrayTraits.h"
#include "ndarray/ArrayBaseN.h"
#include "ndarray/Vector.h"
#include "ndarray/detail/Core.h"
#include "ndarray/views.h"

namespace ndarray {

/**
 *  @brief A multidimensional strided array.
 *
 *  Array is the workhorse class of the ndarray library.
 */
template <typename T, int N, int C>
class Array : public ArrayBaseN< Array<T,N,C> > {
    typedef ArrayBaseN<Array> Super;
    typedef typename Super::Core Core;
    typedef typename Super::CorePtr CorePtr;
public:

    /** 
     *  @brief Default constructor. 
     *
     *  Creates an empty array with zero dimensions and null memory.
     */
    Array() : Super(0, Core::create()) {}

    /**
     *  @brief Non-converting copy constructor.
     */
    Array(Array const & other) : Super(other._data, other._core) {}

    /**
     *  @brief Converting copy constructor. 
     *
     *  Implicit conversion is allowed for non-const to const and for 
     *  more guaranteed RMC to less guaranteed RMC (see \ref index).
     */
    template <typename T_, int C_>
    Array(
        Array<T_,N,C_> const & other
#ifndef DOXYGEN
        , typename boost::enable_if<detail::Convertible<N,T_,C_,T,C>,void*>::type=0
#endif
    ) : Super(other._data, other._core) {}

    /**
     *  @brief Converting copy constructor. 
     *
     *  Implicit conversion is allowed for non-const to const and for 
     *  more guaranteed RMC to less guaranteed RMC (see \ref index).
     */
    template <typename T_, int C_>
    Array(
        ArrayRef<T_,N,C_> const & other
#ifndef DOXYGEN
        , typename boost::enable_if<detail::Convertible<N,T_,C_,T,C>,void*>::type=0
#endif
    ) : Super(other._data, other._core) {}

    /**
     *  @brief Construct an array with the given dimensions and allocated but uninitialized memory.
     *
     *  Unspecified dimensions will have unit size, and if the number of argmuments is greater
     *  than the number of dimensions of the array, the extra arguments will be silently ignored.
     *
     *  This is implemented in initialization.h.
     */
    explicit Array(Size n1, Size n2=1, Size n3=1, Size n4=1, Size n5=1, Size n6=1, Size n7=1, Size n8=1);

    /**
     *  @brief Construct an array with the given dimensions and allocated but uninitialized memory.
     *
     *  This is implemented in initialization.h.
     */
    template <typename U>
    explicit Array(Vector<U,N> const & shape);

    /**
     *  @brief Non-converting shallow assignment.
     */
    Array & operator=(Array const & other) {
        if (&other != this) {
            this->_data = other._data;
            this->_core = other._core;
        }
        return *this;
    }

    /**
     *  @brief Converting shallow assignment. 
     *
     *  Implicit conversion is allowed for non-const -> const and for 
     *  more guaranteed RMC -> less guaranteed RMC (see \ref index).
     */
    template <typename T_, int C_>
#ifndef DOXYGEN
    typename boost::enable_if<detail::Convertible<N,T_,C_,T,C>, Array &>::type
#else
    Array &
#endif
    operator=(Array<T_,N,C_> const & other) {
        this->_data = other._data;
        this->_core = other._core;
        return *this;
    }

    /**
     *  @brief Converting shallow assignment. 
     *
     *  Implicit conversion is allowed for non-const -> const and for 
     *  more guaranteed RMC -> less guaranteed RMC (see \ref index).
     */
    template <typename T_, int C_>
#ifndef DOXYGEN
    typename boost::enable_if<detail::Convertible<N,T_,C_,T,C>, Array &>::type
#else
    Array &
#endif
    operator=(ArrayRef<T_,N,C_> const & other) {
        this->_data = other._data;
        this->_core = other._core;
        return *this;
    }

    /**
     *  @brief Shallow equality comparison: return true if the arrays share data and
     *         have the same shape and strides. 
     */
    template <typename T_, int N_, int C_>
    bool operator==(Array<T_,N_,C_> const & other) const {
        return this->getData() == other.getData()
            && this->getShape() == other.getShape()
            && this->getStrides() == other.getStrides();
    }

    /**
     *  @brief Shallow inequality comparison. 
     */
    template <typename T_, int N_, int C_>
    bool operator!=(Array<T_,N_,C_> const & other) const {
        return !this->operator==(other);
    }

    /// @brief Lightweight shallow swap.
    void swap(Array & other) {
        std::swap(this->_data, other._data);
        this->_core.swap(other._core);
    }

    /**
     *  @brief Return true if the Array is definitely unique.
     *
     *  This will only return true if the manager overrides Manager::isUnique();
     *  this is true for the SimpleManager used by ndarray::allocate, but it is
     *  not true for ExternalManager.
     */
    bool isUnique() const { return this->_core->isUnique(); }

private:
    template <typename T_, int N_, int C_> friend class Array;
    template <typename T_, int N_, int C_> friend class ArrayRef;
    template <typename T_, int N_, int C_> friend struct ArrayTraits;
    template <typename Derived> friend class ArrayBase;
    template <typename Array_> friend class detail::ArrayAccess;

    /// @internal @brief Construct an Array from a pointer and Core.
    Array(T * data, CorePtr const & core) : Super(data, core) {}
};

} // namespace ndarray

#endif // !NDARRAY_Array_h_INCLUDED
