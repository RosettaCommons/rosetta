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
#ifndef NDARRAY_initialization_h_INCLUDED
#define NDARRAY_initialization_h_INCLUDED

/** 
 *  \file ndarray/initialization.h @brief Construction functions for array.
 */

#include "ndarray/Array.h"
#include "ndarray/ArrayRef.h"
#include "ndarray/Manager.h"

namespace ndarray {
namespace detail {

struct NullOwner {};

template <int N, typename Derived>
class Initializer {
public:

    template <typename T, int C>
    operator Array<T,N,C> () const {
        return static_cast<Derived const *>(this)->template apply< Array<T,N,C> >();
    }

    template <typename T, int C>
    operator ArrayRef<T,N,C> () const {
        return static_cast<Derived const *>(this)->template apply< ArrayRef<T,N,C> >();
    }

};

template <int N>
class SimpleInitializer : public Initializer< N, SimpleInitializer<N> > {
public:

    template <typename Target>
    Target apply() const {
        typedef detail::ArrayAccess< Target > Access;
        typedef typename Access::Core Core;
        typedef typename Access::Element Element;
        DataOrderEnum order = (ExpressionTraits< Target >::RMC::value < 0) ? COLUMN_MAJOR : ROW_MAJOR;
        Size total = _shape.product();
        std::pair<Manager::Ptr,Element*> p = SimpleManager<Element>::allocate(total);
        return Access::construct(p.second, Core::create(_shape, order, p.first));
    }

    explicit SimpleInitializer(Vector<Size,N> const & shape) : _shape(shape) {}

private:
    Vector<Size,N> _shape;
};

template <typename T, int N, typename Owner>
class ExternalInitializer : public Initializer< N, ExternalInitializer<T,N,Owner> > {
public:

    template <typename Target>
    Target apply() const {
        typedef detail::ArrayAccess< Target > Access;
        typedef typename Access::Core Core;
        Manager::Ptr manager;
        if (!boost::is_same<Owner,NullOwner>::value) {
            manager = makeManager(_owner);
        }
        return Access::construct(_data, Core::create(_shape, _strides, manager));
    }

    ExternalInitializer(
        T * data, 
        Vector<Size,N> const & shape,
        Vector<Offset,N> const & strides,
        Owner const & owner
    ) : _data(data), _owner(owner), _shape(shape), _strides(strides) {}

private:
    T * _data;
    Owner _owner;
    Vector<Size,N> _shape;
    Vector<Offset,N> _strides;
};

} // namespace detail

/// @addtogroup MainGroup
/// @{

/** 
 *  @brief Create an expression that allocates uninitialized memory for an array.
 *
 *  @returns A temporary object convertible to an Array with fully contiguous row-major strides.
 */
template <int N, typename U>
inline detail::SimpleInitializer<N> allocate(Vector<U,N> const & shape) {
    return detail::SimpleInitializer<N>(Vector<Size,N>(shape));
}

/** 
 *  @brief Create an expression that allocates uninitialized memory for a 1-d array.
 *
 *  @returns A temporary object convertible to an Array with fully contiguous row-major strides.
 */
inline detail::SimpleInitializer<1> allocate(Size n) {
    return detail::SimpleInitializer<1>(ndarray::makeVector(n)); 
}

/** 
 *  @brief Create an expression that allocates uninitialized memory for a 2-d array.
 *
 *  @returns A temporary object convertible to an Array with fully contiguous row-major strides.
 */
inline detail::SimpleInitializer<2> allocate(Size n1, Size n2) {
    return detail::SimpleInitializer<2>(ndarray::makeVector(n1, n2)); 
}

/** 
 *  @brief Create an expression that allocates uninitialized memory for a 3-d array.
 *
 *  @returns A temporary object convertible to an Array with fully contiguous row-major strides.
 */
inline detail::SimpleInitializer<3> allocate(Size n1, Size n2, Size n3) {
    return detail::SimpleInitializer<3>(ndarray::makeVector(n1, n2, n3)); 
}

/** 
 *  @brief Create a new Array by copying an Expression.
 */
template <typename Derived>
inline ArrayRef<typename boost::remove_const<typename Derived::Element>::type, 
                Derived::ND::value, Derived::ND::value>
copy(ExpressionBase<Derived> const & expr) {
    ArrayRef<typename boost::remove_const<typename Derived::Element>::type, 
        Derived::ND::value,Derived::ND::value> r(
            allocate(expr.getShape())
        );
    r = expr;
    return r;
}

/// @brief Compute row- or column-major strides for the given shape.
template <int N>
Vector<Offset,N> computeStrides(Vector<Size,N> const & shape, DataOrderEnum order=ROW_MAJOR) {
    Vector<Offset,N> r(1);
    if (order == ROW_MAJOR) {
        for (int n=N-1; n > 0; --n) r[n-1] = r[n] * shape[n];
    } else {
        for (int n=1; n < N; ++n) r[n] = r[n-1] * shape[n-1];
    }
    return r;
}

/** 
 *  @brief Create an expression that initializes an Array with externally allocated memory.
 *
 *  No checking is done to ensure the shape, strides, and data pointers are sensible.
 *
 *  @param[in] data     A raw pointer to the first element of the Array.
 *  @param[in] shape    A Vector of dimensions for the new Array.
 *  @param[in] strides  A Vector of strides for the new Array.
 *  @param[in] owner    A copy-constructable object with an internal reference count
 *                      that owns the memory pointed at by 'data'.
 *
 *  @returns A temporary object convertible to an Array.
 */
template <typename T, int N, typename U, typename V, typename Owner>
inline detail::ExternalInitializer<T,N,Owner> external(
    T * data,
    Vector<U,N> const & shape,
    Vector<V,N> const & strides,
    Owner const & owner
) {
    return detail::ExternalInitializer<T,N,Owner>(
        data,
        Vector<Size,N>(shape.template cast<Size>()),
        Vector<Offset,N>(strides.template cast<Offset>()),
        owner
    );
}

/** 
 *  @brief Create an expression that initializes an Array with externally allocated memory.
 *
 *  No checking is done to ensure the shape, strides, and data pointers are sensible.  Memory will not
 *  be managed at all; the user must ensure the data pointer remains valid for the lifetime of the array.
 *
 *  @param[in] data     A raw pointer to the first element of the Array.
 *  @param[in] shape    A Vector of dimensions for the new Array.
 *  @param[in] strides  A Vector of strides for the new Array.
 *
 *  @returns A temporary object convertible to an Array.
 */
template <typename T, int N, typename U, typename V>
inline detail::ExternalInitializer<T,N,detail::NullOwner> external(
    T * data,
    Vector<U,N> const & shape,
    Vector<V,N> const & strides
) {
    return detail::ExternalInitializer<T,N,detail::NullOwner>(
        data,
        Vector<Size,N>(shape),
        Vector<Offset,N>(strides),
        detail::NullOwner()
    );
}

/** 
 *  @brief Create an expression that initializes an Array with externally allocated memory.
 *
 *  No checking is done to ensure the shape and data pointers are sensible.
 *
 *  @param[in] data     A raw pointer to the first element of the Array.
 *  @param[in] shape    A Vector of dimensions for the new Array.
 *  @param[in] order    Whether the strides are row- or column-major.
 *  @param[in] owner    A copy-constructable object with an internal reference count
 *                      that owns the memory pointed at by 'data'.
 *
 *  @returns A temporary object convertible to an Array.
 */
template <typename T, int N, typename U, typename Owner>
inline detail::ExternalInitializer<T,N,Owner> external(
    T * data,
    Vector<U,N> const & shape,
    DataOrderEnum order,
    Owner const & owner
) {
    return detail::ExternalInitializer<T,N,Owner>(
        data,
        Vector<Size,N>(shape.template cast<Size>()),
        computeStrides(shape.template cast<Size>(), order),
        owner
    );
}

/** 
 *  @brief Create an expression that initializes an Array with externally allocated memory.
 *
 *  No checking is done to ensure the shape and data pointers are sensible.  Memory will not
 *  be managed at all; the user must ensure the data pointer remains valid for the lifetime of the array.
 *
 *  @param[in] data     A raw pointer to the first element of the Array.
 *  @param[in] shape    A Vector of dimensions for the new Array.
 *  @param[in] order    Whether the strides are row- or column-major.
 *
 *  @returns A temporary object convertible to an Array.
 */
template <typename T, int N, typename U>
inline detail::ExternalInitializer<T,N,detail::NullOwner> external(
    T * data,
    Vector<U,N> const & shape,
    DataOrderEnum order = ROW_MAJOR
) {
    return detail::ExternalInitializer<T,N,detail::NullOwner>(
        data,
        Vector<Size,N>(shape.template cast<Size>()),
        computeStrides(shape.template cast<Size>(), order),
        detail::NullOwner()
    );
}

/// @}

template <typename T, int N, int C>
Array<T,N,C>::Array(Size n1, Size n2, Size n3, Size n4, Size n5, Size n6, Size n7, Size n8)
    : Super(0, CorePtr())
{
    typename Super::Index shape;
    if (N > 0) shape[0] = n1;
    if (N > 1) shape[1] = n2;
    if (N > 2) shape[2] = n3;
    if (N > 3) shape[3] = n4;
    if (N > 4) shape[4] = n5;
    if (N > 5) shape[5] = n6;
    if (N > 6) shape[6] = n7;
    if (N > 7) shape[7] = n8;
    this->operator=(ndarray::allocate(shape));
}

template <typename T, int N, int C>
template <typename U>
Array<T,N,C>::Array(Vector<U,N> const & shape)
    : Super(0, CorePtr())
{
    this->operator=(ndarray::allocate(shape.template cast<U>()));
}

template <typename T, int N, int C>
ArrayRef<T,N,C>::ArrayRef(Size n1, Size n2, Size n3, Size n4, Size n5, Size n6, Size n7, Size n8)
    : Super(Array<T,N,C>(n1, n2, n3, n4, n5, n6, n7, n8))
{}

template <typename T, int N, int C>
template <typename U>
ArrayRef<T,N,C>::ArrayRef(Vector<U,N> const & shape)
    : Super(Array<T,N,C>(shape))
{}

} // namespace ndarray

#endif // !NDARRAY_initialization_h_INCLUDED
