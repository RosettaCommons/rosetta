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
#ifndef NDARRAY_ndarray_h_INCLUDED
#define NDARRAY_ndarray_h_INCLUDED

/** 
 *  @file ndarray.h
 *  @brief Main public header file for ndarray.
 */

#include "ndarray/Array.h"
#include "ndarray/ArrayRef.h"
#include "ndarray/initialization.h"
#ifndef GCC_45
#include "ndarray/operators.h"
#include "ndarray/arange.h"
#endif
#include "ndarray/casts.h"
#include "ndarray/formatting.h"

namespace ndarray {

/** 
 * @mainpage ndarray; Multidimensional Arrays in C++
 *
 * %ndarray is a template library that provides multidimensional array objects in C++, with
 * an interface and features designed to mimic the Python 'numpy' package as much as possible.
 *
 * A tutorial can be found \ref ndarrayTutorial "here".
 *
 * @section comparisons Other Multidimensional Array Libraries
 * A number of other public C++ multidimensional array libraries already exist, most notably
 * boost.MultiArray and Blitz++.  Much of the architecture and some of the interface of the
 * %ndarray templates is built on ideas from both of these, particularly boost.MultiArray.
 * Like Blitz++, it supports shared ownership of data, lazy evaluation of mathematical
 * operations using expression templates, and the ability to use externally allocated memory.
 * Like Boost.MultiArray, %ndarray objects support optimized, natural, nested iteration and
 * preserve const-correctness (albeit with different semantics).
 *
 * @section semantics Copy and Constness Semantics
 * The memory used by %ndarray objects is reference counted, and can be allocated
 * using any STL-compatible allocator.  Arrays can also be constructed from external
 * memory buffers, with full reference counting for any external memory owned by
 * an object that participates in any reference-counting scheme (most notably
 * memory belonging to Python Numpy arrays in C++ Python extensions).  Reference counting
 * can also be disabled for individual arrays that reference external memory that
 * is not reference-counted.
 *
 * Array objects in %ndarray preserve constness in the same way as most C++ smart
 * pointers: there is a distinction between an array with const elements (Array<T const,N>)
 * and a const array value or refernce (Array<T,N> const).  An array with const elements
 * does not support operations that change element values ("deep" operations), while a
 * const array object or reference does not allow the array's data pointer, shape, or strides
 * to be changed ("shallow" operations).  As
 * a result, functions should generally use const references for both input and
 * output array arguments (with const and non-const elements, respectively).
 * Arrays with non-const elements are implicitly (and efficiently) convertible to arrays with
 * const elements, while an array with const elements can be cast to one with non-const elements
 * using the const_array_cast function.
 *
 * @section contiguousness The Row-Major Contiguous Template Parameter
 * In addition to being templated on its element type and total number of dimensions,
 * Array is parameterized by the number of guaranteed row-major contiguous
 * (RMC) dimensions, starting from the end.  This parameter defaults to zero.  For example:
 * <ul>
 * <li> Array<T,2,2> is a fully row-major contiguous matrix.  Each row of the matrix
 *      is a vector with contiguous elements, and there is no space between rows.
 * <li> Array<T,2,1> is a row-major matrix in which each row is a vector with contiguous
 *      elements, but there may be space between the rows (or there may not be; the template
 *      parameter indicates only the guaranteed number of RMC dimensions).
 * <li> Array<T,2,0> is a matrix that can have any strides, including column-major.
 * </ul>
 * An Array with M RMC dimensions can be implicitly converted to an Array with N<=M RMC
 * dimensions.  The static_dimension_cast and dynamic_dimension_cast functions can be used to create
 * arrays with M<N RMC dimensions (dynamic_dimension_cast returns an empty Array if the
 * strides are not appropriately RMC, while static_dimension_cast does no checking).
 *
 * @section ndarrayViews Indexing and Views
 * An Array with dimension N>1 behaves like largely like 
 * an STL container of Array of dimension N-1.  An Array with dimension 1 behaves like a
 * simple container of elements.  Iterators over an Array with N>1 thus dereference
 * to Arrays with dimension N-1, and standard [] indexing also yields
 * the expected lower-dimensional array.
 *
 * Arbitrary views can be retrieved by passing a view definition to an Array's bracket
 * indexing operators.  These view definitions are temporaries created by the view() function:
 * @code
 * // let 'a' be a 3 dimensional array with dimensions (3,5,4)
 *
 * // equivalent to a[1:3,:,2] in numpy:
 * Array<double,2,1> subset1 = a[view(1,3)()(2)];
 *
 * // equivalent to a[:,:,0:4:2] in numpy in numpy:
 * Array<double,3> subset2 = a[view()()(0,4,2)];
 * @endcode
 * Supplying a single integer indexes a dimension by a scalar (which reduces the dimension of
 * output), supplying a pair of integers indicates a contiguous range, and supplying three
 * integers indicates a slice.  Specifying no arguments for a dimension includes the entire
 * dimension.  Not indexing all dimensions is equivalent to including empty parentheses for
 * the remaining dimensions.
 *
 * @section downloads Downloads
 * http://code.google.com/p/ndarray/downloads/list
 *
 * @section License
 * %ndarray is distributed under a simple BSD-like license:
 * http://ndarray.googlecode.com/svn/trunk/COPYING
 */

/** 
 * @page ndarrayTutorial Tutorial
 *
 * @section environment Environment
 * @subsection dependencies Dependencies
 * The versions below indicate the libraries !ndarray is currently developed and tested with.
 * New minor versions should generally work as well, but have not been tested.
 * %ndarray makes extensive use of C++ template metaprogramming, and may
 * not work with older or non-standard-compliant compilers.
 * <ul>
 * <li> Core %ndarray library (ndarray.h): boost 1.38-1.42
 * <li> Python conversion module (ndarray/python.h): python 2.6, numpy 1.2
 * <li> Boost.Python module (ndarray/python/boost/): boost.python 1.38-1.42, python 2.6, numpy 1.2
 * <li> Eigen interface (ndarray/eigen.h): Eigen 2.0
 * <li> Fast Fourier transforms (ndarray/fft.h): FFTW 3.2
 * </ul>
 * 
 * @subsection installation Installation
 * %ndarray is a header-only library; after downloading and unpacking the source,
 * you can start using it immediately simply by adding it to your compiler's include
 * path.
 *
 * For tests, we use the SCons build system, but SCons is not necessary to make use
 * of the library.
 *
 * @section construction Creating New Arrays
 * Array objects have two normal constructors intended for public use, the default
 * constructor and a converting copy constructor.  The default constructor creates
 * an "empty" array, with a null data pointer and zero shape.  The copy constructor
 * creates a shared view into an existing array.
 *
 * @subsection new_arrays New Memory
 * To create a new non-trivial array, one can use the allocate function, which
 * returns a temporary object implicitly convertible to Array:
 * @code
 * Array<double,3,3> a = allocate(makeVector(3,4,5));
 * // equivalent to
 * // >>> a = numpy.empty((3,4,5),dtype=float)
 * // in Python
 * @endcode
 * The makeVector function here is a variable-argument-length constructor for
 * the Vector object, a fixed-size array class whose int variant is
 * used to specify shape and strides for Array.  The appropriate Vector template
 * for a particular Array template is available as the Index typedef within the
 * Array class.
 *
 * The allocate function can also take an STL allocator as a template argument
 * and/or regular argument:
 * @code
 * Array<double,3,3> a = allocate< std::allocator<void> >(makeVector(3,4,5));
 * @endcode
 * @code
 * Array<double,3,3> a = allocate(makeVector(3,4,5), std::allocator<void>());
 * @endcode
 * Note that the type of the allocator does not have to match the type of the
 * array; the allocator's "rebind" member will be used to create the correct
 * allocator when the array is constructed.  Furthermore, unlike standard
 * library containers, Array is not templated on its allocator type; after
 * construction, it is impossible to determine how an Array's memory was
 * allocated.  An Array constructed by allocate is not generally initialized
 * to any value (do not assume it contains zeros).
 * @subsection external_memory External Memory
 * Arrays can also be constructed that point to external data:
 * @code
 * #include <cstdlib>
 * Array<double,1,1>::Owner owner(std::malloc(sizeof(double)*5), std::free);
 * Array<double,1,1>::Index shape = makeVector(5);
 * Array<double,1,1>::Index strides = makeVector(1);
 * Array<double,1,1> a = external(owner.get(), shape, strides, owner);
 * @endcode
 * The 'strides' vector here specifies the space between elements in each
 * dimension; the dot product of the strides vector with an index vector
 * should give the offset of the element with that index from the first
 * element of the array.
 * The 'Owner' type here is a typedef to a boost::shared_ptr, which can
 * take an arbitrary functor as a custom deleter (here, std::free).  By
 * defining an appropriate deleter, an array can manage virtually any kind
 * of memory.  However, it is also possible to create an array with
 * no reference counting by passing an empty owner (or passing none at all):
 * @code
 * #include <cstdlib>
 * double data[] = { 5.3, 1.2, 6.3, 2.8, 7.0 };
 * Array<double,1,1>::Index shape = makeVector(5);
 * Array<double,1,1>::Index strides = makeVector(1);
 * Array<double,1,1> a = external(data, shape, strides);
 * @endcode
 * In this case, the user is responsible for ensuring that the data pointer
 * provided to the array remains valid during the array's lifetime, and
 * is eventually deallocated later.
 *
 * @section assignment Assignment
 * Direct Array assignment is shallow:
 * @code
 * Array<double,1,1> a = allocate(makeVector(5));
 * Array<double,1,1> b;
 * b = a; // the two arrays now share data.
 * @endcode
 * To actually set the elements of an array, we can use Array::deep():
 * @code
 * Array<double,11,> b = allocate(a.getShape());
 * b.deep() = a; // 'b' is now a deep copy of 'a'.
 * @endcode
 * Scalar assignment and augmented assignment operations are also supported:
 * @code
 * b.deep() = 5.2;
 * b.deep() += a;
 * @endcode
 * The deep() method returns a proxy ArrayRef object, which behaves just
 * like an Array aside from its assignment operators, and is
 * implicitly convertible to Array.
 *
 * @section indexing_and_iteration Indexing and Iteration
 * A multidimensional Array behaves like a container of Arrays with lower
 * dimensions, while a one-dimensional Array behaves like a container of
 * elements:
 * @code
 * int data = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
 * Array<int,2,2> a = external(data, makeVector(4,3));
 * // make 'b' a view into the second row of 'a' (Reference is a typedef to Array<int,1,1>):
 * Array<int,2,2>::Reference br = a[1]; // br is an ArrayRef<int,1,1>.
 * Array<int,2,2>::Value b = a[1];      // b is an Array<int,1,1>
 * Array<int,2,2>::Element b0 = b[0];   // b0 == 6; Element is a typedef to int.
 * Array<int,2,2>::Reference::Value b1 = b[1];     // b1 == 7; Reference::Value is also int.
 * Array<int,2,2>::Reference::Reference b2 = b[2]; // b2 == 8; Reference::Reference is int&.
 * @endcode
 * Indexing operations return ArrayRef objects, not Arrays.  This allows them to be
 * assigned to without manually calling deep():
 * @code
 * a[1] -= 3; // subtract three from the entire second row.
 * @endcode
 *
 * For one dimensional arrays, the "Reference" typedef is equivalent to "Element &", while
 * the "Value" typedef is equivalent to "Element".  For multidimensional arrays, "Reference"
 * is a lower-dimensional ArrayRef, while "Value" is a lower-dimensional Array.
 *
 * Array is designed to have lightweight nested iterators, with types provided
 * by the Array::Iterator typedef.  For contiguous one-dimensional arrays (Array<T,1,1>), this
 * is a typedef to a simple pointer.  The typical pattern to iterate over a 3-dimensional array
 * looks like the following:
 * @code
 * Array<double,3,3> a = allocate(makeVector(5,6,8));
 * for (Array<double,3,3>::Iterator i = a.begin(); i != a.end(); ++i) {
 *     for (Array<double,3,3>::Reference::Iterator j = i->begin(); j != i->end(); ++j) {
 *         for (Array<double,3,3>::Reference::Reference::Iterator k = j->begin(); k != j->end(); ++k) {
 *              // *k == a[i - a.begin()][j - i->begin()][k - j->begin()];
 *         }
 *     }
 * }
 * @endcode
 * As expected, the iterators of multidimensional arrays dereference to lower-dimensional arrays,
 * and the iterators of one-dimensional arrays dereference to elements.  With some compilers,
 * it may be advantageous to move the calls to end() outside their loops.
 *
 * Just like direct indexing, multidimensional array iterators dereference to ArrayRef, not Array.
 * 
 * STL-compliant typedefs "iterator", "const_iterator", "reference", "const_reference", and "value"
 * are also provided, though the const variants are not actually const (because Array provides
 * smart-pointer const-correctness rather than container const-correctness).
 *
 * Single elements can be extracted from multidimensional arrays by indexing with ndarray::Vector:
 * @code
 * a[makeVector(3,2,1)] == a[3][2][1];
 * @endcode
 *
 * @section ndarrayTutorialViews Views
 * General views into a an Array are created by passing a ViewDef object to the [] operators
 * of Array, returning a new Array that shares data and owns a reference to the original.
 *
 * ViewDef involves a lot of template metaprogramming, so the actual template class is an
 * internal %detail, and ViewDefs are constructed by calls to the view() function function
 * followed by chained calls to the function call operator, resulting in a syntax that looks
 * like this:
 * @code
 * Array<double,5> a = allocate(makeVector(3,5,2,6,4));
 * b = a[view(1)(0,3)()(4)];
 * @endcode
 * which is equivalent to the Python code:
 * @code
 * a = numpy.empty((3,5,2,6,4),dtype=float)
 * b = a[1,0:3,:,4]
 * @endcode
 * The value passed to each call specifies how to extract values from that dimension:
 * <ul>
 * <li> A single integer selects a single subarray from that dimension, reducing the overall
 *      dimensionality of the view relative to the parent array by one. </li>
 * <li> An empty call selects the entire dimension. </li>
 * <li> A pair of integers selects a contiguous subset of subarrays from that dimension. </li>
 * <li> A triplet of integers selects a noncontiguous subset of subarrays from that dimension. </li>
 * </ul>
 * Any dimensions which are not specified because the length of the ViewDef expression is smaller
 * than the dimensionality of the parent array will be considered full-dimension selections:
 * @code
 * a[view(3)] == a[view(3)()()()()];
 * @endcode
 *
 * @section operators Arithmetic Operators and Comparison
 * Arrays provide element-wise arithmetic operations that make use of expression templates:
 * @code
 * Array<double,2,2> a = allocate(makeVector(3,4));
 * // initialize the elements of 'a'
 * Array<double,2,2> b = allocate(a.getShape());
 * b = a * 3 + 2; // this expression is lazy, and never allocates a temporary array
 * @endcode
 * We can simplify the previous example by initializing 'b' with the copy() function, which is simply
 * a shortcut for allocate and assign:
 * @code
 * Array<double,2,2> a = allocate(makeVector(3,4));
 * // initialize the elements of 'a'
 * Array<double,2,2> b = copy(a * 3 + 2);
 * @endcode
 * As a rule, %ndarray never allocates memory for a new unless you explicitly tell it to with
 * the allocate() or copy() functions.
 *
 * Element-wise comparisons are also supported, but not via overloaded operators:
 * @code
 * Array<double,2,2> a = allocate(makeVector(3,4));
 * // initialize the elements of 'a'
 * Array<double,2,2> b = allocate(makeVector(3,4));
 * // initialize the elements of 'b'
 * Array<bool,2,2> c = copy(equal(a, b));
 * @endcode
 * The element-wise comparison functions (equal, not_equal, less, less_equal, greater, greater_equal)
 * and logical operators (logical_and, logical_or, logical_not) are also lazy, and are most useful
 * when used in conjunction with the reduction functions all() and any():
 * @code
 * Array<double,2,2> a = allocate(makeVector(3,4));
 * // initialize the elements of 'a'
 * Array<double,2,2> b = allocate(makeVector(3,4));
 * // initialize the elements of 'b'
 * bool v1 = all(logical_or(greater(a, b), greater(a, 3.0)));
 * bool v2 = any(less(b, a));
 * @endcode
 * 
 * Array does overload the equality and inequality operators, but these compare for "shallow" equality,
 * not element-wise equality:
 * @code
 * Array<double,2,2> a = allocate(makeVector(3,4));
 * Array<double,2,2> b = copy(a);
 * bool v1 = (a == b); // false, because 'a' and 'b' do not share data.
 * Array<double const,2,1> c(a);
 * bool v2 = (a == c); // true, because 'a' and 'c' have the same data, shape, and strides.
 * @endcode
 */

} // namespace ndarray

#endif // !NDARRAY_ndarray_h_INCLUDED
