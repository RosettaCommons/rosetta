// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file utility/cxx_versioning_macros.hh
/// @brief macros for using function attributes in a portable way
/// @author Jack Maguire, jackmaguire1444@gmail.com

/* __has_cpp_attribute
__has_cpp_attribute is a standard C++ macro for testing if a given attribute is available or not.
Just in case we are dealing with an old enough compiler, or one that isn't standards conformant, make sure that __has_cpp_attribute is defined.
Technically, __has_cpp_attribute is only required by C++20, but both GCC and Clang supported it much earlier. (Certainly by the time of C++17)
*/
#ifndef __has_cpp_attribute
# define __has_cpp_attribute( x ) 0
#endif

/* NODISCARD_ATTR
If a function declared nodiscard or a function returning an enumeration or class declared nodiscard by value is called from a discarded-value expression other than a cast to void, the compiler is encouraged to issue a warning.

*/
#if __has_cpp_attribute(nodiscard) && CXX_17_OR_LATER
#define NODISCARD_ATTR [[nodiscard]]
#else
#define NODISCARD_ATTR
#endif

/* UNUSED_ATTR
This attribute, attached to a function, means that the function is meant to be possibly unused. GCC does not produce a warning for this function.
*/
#if __has_cpp_attribute(unused) && CXX_17_OR_LATER
#define UNUSED_ATTR [[unused]]
#else
#define UNUSED_ATTR
#endif

/* PURE_ATTR
Calls to functions that have no observable effects on the state of the program other than to return a value may lend themselves to optimizations such as common subexpression elimination. Declaring such functions with the pure attribute allows GCC to avoid emitting some calls in repeated invocations of the function with the same argument values.

The pure attribute prohibits a function from modifying the state of the program that is observable by means other than inspecting the function’s return value. However, functions declared with the pure attribute can safely read any non-volatile objects, and modify the value of objects in a way that does not affect their return value or the observable state of the program.

For example,

int hash (char *) __attribute__ ((pure));
tells GCC that subsequent calls to the function hash with the same string can be replaced by the result of the first call provided the state of the program observable by hash, including the contents of the array itself, does not change in between. Even though hash takes a non-const pointer argument it must not modify the array it points to, or any other object whose value the rest of the program may depend on. However, the caller may safely change the contents of the array between successive calls to the function (doing so disables the optimization). The restriction also applies to member objects referenced by the this pointer in C++ non-static member functions.

Some common examples of pure functions are strlen or memcmp. Interesting non-pure functions are functions with infinite loops or those depending on volatile memory or other system resource, that may change between consecutive calls (such as the standard C feof function in a multithreading environment).

The pure attribute imposes similar but looser restrictions on a function’s definition than the const attribute: pure allows the function to read any non-volatile memory, even if it changes in between successive invocations of the function. Declaring the same function with both the pure and the const attribute is diagnosed. Because a pure function cannot have any observable side effects it does not make sense for such a function to return void. Declaring such a function is diagnosed.
*/
#if __has_cpp_attribute(pure) && CXX_17_OR_LATER
#define PURE_ATTR [[pure]]
#else
#define PURE_ATTR
#endif

/* CONST_ATTR
Calls to functions whose return value is not affected by changes to the observable state of the program and that have no observable effects on such state other than to return a value may lend themselves to optimizations such as common subexpression elimination. Declaring such functions with the const attribute allows GCC to avoid emitting some calls in repeated invocations of the function with the same argument values.

For example,
`int square (int) __attribute__ ((const));`
tells GCC that subsequent calls to function square with the same argument value can be replaced by the result of the first call regardless of the statements in between.

The const attribute prohibits a function from reading objects that affect its return value between successive invocations. However, functions declared with the attribute can safely read objects that do not change their return value, such as non-volatile constants.

The const attribute imposes greater restrictions on a function’s definition than the similar pure attribute. Declaring the same function with both the const and the pure attribute is diagnosed. Because a const function cannot have any observable side effects it does not make sense for it to return void. Declaring such a function is diagnosed.

Note that a function that has pointer arguments and examines the data pointed to must not be declared const if the pointed-to data might change between successive invocations of the function. In general, since a function cannot distinguish data that might change from data that cannot, const functions should never take pointer or, in C++, reference arguments. Likewise, a function that calls a non-const function usually must not be const itself.

Since "const" is also a reserved keyword, it doesn't work with the __has_cpp_attribute test macro, so instead we'll enable it if the compiler is running C++17 or later (which requires that unrecognized attributes do not raise errors).
*/
#if CXX_17_OR_LATER
#define CONST_ATTR [[const]]
#else
#define CONST_ATTR
#endif

/* NORETURN_ATTR
This attribute, attached to a function, means that the function never returns - this could be because it always throws an exception, aborts the current thread or immediately exits the current process.

This is defined by C++11, so we don't need to put it conditionally. Though we do use a macro to define it, as the beautifier doesn't like it otherwise.
*/
#define NORETURN_ATTR [[noreturn]]
