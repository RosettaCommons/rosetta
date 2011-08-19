#ifndef INCLUDED_ObjexxFCL_rvalue_cast_hh
#define INCLUDED_ObjexxFCL_rvalue_cast_hh


// rvalue_cast: Wrapper for Passing an rvalue as a non-const Reference Argument
//
// Project: Objexx Fortran Compatibility Library (ObjexxFCL)
//
// Version: 3.0.0
//
// Language: C++
//
// Copyright (c) 2000-2009 Objexx Engineering, Inc. All Rights Reserved.
// Use of this source code or any derivative of it is restricted by license.
// Licensing is available from Objexx Engineering, Inc.:  http://objexx.com  Objexx@objexx.com


namespace ObjexxFCL {


/// @brief rvalue_cast: Wrapper for Passing an rvalue as a non-const Reference Argument
template< typename T >
inline
T &
rvalue_cast( T const & t ) {
	return const_cast< T & >( t );
}


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_rvalue_cast_HH
