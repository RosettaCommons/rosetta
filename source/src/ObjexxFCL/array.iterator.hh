#ifndef INCLUDED_ObjexxFCL_array_iterator_hh
#define INCLUDED_ObjexxFCL_array_iterator_hh


// C Array Iterator Functions
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


// C++ Headers
#include <cstddef>


namespace ObjexxFCL {


/// @brief Begin Iterator for C Array
template < typename T, std::size_t N >
inline
T *
begin( T (&array)[N] )
{
	return array + 0;
}


/// @brief End Iterator for C Array
template < typename T, std::size_t N >
inline
T *
end( T (&array)[N] )
{
	return array + N;
}


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_array_iterator_HH
