#ifndef INCLUDED_ObjexxFCL_FArrayTraits_hh
#define INCLUDED_ObjexxFCL_FArrayTraits_hh


// FArrayTraits: FArray Traits Template
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


// ObjexxFCL Headers
#include <ObjexxFCL/FArrayTraits.fwd.hh>


namespace ObjexxFCL {


/// @brief FArrayTraits: FArray Traits Template
/// @note Specialize for types without a default constructor or if a different initial value is desired
template< typename T >
struct FArrayTraits
{
	typedef  T  traits_type;


	/// @brief Initial Value
	inline
	static
	traits_type
	initial_value()
	{
		return traits_type(); // Use default constructor
	}


}; // FArrayTraits


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArrayTraits_HH
