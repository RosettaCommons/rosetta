#ifndef INCLUDED_ObjexxFCL_ChunkExponent_hh
#define INCLUDED_ObjexxFCL_ChunkExponent_hh


// ChunkExponent: ChunkVector Exponent Wrapper for Function Disambiguation and Range Clipping
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
#include <algorithm>
#include <cstddef>
#include <limits>


namespace ObjexxFCL {


/// @brief ChunkExponent: ChunkVector Exponent Wrapper for Function Disambiguation and Range Clipping
///
/// @remarks
///  The exponent is clipped to be less than the number of bits in its type so that 2^exponent can
///   be stored in that type
class ChunkExponent
{


public: // Types


	typedef  std::size_t  T;
	typedef  T  value_type;


public: // Creation


	/// @brief Constructor (Implicit): Clips Exponent to Valid Range
	inline
	ChunkExponent( T const exponent_a ) :
		exponent_( std::min( exponent_a, static_cast< T >( std::numeric_limits< T >::digits - 1 ) ) )
	{}


	/// @brief Destructor
	inline
	~ChunkExponent()
	{}


public: // Conversion


	/// @brief Exponent Value Conversion
	inline
	operator T() const
	{
		return exponent_;
	}


private: // Data


	/// @brief Exponent value
	T exponent_;


}; // ChunkExponent


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_ChunkExponent_HH
