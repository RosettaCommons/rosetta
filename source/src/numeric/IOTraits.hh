// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/IOTraits.hh
/// @brief  Numerics input/output type traits
/// @author Frank M. D'Ippolito (Objexx@objexx.com)
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_numeric_IOTraits_hh
#define INCLUDED_numeric_IOTraits_hh


namespace numeric {


/// @brief Numerics input/output type traits
template< typename T >
struct IOTraits
{
	typedef T Type;


	/// @brief Precision
	inline
	static
	int
	precision()
	{
		return 0; // No precision for generic types
	}


	/// @brief Width
	inline
	static
	int
	width()
	{
		return 0; // No minimum width for generic types
	}


}; // IOTraits


/// @brief Numerics input/output type traits short int specialization
template<>
struct IOTraits< short int >
{
	typedef short int Type;


	/// @brief Precision
	inline
	static
	int
	precision()
	{
		return 0;
	}


	/// @brief Width
	inline
	static
	int
	width()
	{
		return 7;
	}


}; // IOTraits


/// @brief Numerics input/output type traits int specialization
template<>
struct IOTraits< int >
{
	typedef int Type;


	/// @brief Precision
	inline
	static
	int
	precision()
	{
		return 0;
	}


	/// @brief Width
	inline
	static
	int
	width()
	{
		return 12;
	}


}; // IOTraits


/// @brief: Numerics input/output type traits long int specialization
template<>
struct IOTraits< long int >
{
	typedef long int Type;


	/// @brief Precision
	inline
	static
	int
	precision()
	{
		return 0;
	}


	/// @brief Width
	inline
	static
	int
	width()
	{
		return 22; // Big enough for 64-bit LP64 representation
	}


}; // IOTraits


/// @brief: Numerics input/output type traits unsigned short int specialization
template<>
struct IOTraits< unsigned short int >
{
	typedef unsigned short int Type;


	/// @brief Precision
	inline
	static
	int
	precision()
	{
		return 0;
	}


	/// @brief Width
	inline
	static
	int
	width()
	{
		return 7;
	}


}; // IOTraits


/// @brief: Numerics input/output type traits unsigned int specialization
template<>
struct IOTraits< unsigned int >
{
	typedef unsigned int Type;


	/// @brief Precision
	inline
	static
	int
	precision()
	{
		return 0;
	}


	/// @brief Width
	inline
	static
	int
	width()
	{
		return 12;
	}


}; // IOTraits


/// @brief Numerics input/output type traits unsigned long int specialization
template<>
struct IOTraits< unsigned long int >
{
	typedef unsigned long int Type;


	/// @brief Precision
	inline
	static
	int
	precision()
	{
		return 0;
	}


	/// @brief Width
	inline
	static
	int
	width()
	{
		return 22; // Big enough for 64-bit LP64 representation
	}


}; // IOTraits


/// @brief Numerics input/output type traits float Specialization
template<>
struct IOTraits< float >
{
	typedef float Type;


	/// @brief Precision
	inline
	static
	int
	precision()
	{
		return 8;
	}


	/// @brief Width
	inline
	static
	int
	width()
	{
		return 15;
	}


}; // IOTraits


/// @brief Numerics input/output type traits double specialization
template<>
struct IOTraits< double >
{
	typedef double Type;


	/// @brief Precision
	inline
	static
	int
	precision()
	{
		return 16;
	}


	/// @brief Width
	inline
	static
	int
	width()
	{
		return 23;
	}


}; // IOTraits


/// @brief Numerics input/output type traits long double specialization
template<>
struct IOTraits< long double >
{
	typedef long double Type;


	/// @brief Precision
	inline
	static
	int
	precision()
	{
		return 33; // Big enough for 128-bit representation
	}


	/// @brief Width
	inline
	static
	int
	width()
	{
		return 40; // Big enough for 128-bit representation
	}


}; // IOTraits


} // namespace numeric


#endif // INCLUDED_numeric_IOTraits_HH
