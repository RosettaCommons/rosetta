// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/vectorL_Selector.hh
/// @brief  vectorL lower index selectors
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_vectorL_Selector_hh
#define INCLUDED_utility_vectorL_Selector_hh


// Unit headers
#include <utility/vectorL.fwd.hh>

// C++ headers
#include <cstddef>


namespace utility {


/// @brief vectorL index type selector: Nonnegative lower index non-specialization
template< bool >
struct vectorL_IndexSelector
{
	typedef  std::size_t  index_type;
	typedef  std::size_t  Index;
};


/// @brief vectorL index type selector: Negative lower index specialization
template<>
struct vectorL_IndexSelector< false >
{
	typedef  ssize_t  index_type;
	typedef  ssize_t  Index;
};


/// @brief vectorL lower index zero selector: Nonzero lower index non-specialization
template< bool >
struct vectorL_ZeroSelector
{


	/// @brief ( i >= l )
	inline
	static
	bool
	ge( std::size_t const i, std::size_t const l )
	{
		return ( i >= l );
	}


	/// @brief ( i >= l ) and ( i - l ) is within ssize_t range
	inline
	static
	bool
	ge( ssize_t const i, ssize_t const l )
	{
		return ( ( i >= l ) && ( i - l >= 0 ) ); // Second clause catches range errors for the difference
	}


};


/// @brief vectorL lower index zero selector: Zero lower index specialization
template<>
struct vectorL_ZeroSelector< false >
{


	/// @brief ( i >= l )
	inline
	static
	bool
	ge( std::size_t const, std::size_t const )
	{
		return true; // Always true when l == 0
	}


};


} // namespace utility


#endif // INCLUDED_utility_vectorL_Selector_HH
