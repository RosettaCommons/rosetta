// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/pointer/pointer.functions.hh
/// @brief  Pointer utility functions
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_pointer_pointer_functions_hh
#define INCLUDED_utility_pointer_pointer_functions_hh


namespace utility {
namespace pointer {


///@brief  Delete a pointer and assign a new pointer
template< typename T >
inline
void
delete_and_assign( T * & p, T * n )
{
	delete p; p = n;
}


///@brief  Delete a pointer and assign and set it to zero
template< typename T >
inline
void
delete_and_zero( T * & p )
{
	delete p; p = 0;
}


///@brief  Delete a pointer to an array and assign a new pointer
template< typename T >
inline
void
delete_and_assign_array( T * & p, T * n )
{
	delete[] p; p = n;
}


///@brief  Delete a pointer to an array and set it to zero
template< typename T >
inline
void
delete_and_zero_array( T * & p )
{
	delete[] p; p = 0;
}


} // namespace pointer
} // namespace utility


#endif // INCLUDED_utility_pointer_pointer_functions_HH
