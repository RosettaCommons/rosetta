// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/pointer/owning_ptr.functions.hh
/// @brief  owning_ptr acquire and release functions
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @remarks
///  @li Objects owned by owning_ptr should inherit from ReferenceCount to get the
///      reference counting mechanism and the add_ref and remove_ref functions.


#ifndef INCLUDED_utility_pointer_refcount_owning_ptr_functions_hh
#define INCLUDED_utility_pointer_refcount_owning_ptr_functions_hh


namespace utility {
namespace pointer {


/// @brief Add a reference to the object acquired by an owning_ptr
template< typename T >
inline
void
owning_ptr_acquire( T * p )
{
	p->add_ref();
}


/// @brief Remove a reference from the object released by an owning_ptr
template< typename T >
inline
void
owning_ptr_release( T * p )
{
	p->remove_ref();
}


} // namespace pointer
} // namespace utility


#endif // INCLUDED_utility_pointer_refcount_owning_ptr_functions_HH
