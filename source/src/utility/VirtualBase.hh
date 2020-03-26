// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/VirtualBase.hh
/// @brief  A Base class for virtual inheritance hierarchies
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_utility_VirtualBase_hh
#define INCLUDED_utility_VirtualBase_hh

// Unit headers
#include <utility/VirtualBase.fwd.hh>

namespace utility {

/// @brief Base class for Polymorphic classes
/// All classes in Rosetta used polymorphically
/// (that is, where a pointer/reference to a base type can point/reference a derived type)
/// should inherit from this class
///
/// This class serves two purposes:
/// 1) It serves as a common base class for Rosetta object,
///   such that a range of types can be stored as VirtualBaseOPs
///   (and then dynamically cast to the proper type)
/// 2) It ensures that a class hierarchy used polymorphically has a virtual destructor
///   (which is necessary if a derived class can be stored/referenced by a base class pointer/reference).
///
/// @details This class used to be at utility::pointer::ReferenceCount (utility/pointer/ReferenceCount.hh)
/// but was renamed, as it's no longer directly related to pointers/reference counting.
class VirtualBase
{

public: // Creation

	/// @brief Default constructor
	VirtualBase() = default;

	/// @brief The virtual destructor is one of the main reasons for the VirtualBase class
	virtual ~VirtualBase() = default;

	// Because we have a user-defined destructor (to ensure virtual)
	// we need to fill out the rest of the rule-of-five to make sure they exist
	// and aren't implicitly deleted
	VirtualBase( VirtualBase const & ) = default;
	VirtualBase( VirtualBase && ) = default;
	VirtualBase & operator =( VirtualBase const &) = default;
	VirtualBase & operator =( VirtualBase && ) = default;

}; // VirtualBase


} // namespace utility

#endif // INCLUDED_utility_VirtualBase_hh
