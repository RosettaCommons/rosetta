// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/AtomWithDOFChange.hh
/// @brief  Data structure for output-sensitie refold data class declaration
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_core_kinematics_AtomWithDOFChange_hh
#define INCLUDED_core_kinematics_AtomWithDOFChange_hh

// Unit headers
#include <core/kinematics/AtomWithDOFChange.fwd.hh>

/// Project Headers
#include <core/id/AtomID.hh>

namespace core {
namespace kinematics {

/// @brief simple class for use in output-sensitive refold subroutine.
class AtomWithDOFChange
{
public:
	AtomWithDOFChange(
		id::AtomID atomid
	) :
		atomid_( atomid ),
		reached_( false )
	{}

	AtomWithDOFChange() :
		atomid_(),
		reached_( false )
	{}

	id::AtomID atomid_;
	bool   reached_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // namespace kinematics
} // namespace core


#endif // INCLUDED_core_kinematics_AtomWithDOFChange_HH
