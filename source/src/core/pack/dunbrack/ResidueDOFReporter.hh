// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/scoring/dunbrack/ResidueDOFReporter.hh
/// @brief   Class to measure the DOFs used by a RotamerLibrary
/// @author  Andrew Leaver-Fay

#ifndef INCLUDED_core_pack_dunbrack_ResidueDOFReporter_hh
#define INCLUDED_core_pack_dunbrack_ResidueDOFReporter_hh

// Unit Headers
#include <core/pack/dunbrack/ResidueDOFReporter.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/id/PartialAtomID.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// utility headers
#include <utility/VirtualBase.hh>

// C++ headers
#include <set>

namespace core {
namespace pack {
namespace dunbrack {

class ResidueDOFReporter : public utility::VirtualBase {
public:
	ResidueDOFReporter();
	~ResidueDOFReporter() override;

	virtual
	Real
	get_dof( conformation::Residue const & rsd, pose::Pose const & pose ) const = 0;

	virtual
	void
	insert_atoms_defining_dof(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		std::set< id::PartialAtomID > & atoms
	) const = 0;

};

}
}
}


#endif
