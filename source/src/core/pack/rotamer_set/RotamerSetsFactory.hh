// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rotamer_set/RotamerSetsFactory.hh
/// @brief  (Symmetry agnostic) factory for the RotamerSets class
/// @author Rocco Moretti (rmorettase@gmail.com)

#ifndef INCLUDED_core_pack_rotamer_set_RotamerSetsFactory_hh
#define INCLUDED_core_pack_rotamer_set_RotamerSetsFactory_hh

// Unit headers
#include <core/pack/rotamer_set/RotamerSetsFactory.fwd.hh>

// Package headers
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>


// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace rotamer_set {

class RotamerSetsFactory : public utility::pointer::ReferenceCount
{
public:
	virtual ~RotamerSetsFactory();

	/// @brief Create a generic RotamerSet object for the pose.
	static RotamerSetsOP create_rotamer_sets( core::pose::Pose const & pose );
};

}
}
}

#endif
