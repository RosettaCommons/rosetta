// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/utilities/PrepareForCentroid.hh>
#include <protocols/loop_modeling/utilities/PrepareForCentroidCreator.hh>

// Core headers
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/pose/Pose.hh>
#include <core/util/SwitchResidueTypeSet.hh>

// Protocol headers
#include <protocols/moves/Mover.hh>

namespace protocols {
namespace loop_modeling {
namespace utilities {

moves::MoverOP PrepareForCentroidCreator::create_mover() const {
	return moves::MoverOP( new PrepareForCentroid );
}

string PrepareForCentroidCreator::keyname() const {
	return "PrepareForCentroid";
}

PrepareForCentroid::PrepareForCentroid() {}

bool PrepareForCentroid::do_apply(Pose & pose) {
	using core::util::switch_to_residue_type_set;
	using core::chemical::CENTROID;

	if ( ! pose.is_centroid() ) {
		switch_to_residue_type_set(pose, CENTROID);
	}

	return true;
}

}
}
}

