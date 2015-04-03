// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/comparative_modeling/features/ResidueFeature.hh
/// @brief abstract base class for per-residue features used in comparative
/// modeling.
/// @author James Thompson

#include <core/id/SequenceMapping.hh>
#include <protocols/comparative_modeling/features/ResidueFeature.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace comparative_modeling {
namespace features {

/// @details Auto-generated virtual destructor
ResidueFeature::~ResidueFeature() {}

core::Size ResidueFeature::resnum() const {
	return resnum_;
}

void ResidueFeature::resnum( core::Size const resnum ) {
	resnum_ = resnum;
}

void ResidueFeature::remap( core::id::SequenceMapping const & mapping ) {
	resnum( mapping[ resnum() ] );
}

} // protocols
} // comparative_modeling
} // features
