// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/evaluation/StructuralSimilarityEvaluator.cc
/// @brief
/// @author James Thompson

#include <core/types.hh>
// AUTO-REMOVED #include <core/pose/util.hh>
#include <core/pose/Pose.hh>

#include <core/scoring/rms_util.hh>

#include <utility/vector1.hh>
#include <numeric/util.hh>

#include <protocols/evaluation/StructuralSimilarityEvaluator.hh>

#include <string>

//Auto Headers
#include <core/conformation/Residue.hh>
namespace protocols  {
namespace evaluation {

StructuralSimilarityEvaluator::StructuralSimilarityEvaluator(
	utility::vector1< core::pose::Pose > const & poses,
	std::string const & atom_name,
	std::string const & tag
) :
	SingleValuePoseEvaluator< core::Real >(tag),
	poses_(poses),
	atom_name_(atom_name)
{}

StructuralSimilarityEvaluator::~StructuralSimilarityEvaluator() {}

void StructuralSimilarityEvaluator::apply(
	core::pose::Pose & pose,
	std::string /*tag*/,
	core::io::silent::SilentStruct & ss
) const {
	using core::Size;
	using core::Real;
	using core::pose::Pose;
	using utility::vector1;

	typedef vector1< Pose >::const_iterator iter;

	vector1< Real > gdtmms( poses_.size(), 0.0 );;
	for ( Size ii = 1; ii <= poses_.size(); ++ii ) {
		gdtmms[ii] = core::scoring::CA_gdtmm( pose, poses_[ii] );
	}

	Real const median_gdtmm( numeric::median( gdtmms ) );
	ss.add_energy( "median_sim", median_gdtmm );
}

} // evaluation
} // protocols
