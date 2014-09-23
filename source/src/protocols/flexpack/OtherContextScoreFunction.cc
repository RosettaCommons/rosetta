// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/flexpack/OtherContextScoreFunction.cc
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), Florian Richter (floric@u.washington.edu)

// Unit headers
#include <protocols/flexpack/OtherContextScoreFunction.hh>

// Project headers
#include <core/pose/Pose.hh>  // #INCLUDING POSE.HH IN CC FILES ONLY

#include <utility/vector1.hh>


namespace protocols {
namespace flexpack {


OtherContextScoreFunction::OtherContextScoreFunction() :
	ScoreFunction(),
	scored_context_pose_( false )
{
}

OtherContextScoreFunction::~OtherContextScoreFunction() {}

OtherContextScoreFunction::OtherContextScoreFunction(
	core::pose::Pose const & context_pose
) :
	ScoreFunction(),
	context_pose_( core::pose::PoseOP( new core::pose::Pose( context_pose ) ) ),
	scored_context_pose_( false )
{
}

void
OtherContextScoreFunction::pre_scoring()
{
	(*this)( *context_pose_ ); // score the context pose.
}

void
OtherContextScoreFunction::set_context_pose( core::pose::Pose const & pose )
{
	context_pose_ = core::pose::PoseOP( new core::pose::Pose( pose ) );
	scored_context_pose_ = false;
}

void
OtherContextScoreFunction::eval_cd_1b(
	core::conformation::Residue const & rsd,
	core::pose::Pose const & /*pose*/,
	core::scoring::EnergyMap & emap
) const
{
	parent::eval_cd_1b( rsd, *context_pose_, emap );
}


void
OtherContextScoreFunction::eval_cd_2b(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	core::pose::Pose const & /*pose*/,
	core::scoring::EnergyMap & emap
) const
{
	parent::eval_cd_2b( rsd1, rsd2, *context_pose_, emap );
}


void
OtherContextScoreFunction::eval_cd_intrares_energy(
	core::conformation::Residue const & rsd,
	core::pose::Pose const & /*pose*/,
	core::scoring::EnergyMap & emap
) const
{
	parent::eval_cd_intrares_energy( rsd, *context_pose_, emap );
}


}
}

