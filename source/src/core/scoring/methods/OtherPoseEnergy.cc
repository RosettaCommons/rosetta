// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/methods/OtherPoseEnergy.cc
/// @brief Way to handle scoring of a collection of poses, allowing shift of focus between poses.
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/scoring/methods/OtherPoseEnergy.hh>
#include <core/scoring/methods/OtherPoseEnergyCreator.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh>

namespace core {
namespace scoring {
namespace methods {

/// @details This must return a fresh instance of the OtherPoseEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
OtherPoseEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new OtherPoseEnergy );
}

ScoreTypes
OtherPoseEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( other_pose );
	return sts;
}


/// c-tor
OtherPoseEnergy::OtherPoseEnergy() :
	parent( methods::EnergyMethodCreatorOP( new OtherPoseEnergyCreator ) )
{}


/// clone
EnergyMethodOP
OtherPoseEnergy::clone() const
{
	return EnergyMethodOP( new OtherPoseEnergy );
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

void
OtherPoseEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const & scorefxn,
	EnergyMap & totals
) const {

	using namespace pose;
	using namespace pose::full_model_info;

	FullModelInfo & full_model_info = nonconst_full_model_info( pose );
	utility::vector1< PoseOP > const & other_pose_list = full_model_info.other_pose_list();
	if ( other_pose_list.size() == 0 ) return; // no op.

	// watch out for double-counting -- some score terms look at full model, including other poses, already.
	ScoreFunctionOP other_pose_scorefxn = scorefxn.clone();
	other_pose_scorefxn->set_weight( intermol, 0.0 );
	other_pose_scorefxn->set_weight( loop_close, 0.0 );
	other_pose_scorefxn->set_weight( missing_res, 0.0 );
	other_pose_scorefxn->set_weight( free_res, 0.0 );

	totals[ other_pose ] = 0.0;

	for ( Size n = 1; n <= other_pose_list.size(); n++ ){

		PoseOP other_pose_op = other_pose_list[n];
		Real other_pose_score = ( *other_pose_scorefxn )( *other_pose_op );
		totals[ other_pose ] += other_pose_score;

	}

} // finalize_total_energy


core::Size
OtherPoseEnergy::version() const
{
	return 1; // Initial versioning
}

} //methods
} //scoring
} //core
