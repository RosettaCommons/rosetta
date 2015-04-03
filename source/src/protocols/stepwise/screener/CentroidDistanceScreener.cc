// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/CentroidDistanceScreener.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/CentroidDistanceScreener.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSampler.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerWithResidueList.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerWithResidueAlternatives.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/rna/util.hh>
#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "protocols.stepwise.screener.CentroidDistanceScreener" );

using namespace core;
using namespace protocols::stepwise::sampler;
using namespace protocols::stepwise::sampler::rigid_body;

namespace protocols {
namespace stepwise {
namespace screener {

	//Constructor
	CentroidDistanceScreener::CentroidDistanceScreener(  pose::Pose & screening_pose,
																											 Size const moving_res,
																											 Vector const & reference_centroid,
																											 core::Real const max_distance_squared ):
		screening_pose_( screening_pose ),
		moving_res_( moving_res ),
		reference_centroid_( reference_centroid ),
		max_distance_squared_( max_distance_squared )
	{}

	//Destructor
	CentroidDistanceScreener::~CentroidDistanceScreener()
	{}


	//////////////////////////////////////////
	bool
	CentroidDistanceScreener::check_screen() {
		Vector const moving_centroid = core::chemical::rna::get_rna_base_centroid( screening_pose_.residue( moving_res_ ) );
		return ( ( moving_centroid - reference_centroid_ ).length_squared() <= max_distance_squared_ );
	}

	//////////////////////////////////////////
	void
	CentroidDistanceScreener::fast_forward( sampler::StepWiseSamplerBaseOP sampler ){
		if ( sampler->type() == RIGID_BODY_WITH_RESIDUE_LIST ){
			RigidBodyStepWiseSamplerWithResidueList & rigid_body_rotamer_with_copy_dofs = *( static_cast< RigidBodyStepWiseSamplerWithResidueList * >( sampler.get() ) );
			rigid_body_rotamer_with_copy_dofs.fast_forward_to_next_translation();
		} else if ( sampler->type() == RIGID_BODY_WITH_RESIDUE_ALTERNATIVES ){
			RigidBodyStepWiseSamplerWithResidueAlternatives & rigid_body_rotamer_with_residue_alternatives = *( static_cast< RigidBodyStepWiseSamplerWithResidueAlternatives * >( sampler.get() ) );
			rigid_body_rotamer_with_residue_alternatives.fast_forward_to_next_translation();
		} else if ( sampler->type() == RIGID_BODY ){
			RigidBodyStepWiseSampler & rigid_body_rotamer = *( static_cast< RigidBodyStepWiseSampler * >( sampler.get() ) );
			rigid_body_rotamer.fast_forward_to_next_translation();
		}
	}


} //screener
} //stepwise
} //protocols
