// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/BaseCentroidScreener.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/BaseCentroidScreener.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_BaseCentroidChecker.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerBase.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSampler.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerWithResidueList.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerWithResidueAlternatives.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.screener.BaseCentroidScreener" );

namespace protocols {
namespace stepwise {
namespace screener {

	//Constructor
	BaseCentroidScreener::BaseCentroidScreener( modeler::rna::checker::RNA_BaseCentroidCheckerOP base_centroid_checker,
																							core::pose::PoseOP screening_pose,
																							bool const force_centroid_interaction /* = true */):
		base_centroid_checker_( base_centroid_checker ),
		screening_pose_( screening_pose ),
		force_centroid_interaction_( force_centroid_interaction ),
		using_stub_( false ),
		moving_res_base_stub_( core::kinematics::default_stub )
	{
	}

	//constructor
	BaseCentroidScreener::BaseCentroidScreener( modeler::rna::checker::RNA_BaseCentroidCheckerOP base_centroid_checker,
																							core::kinematics::Stub const & moving_res_base_stub ):
		base_centroid_checker_( base_centroid_checker ),
		screening_pose_( 0 ),
		force_centroid_interaction_( true ),
		using_stub_( true ),
		moving_res_base_stub_( moving_res_base_stub )
	{
	}

	//Destructor
	BaseCentroidScreener::~BaseCentroidScreener()
	{}


	////////////////////////////////////////////////////////////////////////////////
	bool
	BaseCentroidScreener::check_screen(){

		///////////////////////////////////
		//Reminder note of dependency: Important that base_stub_list is updated even in the case where
		// gap_size_ == 0 (and bulge is allowed)
		// since base_stub_list is used later below in the chain_break_screening section Jan 28, 2010 Parin S.

		if ( using_stub_ ){
			if ( !base_centroid_checker_->check_centroid_interaction( moving_res_base_stub_, count_data_ ) ) return false;
		} else {
			bool const found_a_centroid_interaction_partner = base_centroid_checker_->update_base_stub_list_and_check_centroid_interaction( *screening_pose_, count_data_ );
			//Essentially this doesn't screen for centroid interaction at chain_break.
			if ( force_centroid_interaction_ && !found_a_centroid_interaction_partner ) return false;
		}

		// Note that is does not update base_stub_list. To do that, use update_base_stub_list_and_check_that_terminal_res_are_unstacked
		if ( !base_centroid_checker_->check_that_terminal_res_are_unstacked() ) return false;

		return true;
	}

	////////////////////////////////////////////////////////////////////////////
	void
	BaseCentroidScreener::fast_forward( sampler::StepWiseSamplerBaseOP sampler ){
		using namespace sampler;
		using namespace sampler::rigid_body;

		if ( sampler->type() == RIGID_BODY_WITH_RESIDUE_LIST ){
			RigidBodyStepWiseSamplerWithResidueList & rigid_body_rotamer_with_copy_dofs = *( static_cast< RigidBodyStepWiseSamplerWithResidueList * >( sampler.get() ) );
			rigid_body_rotamer_with_copy_dofs.fast_forward_to_next_euler_gamma();
			return;
		} else if ( sampler->type() == RIGID_BODY_WITH_RESIDUE_ALTERNATIVES ){
			RigidBodyStepWiseSamplerWithResidueAlternatives & rigid_body_rotamer_with_residue_alternatives = *( static_cast< RigidBodyStepWiseSamplerWithResidueAlternatives * >( sampler.get() ) );
			rigid_body_rotamer_with_residue_alternatives.fast_forward_to_next_euler_gamma();
			return;
		}	 else if ( sampler->type() == RIGID_BODY ){
			RigidBodyStepWiseSampler & rigid_body_rotamer = *( static_cast< RigidBodyStepWiseSampler * >( sampler.get() ) );
			rigid_body_rotamer.fast_forward_to_next_euler_gamma();
		}

	}

} //screener
} //stepwise
} //protocols
