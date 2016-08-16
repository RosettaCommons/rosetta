// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/screener/VDW_BinScreener.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/VDW_BinScreener.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_VDW_BinChecker.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerWithResidueList.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerWithResidueAlternatives.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/Stub.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.screener.VDW_BinScreener" );

namespace protocols {
namespace stepwise {
namespace screener {

//Constructor
VDW_BinScreener::VDW_BinScreener( modeler::rna::checker::RNA_VDW_BinCheckerOP vdw_bin_checker,
	core::pose::Pose & screening_pose,
	Size const moving_res ):
	vdw_bin_checker_( vdw_bin_checker),
	screening_pose_( screening_pose ),
	moving_res_( moving_res ),
	using_stub_( false ),
	screening_moving_rsd_at_origin_( /* 0 */ ),
	moving_res_base_stub_( core::kinematics::default_stub )
{}

VDW_BinScreener::VDW_BinScreener( modeler::rna::checker::RNA_VDW_BinCheckerOP vdw_bin_checker,
	pose::Pose & screening_pose,
	Size const moving_res,
	core::conformation::ResidueCOP screening_moving_rsd_at_origin,
	core::kinematics::Stub const & moving_res_base_stub ):
	vdw_bin_checker_( vdw_bin_checker),
	screening_pose_( screening_pose ),
	moving_res_( moving_res ),
	using_stub_( true ),
	screening_moving_rsd_at_origin_( screening_moving_rsd_at_origin ),
	moving_res_base_stub_( moving_res_base_stub )
{}

//Destructor
VDW_BinScreener::~VDW_BinScreener()
{}

bool
VDW_BinScreener::check_screen(){
	if ( using_stub_ ) {
		return vdw_bin_checker_->VDW_rep_screen( screening_pose_, moving_res_,
			*screening_moving_rsd_at_origin_, moving_res_base_stub_ );
	}
	return vdw_bin_checker_->VDW_rep_screen( screening_pose_, moving_res_ );
}

////////////////////////////////////////////////////////////////////////////
void
VDW_BinScreener::fast_forward( sampler::StepWiseSamplerBaseOP sampler ){
	using namespace sampler;
	using namespace sampler::rigid_body;
	if ( using_stub_ ) {
		if ( sampler->type() == RIGID_BODY_WITH_RESIDUE_LIST ) {
			RigidBodyStepWiseSamplerWithResidueList & rigid_body_rotamer_with_copy_dofs = *( static_cast< RigidBodyStepWiseSamplerWithResidueList * >( sampler.get() ) );
			rigid_body_rotamer_with_copy_dofs.fast_forward_to_next_rigid_body();
		} else if ( sampler->type() == RIGID_BODY_WITH_RESIDUE_ALTERNATIVES ) {
			RigidBodyStepWiseSamplerWithResidueAlternatives & rigid_body_rotamer_with_residue_alternatives = *( static_cast< RigidBodyStepWiseSamplerWithResidueAlternatives * >( sampler.get() ) );
			rigid_body_rotamer_with_residue_alternatives.fast_forward_to_next_rigid_body();
		}
	}
}

} //screener
} //stepwise
} //protocols
