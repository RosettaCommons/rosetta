// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/ChainClosableGeometryStubBasedScreener.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/ChainClosableGeometryStubBasedScreener.hh>
#include <protocols/stepwise/sampling/rna/checker/ChainClosableGeometryChecker.hh>
#include <protocols/rotamer_sampler/copy_dofs/ResidueListRotamer.hh>
#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamerWithResidueList.hh>
#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamerWithResidueAlternatives.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.hh>
#include <utility/tools/make_vector1.hh>
#include <basic/Tracer.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.screener.ChainClosableGeometryStubBasedScreener" );

using namespace core;

namespace protocols {
namespace stepwise {
namespace screener {

	//constructor
	ChainClosableGeometryStubBasedScreener::ChainClosableGeometryStubBasedScreener( sampling::rna::checker::ChainClosableGeometryCheckerOP chain_closable_geometry_checker,
																																	utility::vector1< core::pose::PoseOP > screening_pose_list,
																																	core::kinematics::Stub const & moving_res_base_stub,
																																	Size const reference_res,
																																	utility::vector1< core::conformation::ResidueOP > moving_rsd_at_origin_list,
																																	bool const using_predefined_moving_rsd_list /* = true */ ):
		StepWiseResiduePairScreener( chain_closable_geometry_checker->five_prime_chain_break_res(),
																 chain_closable_geometry_checker->three_prime_chain_break_res() ),
		chain_closable_geometry_checker_( chain_closable_geometry_checker ),
		screening_pose_list_( screening_pose_list ),
		moving_rsd_at_origin_list_( moving_rsd_at_origin_list ),
		moving_res_base_stub_( moving_res_base_stub ),
		reference_res_( reference_res ),
		using_predefined_moving_rsd_list_( using_predefined_moving_rsd_list )
	{}

	//constructor
	ChainClosableGeometryStubBasedScreener::ChainClosableGeometryStubBasedScreener( sampling::rna::checker::ChainClosableGeometryCheckerOP chain_closable_geometry_checker,
																																	utility::vector1< core::pose::PoseOP > screening_pose_list,
																																	core::kinematics::Stub const & moving_res_base_stub,
																																	Size const reference_res ):
		StepWiseResiduePairScreener( chain_closable_geometry_checker->five_prime_chain_break_res(),
																 chain_closable_geometry_checker->three_prime_chain_break_res() ),
		chain_closable_geometry_checker_( chain_closable_geometry_checker ),
		screening_pose_list_( screening_pose_list ),
		moving_res_base_stub_( moving_res_base_stub ),
		reference_res_( reference_res ),
		using_predefined_moving_rsd_list_( false )
	{}

	//Destructor
	ChainClosableGeometryStubBasedScreener::~ChainClosableGeometryStubBasedScreener()
	{}

	///////////////////////////////////////////////////////////////////
	void
	ChainClosableGeometryStubBasedScreener::get_update( rotamer_sampler::RotamerBaseOP sampler ){

		using namespace rotamer_sampler;
		using namespace rotamer_sampler::rigid_body;

		if ( sampler->type() == RIGID_BODY_WITH_RESIDUE_LIST ){
			RigidBodyRotamerWithResidueList & rigid_body_rotamer_with_residue_list = *( static_cast< RigidBodyRotamerWithResidueList * >( sampler.get() ) );
			moving_rsd_at_origin_ = rigid_body_rotamer_with_residue_list.get_residue_at_origin();
			return;
		} else 		if ( sampler->type() == RIGID_BODY_WITH_RESIDUE_ALTERNATIVES ){
			RigidBodyRotamerWithResidueAlternatives & rigid_body_rotamer_with_residue_alternatives = *( static_cast< RigidBodyRotamerWithResidueAlternatives * >( sampler.get() ) );
			moving_rsd_at_origin_ = rigid_body_rotamer_with_residue_alternatives.get_residue_at_origin().clone();
			return;
		}
		runtime_assert ( sampler->type() == RESIDUE_LIST );
		ResidueListRotamer & copy_dofs_rotamer = *( static_cast< ResidueListRotamer * >( sampler.get() ) );
		moving_rsd_at_origin_ = copy_dofs_rotamer.get_residue_at_origin();
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	///////////////Chain_break_screening -- distance cut                     /////////////////
	//////////////////////////////////////////////////////////////////////////////////////////
	bool
	ChainClosableGeometryStubBasedScreener::check_screen() {
		if ( !using_predefined_moving_rsd_list_ ) moving_rsd_at_origin_list_ = utility::tools::make_vector1( moving_rsd_at_origin_ );
		return chain_closable_geometry_checker_->check_screen( screening_pose_list_, moving_rsd_at_origin_list_,
																									moving_res_base_stub_, reference_res_ );
	}

} //screener
} //stepwise
} //protocols
