// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/working_parameters/StepWiseProteinWorkingParametersUtil.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampling/working_parameters/util.hh>
#include <protocols/stepwise/sampling/working_parameters/StepWiseProteinWorkingParameters.hh>
#include <protocols/stepwise/sampling/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelParameters.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <basic/Tracer.hh>
#include <utility/stream_util.hh>

static basic::Tracer TR( "protocols.stepwise.sampling.protein.working_parameters.util" );

using namespace core;

namespace protocols {
namespace stepwise {
namespace sampling {
namespace protein {
namespace working_parameters {

	StepWiseProteinWorkingParametersOP
	setup_protein_working_parameters_for_swa( utility::vector1< Size > const & moving_res_list,
																				pose::Pose const & pose,
																				pose::PoseCOP native_pose,
																				utility::vector1< Size > const & bridge_res,
																				utility::vector1< Size > const & working_minimize_res ){

		using namespace core::pose::full_model_info;

		FullModelInfo const & full_model_info = const_full_model_info( pose );
		std::string const & full_sequence = full_model_info.full_sequence();
		utility::vector1< Size > const & res_list = full_model_info.res_list();
		utility::vector1< Size > const & fixed_domain_map = full_model_info.fixed_domain_map();

		utility::vector1< Size > const & working_moving_res_list = moving_res_list;
		utility::vector1< Size > working_moving_suite_list;
		for ( Size n = 1; n <= moving_res_list.size(); n++ ){
			Size const & moving_res = moving_res_list[n];
			Size const parent_res = pose.fold_tree().get_parent_residue( moving_res );
			if ( pose.fold_tree().jump_nr( moving_res, parent_res ) == 0 ){
				working_moving_suite_list.push_back( std::min( moving_res, parent_res ) );
			}
		}

		utility::vector1< Size > is_working_res( full_sequence.size(), Size(0) );
		std::string working_sequence;
		for ( Size n = 1; n <= res_list.size(); n++ ) {
			is_working_res[ res_list[ n ] ] = 1;
			// following magic numbers match what is in StepWiseProteinPoseSetup:
			if ( moving_res_list.has_value( n ) ) is_working_res[ res_list[ n ] ] = 999 /*horrible HACK!*/;
			if ( bridge_res.has_value( res_list[ n ] ) ) is_working_res[ res_list[ n ] ] = 123 /*horrible HACK!*/;
			working_sequence += full_sequence[ n-1 ];
		}

		utility::vector1< Size > calc_rms_res_ = full_model_info.full_model_parameters()->get_res_list( CALC_RMS );
		if ( calc_rms_res_.size() == 0 ) calc_rms_res_ = full_model_info.sub_to_full( moving_res_list );

		utility::vector1< Size > fixed_res;
		utility::vector1< Size > const & extra_minimize_res = const_full_model_info( pose ).extra_minimize_res();
		for ( Size n = 1; n <= full_sequence.size(); n++ )	{
			if ( extra_minimize_res.has_value( n ) ) continue;
			if ( fixed_domain_map[ n ] > 0 ) {
				fixed_res.push_back( n );
				continue;
			}
			if ( working_minimize_res.size() > 0 ){ // finer specification  of what to minimize
				if( res_list.has_value( n ) && !working_minimize_res.has_value( res_list.index( n ) ) ) fixed_res.push_back( n );
			}
		}

		pose::PoseOP working_native_pose;
		if ( native_pose ) {
			working_native_pose = new pose::Pose;
			protocols::stepwise::sampling::pdbslice( *working_native_pose, *native_pose, res_list );
		}
		// ALSO WILL NEED A DEFAULT BRIDGE RES SETTER!
		StepWiseProteinWorkingParametersOP working_parameters = new StepWiseProteinWorkingParameters;
		working_parameters->set_full_to_sub( full_model_info.full_to_sub() ); // will update sub_to_full
		working_parameters->set_is_working_res( is_working_res );
		working_parameters->set_sequence( full_model_info.full_sequence() );
		//		working_parameters->set_working_sequence( working_sequence ); // auto-updated with set_sequence
		//		working_parameters->set_working_res_list( res_list ); // figured out from full_to_sub
		working_parameters->set_working_moving_res_list( working_moving_res_list );
		//working_parameters->set_working_moving_suite_list( working_moving_suite_list );  // set at the same time as working_moving_res_list
		//working_parameters->set_chain_boundaries(); // Not in use later. legacy of PoseSetup
		//working_parameters->set_which_chain_has_moving_res(); // Not in use later. legacy of PoseSetup
		//working_parameters->gap_size(); // Not in use later. legacy of PoseSetup
		//working_parameters->first_chain_break_res(); // Not in use later. legacy of PoseSetup
		//working_parameters->is_prepend(); // Not in use later. legacy of PoseSetup
		//working_parameters->is_internal(); // Not in use later. legacy of PoseSetup
		//working_parameters->partition_definition(); // Not in use later. legacy of PoseSetup
		working_parameters->set_working_native_pose( working_native_pose );

		working_parameters->set_working_fixed_res( full_model_info.full_to_sub( fixed_res ) );
		//working_parameters->set_working_terminal_res( full_model_info.full_to_sub( terminal_res ) ); // not in use for proteins.
		//		working_parameters->set_working_superimpose_res( full_model_info.full_to_sub( superimpose_res ) ); // not in use! assumes pose is already superimposed.
		working_parameters->set_working_calc_rms_res( full_model_info.full_to_sub( calc_rms_res_ ) );
		working_parameters->set_working_bridge_res( full_model_info.full_to_sub( bridge_res ) );
		//		working_parameters->set_working_moving_pos( full_model_info.full_to_sub( moving_pos ) ); // not in use! legacy of PoseSetup

		return working_parameters;

	}


} //working_parameters
} //protein
} //sampling
} //stepwise
} //protocols
