// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/working_parameters/StepWiseRNA_WorkingParametersUtil.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/working_parameters/util.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/modeler/protein/util.hh>
#include <protocols/stepwise/legacy/modeler/rna/StepWiseRNA_WorkingParametersSetup.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/movemap/util.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/toolbox/AtomLevelDomainMap.hh>
#include <protocols/stepwise/sampler/rna/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelParameters.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/util.hh>
#include <utility/tools/make_vector1.hh>
#include <ObjexxFCL/string.functions.hh>

#include <utility/string_util.hh>
#include <utility/stream_util.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.modeler.working_parameters.util" );

using namespace core;
using utility::tools::make_vector1;
using namespace protocols::stepwise::modeler::rna;
using namespace protocols::stepwise::modeler::working_parameters;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace working_parameters {

///////////////////////////////////////////////////////////////////////////////////////////////////
// Originally developed within StepWiseRNA_Modeler to allow wrapping around
//  Parin's SWA functions. Now moved out to enable simplification and/or refactoring.
//
// Briefly, we need to make a StepWiseRNA_WorkingParameters object that will
// be fed into various StepWiseRNA movers.
// This carried a grab bag of residue lists referring to the global pose, the working pose,
// sequence mappings, "is_prepend_map", etc.
//    -- rhiju
//
/////////////////////////////////////////////////////////////////////////
StepWiseWorkingParametersOP
setup_working_parameters_for_swa( utility::vector1< Size > const & moving_res_list,
	pose::Pose const & pose,
	pose::PoseCOP native_pose,
	utility::vector1< Size > const & user_defined_bridge_res,
	utility::vector1< Size > const & working_minimize_res ){

	using namespace core::pose::full_model_info;

	StepWiseWorkingParametersOP working_parameters( new working_parameters::StepWiseWorkingParameters );
	if ( moving_res_list.size() == 1 ) working_parameters = setup_working_parameters_explicit( moving_res_list[1], pose, native_pose );

	FullModelInfo const & full_model_info = const_full_model_info( pose );
	std::string const & full_sequence = full_model_info.full_sequence();
	utility::vector1< Size > const & res_list = full_model_info.res_list();
	utility::vector1< Size > const & fixed_domain_map = full_model_info.fixed_domain_map();

	utility::vector1< Size > const & working_moving_res_list = moving_res_list;
	utility::vector1< Size > working_moving_suite_list;
	bool is_prepend( false );
	std::map< Size, bool > is_prepend_map;
	for ( Size n = moving_res_list.size(); n >= 1; n-- ) { // backwards to set is_prepend to value at initial moving_res.
		Size const & moving_res = moving_res_list[n];
		Size const parent_res = pose.fold_tree().get_parent_residue( moving_res );
		if ( pose.fold_tree().jump_nr( moving_res, parent_res ) == 0 ) working_moving_suite_list.push_back( std::min( moving_res, parent_res ) );
		is_prepend = ( moving_res < parent_res );
		is_prepend_map[ res_list[ moving_res ] ] = is_prepend;
	}

	utility::vector1< Size > bridge_res = user_defined_bridge_res; // for protein loop closure.
	if ( bridge_res.size() == 0 ) bridge_res = protein::get_bridge_res( pose, moving_res_list );

	utility::vector1< Size > is_working_res( full_sequence.size(), Size(0) );
	//std::string working_sequence;
	for ( Size n = 1; n <= res_list.size(); n++ ) {
		is_working_res[ res_list[ n ] ] = 1;
		// following magic numbers match what is in StepWiseProteinPoseSetup:
		if ( moving_res_list.has_value( n ) )        is_working_res[ res_list[ n ] ] = MOVING_RES;
		if ( bridge_res.has_value( res_list[ n ] ) ) is_working_res[ res_list[ n ] ] = BRIDGE_RES;
		//working_sequence += full_sequence[ n-1 ];
	}

	utility::vector1< Size > calc_rms_res_ = full_model_info.full_model_parameters()->get_res_list( CALC_RMS );
	if ( calc_rms_res_.size() == 0 ) calc_rms_res_ = full_model_info.sub_to_full( moving_res_list );

	utility::vector1< Size > fixed_res;
	utility::vector1< Size > const & extra_minimize_res = const_full_model_info( pose ).extra_minimize_res();
	utility::vector1< Size > const & sample_res = const_full_model_info( pose ).sample_res();
	for ( Size n = 1; n <= full_sequence.size(); n++ ) {
		if ( extra_minimize_res.has_value( n ) ) continue;
		if ( sample_res.has_value( n ) ) continue;
		if ( fixed_domain_map[ n ] > 0  && !working_minimize_res.has_value( res_list.index( n ) ) ) {
			fixed_res.push_back( n );
			continue;
		}
	}

	utility::vector1< Size > const working_moving_partition_res = figure_out_moving_partition_res( pose, moving_res_list );
	ObjexxFCL::FArray1D< bool > partition_definition( pose.total_residue(), false );
	for ( Size n = 1; n <= pose.total_residue(); n++ ) partition_definition( n ) = working_moving_partition_res.has_value( n );

	pose::PoseOP working_native_pose;
	if ( native_pose != 0 ) {
		working_native_pose = native_pose->clone();
		utility::vector1< Size > const & native_res_list = get_res_list_from_full_model_info( *working_native_pose );
		utility::vector1< Size > res_list_for_slicing;
		for ( Size n = 1; n <= res_list.size(); n++ ) {
			if ( native_res_list.has_value( res_list[n] ) ) res_list_for_slicing.push_back( native_res_list.index( res_list[n] ) );
		}
		core::pose::pdbslice( *working_native_pose, *native_pose, res_list_for_slicing );
	}


	working_parameters->set_full_to_sub( full_model_info.full_to_sub() ); // will update sub_to_full
	working_parameters->set_is_working_res( is_working_res );
	working_parameters->set_sequence( full_model_info.full_sequence() );
	working_parameters->set_is_prepend( is_prepend );
	working_parameters->set_is_prepend_map( is_prepend_map );
	//  working_parameters->set_working_sequence( working_sequence ); // auto-updated with set_sequence
	//  working_parameters->set_working_res_list( res_list ); // figured out from full_to_sub
	working_parameters->set_working_moving_res_list( working_moving_res_list );
	//working_parameters->set_working_moving_suite_list( working_moving_suite_list );  // set at the same time as working_moving_res_list
	//working_parameters->set_chain_boundaries(); // Not in use later. legacy of PoseSetup
	//working_parameters->set_which_chain_has_moving_res(); // Not in use later. legacy of PoseSetup
	//working_parameters->gap_size(); // Not in use later. legacy of PoseSetup
	//working_parameters->first_chain_break_res(); // Not in use later. legacy of PoseSetup
	//working_parameters->is_prepend(); // Not in use later. legacy of PoseSetup
	//working_parameters->is_internal(); // Not in use later. legacy of PoseSetup
	working_parameters->set_working_native_pose( working_native_pose );

	working_parameters->set_working_fixed_res( working_parameters->apply_full_to_sub_mapping( fixed_res ) );
	working_parameters->set_working_best_alignment( full_model_info.full_to_sub( fixed_res ) );
	working_parameters->set_working_native_alignment( full_model_info.full_to_sub( fixed_res ) );
	working_parameters->set_native_alignment( fixed_res );

	//  working_parameters->set_working_superimpose_res( full_model_info.full_to_sub( superimpose_res ) ); // not in use! assumes pose is already superimposed.
	working_parameters->set_working_calc_rms_res( full_model_info.full_to_sub( calc_rms_res_ ) );
	working_parameters->set_working_bridge_res( full_model_info.full_to_sub( bridge_res ) );
	working_parameters->set_working_moving_partition_res( working_moving_partition_res );
	working_parameters->set_partition_definition( partition_definition );

	working_parameters->set_fold_tree( pose.fold_tree() );

	return working_parameters;
}


/////////////////////////////////////////////////////////////////////////////////////////////////
StepWiseWorkingParametersOP
setup_working_parameters_explicit( Size const rebuild_res, /* this must be in moving partition */
	pose::Pose const & pose,
	pose::PoseCOP native_pose ){
	using namespace chemical;
	using namespace pose::full_model_info;

	FullModelInfo const & full_model_info = const_full_model_info( pose );
	utility::vector1< Size > const & calc_rms_res = full_model_info.calc_rms_res();
	utility::vector1< Size > const & terminal_res = full_model_info.rna_terminal_res();
	utility::vector1< Size > const & block_stack_above_res = full_model_info.rna_block_stack_above_res();
	utility::vector1< Size > const & block_stack_below_res = full_model_info.rna_block_stack_below_res();
	utility::vector1< Size > const & syn_chi_res_list = full_model_info.rna_syn_chi_res();
	utility::vector1< Size > const & anti_chi_res_list = full_model_info.rna_anti_chi_res();
	utility::vector1< Size > const & north_sugar_list = full_model_info.full_model_parameters()->get_res_list( RNA_NORTH_SUGAR );
	utility::vector1< Size > const & south_sugar_list = full_model_info.full_model_parameters()->get_res_list( RNA_SOUTH_SUGAR );

	Size moving_res( rebuild_res ), rebuild_suite( 0 );
	bool floating_base( false ), is_internal( false );
	Size const reference_res = pose.fold_tree().get_parent_residue( rebuild_res, floating_base /*is_jump*/ );
	// for internal moves, need to be smart about input_res definitions -- 'domains' that are separated by moving residue.
	// this loop will also determine any chainbreak that requires closure
	utility::vector1< bool > partition_definition;
	Size floating_base_anchor_res( 0 );
	if ( floating_base ) {
		floating_base_anchor_res = reference_res;
		partition_definition = get_partition_definition_floating_base( pose, rebuild_res );
	} else {
		runtime_assert( (rebuild_res == reference_res + 1) || (rebuild_res == reference_res - 1 ) );
		rebuild_suite = ( rebuild_res < reference_res ) ? rebuild_res : reference_res;
		partition_definition = get_partition_definition( pose, rebuild_suite );
	}

	utility::vector1< Size > input_res1, input_res2 /*blank*/, cutpoint_open, cutpoint_closed;
	//Size cutpoint_closed_distal( 0 );
	for ( Size n = 1; n <= pose.total_residue(); n++ ) {
		if ( !partition_definition[ n ] ) input_res1.push_back( n );
		else input_res2.push_back( n );
	}
	utility::vector1< Size > const & moving_partition_res = ( partition_definition[ moving_res ] ) ? input_res2 :  input_res1;
	utility::vector1< Size > five_prime_chain_breaks, three_prime_chain_breaks, chain_break_gap_sizes;
	figure_out_moving_chain_breaks( pose, moving_partition_res,
		cutpoint_closed,
		five_prime_chain_breaks, three_prime_chain_breaks, chain_break_gap_sizes );
	for ( Size n = 1; n <= pose.total_residue(); n++ ) {
		if ( pose.fold_tree().is_cutpoint( n ) && !five_prime_chain_breaks.has_value( n ) ) {
			cutpoint_open.push_back( n );
		}
	}
	Size five_prime_chain_break_res( 0 ), gap_size( 999 );
	for ( Size i = 1; i <= five_prime_chain_breaks.size(); i++ ) {
		Size const n = five_prime_chain_breaks[ i ];
		if ( ( reference_res < rebuild_res && rebuild_res == (n+1) ) ||
				( reference_res > rebuild_res && rebuild_res == n ) ) {
			runtime_assert( floating_base ); continue;
		}
		five_prime_chain_break_res = n;
		gap_size = chain_break_gap_sizes[ i ];
		//if ( cutpoint_closed.has_value( n ) ) cutpoint_closed_distal = n; // same as five_prime_chain_break_res
		break;
	}

	// TR << TR.Blue << "INPUT_RES1 " << make_tag_with_dashes(input_res1) << " INPUT_RES2 " << make_tag_with_dashes(input_res2) << " REBUILD_RES " << rebuild_res << " REBUILD_SUITE " <<  rebuild_suite << " CUTPOINT_CLOSED " << cutpoint_closed << " CUTPOINT_CLOSED_DISTAL " << cutpoint_closed_distal << " FLOATING_BASE " << floating_base << TR.Reset << std::endl;
	// TR << pose.fold_tree() << std::endl;

	// if there's really just a single nucleotide being rebuilt, its nucleoside is not fixed.
	// but if there's a whole chunk of stuff, its sugar & base are assumed fixed.
	//note that working_fixed_res_guess, which is really a list of fixed nucleosides,
	// should now include the 'moving res' [unless re-specified by user down below].

	// check if this is an 'internal' move.
	utility::vector1< Size> working_fixed_res_guess; // this is just a placeholder now, I think.
	if ( moving_partition_res.size() > 1 ) {
		if ( !working_fixed_res_guess.has_value( rebuild_res ) ) working_fixed_res_guess.push_back( rebuild_res );
		if ( !floating_base ) {
			is_internal = true;
			//moving_res = rebuild_suite; // parin-style old convention, where moving_res = moving_suite for internal suites. Sigh.
		}
	}

	std::string full_sequence = full_model_info.full_sequence();
	Size nres = full_sequence.size();
	if ( full_sequence[ nres - 1] == 'X' ) {
		full_sequence = full_sequence.substr( 0, nres - 1 );
		nres -= 1;
	}
	utility::vector1< bool > is_working_res;
	for ( Size n = 1; n <= nres; n++ ) is_working_res.push_back( full_model_info.res_list().has_value( n ) );
	bool const is_prepend = ( rebuild_res < reference_res );

	std::map< Size, bool > is_prepend_map; // used for rmsd calculations -- include phosphate. May not be necessary.
	utility::vector1< Size > calc_rms_res_ = calc_rms_res;
	if ( calc_rms_res_.size() == 0 ) calc_rms_res_.push_back( full_model_info.sub_to_full( rebuild_res ) );
	for ( Size i = 1; i <= calc_rms_res_.size(); i++ ) {
		Size const & n = calc_rms_res_[ i ];
		if ( is_working_res[ n ] ) {
			is_prepend_map[ n ] = is_prepend; //pose.residue_type( full_model_info.full_to_sub( n ) ).has_variant_type( "VIRTUAL_PHOSPHATE" );
		}
	}

	StepWiseWorkingParametersOP working_parameters( new working_parameters::StepWiseWorkingParameters );
	working_parameters->set_output_extra_RMSDs( false );
	working_parameters->set_is_simple_full_length_job_params( false );
	working_parameters->set_full_to_sub( full_model_info.full_to_sub() ); // will update sub_to_full
	working_parameters->set_full_sequence( full_sequence ); // will update working_sequence.
	working_parameters->set_moving_res( full_model_info.sub_to_full( moving_res ) );
	working_parameters->set_is_prepend( is_prepend );
	working_parameters->set_is_prepend_map( is_prepend_map );
	working_parameters->set_is_internal( is_internal );
	working_parameters->set_floating_base( floating_base );
	working_parameters->set_floating_base_anchor_res( full_model_info.sub_to_full( floating_base_anchor_res ) );
	working_parameters->set_rebuild_bulge_mode( figure_out_rebuild_bulge_mode( pose, rebuild_res ) );
	working_parameters->set_sample_both_sugar_base_rotamer( figure_out_sample_both_sugar_base_rotamer( pose, floating_base, rebuild_suite ) );
	working_parameters->set_partition_definition( partition_definition );
	working_parameters->set_is_working_res( is_working_res ); // may not be in use later. anyway, set it, just in case.
	working_parameters->set_working_moving_res_list( make_vector1( moving_res ) ); // also sets up working_moving_suite and working_moving_suite_list
	//working_parameters->set_chain_boundaries(); // Not in use later.
	working_parameters->set_five_prime_chain_break_res( five_prime_chain_break_res );
	working_parameters->set_gap_size( gap_size );
	working_parameters->set_working_native_pose( native_pose );
	working_parameters->set_working_fixed_res( working_fixed_res_guess ); // will get updated later, I think.
	working_parameters->set_calc_rms_res( calc_rms_res_ );
	working_parameters->set_terminal_res( terminal_res );
	working_parameters->set_block_stack_above_res( block_stack_above_res );
	working_parameters->set_block_stack_below_res( block_stack_below_res );
	working_parameters->set_working_moving_partition_res( moving_partition_res );
	working_parameters->set_input_res_vectors( make_vector1( full_model_info.sub_to_full( input_res1 ),
		full_model_info.sub_to_full( input_res2 ) ) );
	working_parameters->set_cutpoint_closed_list( full_model_info.sub_to_full( cutpoint_closed ) );
	working_parameters->set_cutpoint_open_list( full_model_info.sub_to_full( cutpoint_open ) );
	working_parameters->set_working_best_alignment( working_fixed_res_guess );
	working_parameters->set_working_native_alignment( working_fixed_res_guess );
	working_parameters->set_native_alignment( full_model_info.sub_to_full( working_fixed_res_guess ) ); // is this in use?
	working_parameters->set_global_sample_res_list( full_model_info.sub_to_full( make_vector1( moving_res ) ) );
	working_parameters->set_force_syn_chi_res_list( syn_chi_res_list ); // will get automatically converted to working numbering.
	working_parameters->set_force_anti_chi_res_list( anti_chi_res_list ); // will get automatically converted to working numbering.
	working_parameters->set_force_north_sugar_list( north_sugar_list ); // will get automatically converted to working numbering.
	working_parameters->set_force_south_sugar_list( south_sugar_list ); // will get automatically converted to working numbering.
	working_parameters->set_fold_tree( pose.fold_tree() );

	TR.Debug << "past working_parameters initialization " << std::endl;

	return working_parameters;
}

/////////////////////////////////////////////////////////////////////////
bool
figure_out_rebuild_bulge_mode( pose::Pose const & pose, Size const rebuild_res ){
	kinematics::FoldTree const & f = pose.fold_tree();

	// if the rebuild_res is connected to anything by a jump, better not virtualize it.
	for ( Size n = 1; n <= f.num_jump(); n++ ) {
		if ( f.upstream_jump_residue( n )   == static_cast<int>(rebuild_res) ) return false;
		if ( f.downstream_jump_residue( n ) == static_cast<int>(rebuild_res) ) return false;
	}

	if ( rebuild_res > 1 &&
			pose.residue( rebuild_res - 1 ).has_variant_type( core::chemical::VIRTUAL_RIBOSE ) &&
			( !f.is_cutpoint( rebuild_res - 1 ) ||
			is_cutpoint_closed( pose, rebuild_res - 1 ) ) &&
			!f.is_cutpoint( rebuild_res ) &&
			f.jump_nr( rebuild_res - 1, rebuild_res + 1) > 0 ) return true;
	if ( rebuild_res < pose.total_residue() &&
			pose.residue( rebuild_res + 1 ).has_variant_type( core::chemical::VIRTUAL_RIBOSE ) &&
			( !f.is_cutpoint( rebuild_res ) ||
			is_cutpoint_closed( pose, rebuild_res ) ) &&
			!f.is_cutpoint( rebuild_res - 1 ) &&
			f.jump_nr( rebuild_res - 1, rebuild_res + 1 ) > 0 ) return true;
	return false;
}

/////////////////////////////////////////////////////////////////////////
bool
figure_out_sample_both_sugar_base_rotamer( pose::Pose const & pose,
	bool const floating_base,
	Size const rebuild_suite ){
	if ( !floating_base &&
			sampler::rna::modeler_sugar_at_five_prime(  pose, rebuild_suite ) &&
			sampler::rna::modeler_sugar_at_three_prime( pose, rebuild_suite ) ) {
		if ( !pose.residue( rebuild_suite ).has_variant_type( core::chemical::VIRTUAL_RIBOSE ) &&
				!pose.residue( rebuild_suite + 1 ).has_variant_type( core::chemical::VIRTUAL_RIBOSE ) ) {
			return true;
		}
	}

	return false;
}

} //working_parameters
} //modeler
} //stepwise
} //protocols
