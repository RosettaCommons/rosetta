// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/rna/StepWiseRNA_JobParametersUtil.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampling/rna/StepWiseRNA_JobParametersUtil.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/stepwise/sampling/rna/legacy/StepWiseRNA_JobParametersSetup.hh>
#include <protocols/stepwise/StepWiseUtil.hh>
#include <protocols/toolbox/AllowInsert.hh>
#include <protocols/rotamer_sampler/rna/RNA_RotamerSamplerUtil.hh>
#include <core/chemical/VariantType.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelInfoUtil.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/string_util.hh>
#include <ObjexxFCL/string.functions.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.sampling.rna.StepWiseRNA_JobParametersUtil" );

using namespace core;
using utility::tools::make_vector1;

namespace protocols {
namespace stepwise {
namespace sampling {
namespace rna {

///////////////////////////////////////////////////////////////////////////////////////////////////
// Originally developed within StepWiseRNA_Modeler to allow wrapping around
//  Parin's SWA functions. Now moved out to enable simplification and/or refactoring.
//
// Briefly, we need to make a StepWiseRNA_JobParameters object that will
// be fed into various StepWiseRNA movers.
// This carried a grab bag of residue lists referring to the global pose, the working pose,
// sequence mappings, "is_prepend_map", etc.
//    -- rhiju
//
StepWiseRNA_JobParametersOP
setup_job_parameters_for_swa( Size const rebuild_res,
															pose::Pose const & pose,
															pose::PoseCOP native_pose,
															utility::vector1< Size > const & rmsd_res_list,
															utility::vector1< Size > const & terminal_res,
															utility::vector1< Size > const & syn_chi_res_list,
															utility::vector1< Size > const & extra_minimize_res,
															utility::vector1< Size > & fixed_res,
															utility::vector1< Size > & minimize_res,
															toolbox::AllowInsertOP & allow_insert,
															kinematics::MoveMapOP & minimize_move_map ) {
	using namespace pose;
	using namespace id;
	using namespace protocols::stepwise;
	using namespace pose::full_model_info;

	// FullModelInfo can be setup with the command make_sure_full_model_info_is_setup( pose ) before calling this function.
	// It is used to figure out which suites need to be minimized (using info in domain map), and
	// is also necessary if there are jumps to floating bases -- need to know how many intervening residues are skipped.
	FullModelInfo const & full_model_info = const_full_model_info( pose );

	StepWiseRNA_JobParametersOP job_parameters = new StepWiseRNA_JobParameters;
	utility::vector1< Size > suites_that_must_be_minimized;
	std::string full_sequence = pose.sequence();

	// what if there is a virtual residue? need to remove it, actually, before running stepwise_rna_job_parameters_setup.
	Size nres = pose.total_residue();
	utility::vector1< Size > not_rebuild_res;
	for ( Size n = 1; n <= nres; n++ ) if ( n != rebuild_res ) not_rebuild_res.push_back( n );
	utility::vector1< Size > fixed_res_guess = not_rebuild_res; // may be revised below.

	if ( rebuild_res > 0 ){
		// deprecate legacy at end of 2014 if not needed
		//StepWiseRNA_JobParametersOP job_parameters_legacy = setup_job_parameters_legacy( rebuild_res, pose, native_pose, rmsd_res_list,
		//		 																																								 terminal_res, syn_chi_res_list, fixed_res_guess,
		//		 																																								 suites_that_must_be_minimized );
		//		print_JobParameters_info( job_parameters_legacy, "LEGACY", TR );

		StepWiseRNA_JobParametersOP job_parameters_explicit = setup_job_parameters_explicit( rebuild_res, pose, native_pose, rmsd_res_list,
																																												 terminal_res, syn_chi_res_list, fixed_res_guess,
																																												 suites_that_must_be_minimized );
		//		print_JobParameters_info( job_parameters_explicit, "EXPLICIT", TR );
		job_parameters = job_parameters_explicit;
	}


	// a little touch up. there is a re-rooting of the pose that happens later, but the pose is const here.
	// however, this is our opportunity to update job_parameters.
	//	utility::vector1< Size > root_partition_res, moving_partition_res;
	//figure_out_root_partition_res( pose, job_parameters, root_partition_res, moving_partition_res );
	//	job_parameters->set_working_moving_partition_pos( moving_partition_res );

	// If setup_job_parameters_for_stepwise_) is called, then the user has not supplied their own StepWiseRNA_JobParameters object to the modeler. This means that StepWiseRNAMinimizer will not be generating a move map on its own, so we need to supply it with an AllowInsert object in order to handle the possibility of variable bond geometries.
	allow_insert = new toolbox::AllowInsert(pose); // Default constructor that allows everything to move

	// user input minimize_res...
	if ( minimize_res.size() > 0 ) { // specifying more residues which could move during the minimize step -- standard for SWM.
		fixed_res.clear();
		for ( Size n = 1; n <= nres; n++ ) {
			if ( !minimize_res.has_value( n ) )	{
				fixed_res.push_back( n );
				allow_insert->set( n, false );
			}
		}
	} else if ( fixed_res.size() > 0 ){ // how 'standard' SWA specifies moving information.
		runtime_assert( minimize_res.size() == 0 );
		for ( Size n = 1; n <= nres; n++ ) {
			if ( !fixed_res.has_value( n ) )	{
				minimize_res.push_back( n );
			} else {
				allow_insert->set( n, false );
			}
		}
	} else { // 'reasonable' default behavior, inferred above.
		for ( Size n = 1; n <= nres; n++ ) {
			if ( !fixed_res_guess.has_value( n ) )	{
				minimize_res.push_back( n );
			} else {
				allow_insert->set( n, false );
			}
		}
		fixed_res = fixed_res_guess;
	}

	if ( fixed_res.size() > 0 ) job_parameters->set_working_fixed_res( fixed_res ); // is this necessary if we just supply movemap?
	job_parameters->set_working_native_alignment( fixed_res );

	//Now we perform the additional task of updating the AllowInsert object based on any optional additional residues that the user wants minimized, as specified in minimizer_extra_minimize_res. The intended mode of operation is that the user decides on extra_minimize_res at the high level for an entire SWM run, and then changes minimize_res to control a specific minimization event.
	update_allow_insert_with_extra_minimize_res( pose, allow_insert, extra_minimize_res );

	minimize_move_map = new kinematics::MoveMap;
	//figure_out_stepwise_rna_movemap( *minimize_move_map, pose, minimize_res );
	figure_out_stepwise_rna_movemap( *minimize_move_map, pose, allow_insert );

	// last, but not least, there might be some information in the domain map. Note
	// that generally we could instead replace fixed_res with an inputted domain map.
	// that is, get rid of fixed_res & minimize_res and instead have a local fixed_domain_map,
	// which can instead be updated by set_fixed_res.
	// how to tell Modeler to *not* minimize additional suites?
	utility::vector1< Size > const & res_list = full_model_info.res_list();
	utility::vector1< Size > const & fixed_domain_map = full_model_info.fixed_domain_map();
	utility::vector1< Size > const & cutpoint_open_in_full_model = full_model_info.cutpoint_open_in_full_model();
	for ( Size n = 1; n < pose.total_residue(); n++ ){
		if ( !cutpoint_open_in_full_model.has_value( res_list[ n ] ) &&
				( res_list[ n + 1 ]  == res_list[ n ] + 1 ) &&
				 !minimize_res.has_value( n )  && !minimize_res.has_value( n+1 )  &&
				 ( fixed_domain_map[ res_list[ n + 1 ] ] !=  fixed_domain_map[ res_list[ n ] ] ) &&
				 !suites_that_must_be_minimized.has_value( n ) ){
			TR.Debug << "ADDING NEW SUITE TO BE MINIMIZED BASED ON LOCATION AT DOMAIN BOUNDARY: " << n << std::endl;
			suites_that_must_be_minimized.push_back( n );
		}
	}
	TR.Debug << "SUITES_THAT_MUST_BE_MINIMIZED " << suites_that_must_be_minimized << std::endl;
	for ( Size n = 1; n <= suites_that_must_be_minimized.size(); n++ ){
		Size const suite_num = suites_that_must_be_minimized[ n ];
		minimize_move_map->set( TorsionID( suite_num,   id::BB, 5 ), true ); // epsilon
		minimize_move_map->set( TorsionID( suite_num,   id::BB, 6 ), true ); // zeta
		minimize_move_map->set( TorsionID( suite_num+1, id::BB, 1 ), true ); // alpha
		minimize_move_map->set( TorsionID( suite_num+1, id::BB, 2 ), true ); // beta
		minimize_move_map->set( TorsionID( suite_num+1, id::BB, 3 ), true ); // gamma
	}
	if ( job_parameters->floating_base() ){
		TR.Debug << TR.Blue << "WILL MINIMIZE JUMP: " << pose.fold_tree().get_jump_that_builds_residue( rebuild_res ) << TR.Reset << std::endl;
		minimize_move_map->set_jump( pose.fold_tree().get_jump_that_builds_residue( rebuild_res ), true );
	}

	return job_parameters;
}


///////////////////////////////////////////////////////////////////////////////
void
figure_out_root_partition_res( pose::Pose const & pose,
															 StepWiseRNA_JobParametersCOP job_parameters,
															 utility::vector1< Size > & root_partition_res,
															 utility::vector1< Size > & moving_partition_res ) {
	root_partition_res.clear();
	moving_partition_res.clear();
	Size const & moving_res = job_parameters->working_moving_res();
	if ( moving_res > 0 ){
		ObjexxFCL::FArray1D < bool > const & partition_definition = job_parameters->partition_definition();
		utility::vector1< Size > partition_res0, partition_res1;
		for ( Size i = 1; i <= pose.total_residue(); i++ ){
			if ( !partition_definition( i ) ) partition_res0.push_back( i );
			else partition_res1.push_back( i );
		}
		int root_partition( -1 );
		if ( partition_res0.size() > partition_res1.size() ){
			root_partition = 0;
			runtime_assert( partition_res1.has_value( moving_res ) );
		} else if ( partition_res0.size() < partition_res1.size() ){
			root_partition = 1;
			runtime_assert( partition_res0.has_value( moving_res ) );
		} else { // internal. moving_res corresponds to 5' residue of suite that is moved.
			root_partition = ( partition_res1.has_value( moving_res ) );
		}

		if ( root_partition == 0 ){
			for ( Size i = 1; i <= partition_res0.size(); i++ ) root_partition_res.push_back(   partition_res0[i] );
			for ( Size i = 1; i <= partition_res1.size(); i++ ) moving_partition_res.push_back( partition_res1[i] );
		} else {
			for ( Size i = 1; i <= partition_res0.size(); i++ ) moving_partition_res.push_back( partition_res0[i] );
			for ( Size i = 1; i <= partition_res1.size(); i++ ) root_partition_res.push_back(   partition_res1[i] );
		}
	} else {
		for ( Size i = 1; i <= pose.total_residue(); i++ ) root_partition_res.push_back( i );
	}

}


/////////////////////////////////////////////////////////////////////////////////////////////////
StepWiseRNA_JobParametersOP
setup_job_parameters_explicit( Size const rebuild_res, /* this must be in moving partition */
															 pose::Pose const & pose,
															 pose::PoseCOP native_pose,
															 utility::vector1< Size > const & rmsd_res_list,
															 utility::vector1< Size > const & terminal_res,
															 utility::vector1< Size > const & syn_chi_res_list,
															 utility::vector1< Size > & fixed_res_guess,
															 utility::vector1< Size > & suites_that_must_be_minimized ){
	using namespace chemical;
	using namespace pose::full_model_info;

	FullModelInfo const & full_model_info = const_full_model_info( pose );

	TR.Debug << pose.fold_tree() << std::endl;
	TR << "Rebuild residue: " << rebuild_res << std::endl;

	Size moving_res( rebuild_res ), rebuild_suite( 0 );
	bool floating_base( false ), is_internal( false );
	Size const reference_res = pose.fold_tree().get_parent_residue( rebuild_res, floating_base /*is_jump*/ );
	// for internal moves, need to be smart about input_res definitions -- 'domains' that are separated by moving residue.
	// this loop will also determine any chainbreak that requires closure
	utility::vector1< bool > partition_definition;
	Size floating_base_anchor_res( 0 );
	if ( floating_base ){
		floating_base_anchor_res = reference_res;
		partition_definition = get_partition_definition_floating_base( pose, rebuild_res );
	} else {
		runtime_assert( (rebuild_res == reference_res + 1) || (rebuild_res == reference_res - 1 ) );
		rebuild_suite = ( rebuild_res < reference_res ) ? rebuild_res : reference_res;
		partition_definition = get_partition_definition( pose, rebuild_suite );
	}

	utility::vector1< Size > input_res1, input_res2 /*blank*/, cutpoint_open, cutpoint_closed;
	Size cutpoint_closed_distal( 0 );
	for ( Size n = 1; n <= pose.total_residue(); n++ ){
		if ( !partition_definition[ n ] ) input_res1.push_back( n );
		else input_res2.push_back( n );
	}
	utility::vector1< Size > const & moving_partition_res = ( partition_definition[ moving_res ] ) ? input_res2 :  input_res1;
	utility::vector1< Size > five_prime_chain_breaks, three_prime_chain_breaks, chain_break_gap_sizes;
	figure_out_moving_chain_breaks( pose, moving_partition_res,
																	cutpoint_closed,
																	five_prime_chain_breaks, three_prime_chain_breaks, chain_break_gap_sizes );
	for ( Size n = 1; n <= pose.total_residue(); n++ ){
		if ( pose.fold_tree().is_cutpoint( n ) && !five_prime_chain_breaks.has_value( n ) ){
			cutpoint_open.push_back( n );
		}
	}
	Size five_prime_chain_break_res( 0 ), gap_size( 999 );
	for ( Size i = 1; i <= five_prime_chain_breaks.size(); i++ ){
		Size const n = five_prime_chain_breaks[ i ];
		if ( ( reference_res < rebuild_res && rebuild_res == (n+1) ) ||
				 ( reference_res > rebuild_res && rebuild_res == n )	){
			runtime_assert( floating_base ); continue;
		}
		five_prime_chain_break_res = n;
		gap_size = chain_break_gap_sizes[ i ];
		if ( cutpoint_closed.has_value( n ) ) cutpoint_closed_distal = n; // same as five_prime_chain_break_res
		break;
	}

	TR << "INPUT_RES1 " << make_tag_with_dashes(input_res1) << " INPUT_RES2 " << make_tag_with_dashes(input_res2) << " REBUILD_RES " << rebuild_res << " REBUILD_SUITE " <<  rebuild_suite << " CUTPOINT_CLOSED " << cutpoint_closed << " CUTPOINT_CLOSED_DISTAL " << cutpoint_closed_distal << " FLOATING_BASE " << floating_base << std::endl;
	TR << pose.fold_tree() << std::endl;

	// if there's really just a single nucleotide being rebuilt, its nucleoside is not fixed.
	// but if there's a whole chunk of stuff, its sugar & base are assumed fixed.
	//note that fixed_res_guess, which is really a list of fixed nucleosides,
	// should now include the 'moving res' [unless re-specified by user down below].

	// check if this is an 'internal' move.
	if ( moving_partition_res.size() > 1 ) {
		if ( !fixed_res_guess.has_value( rebuild_res ) ) fixed_res_guess.push_back( rebuild_res );
		if ( !floating_base ){
			is_internal = true;
			//moving_res = rebuild_suite; // parin-style old convention, where moving_res = moving_suite for internal suites. Sigh.
		}
	}

	// To specify that the suite moves, we actually need to directly address the movemap... see below.
	suites_that_must_be_minimized.clear();
	if ( rebuild_suite > 0 ) suites_that_must_be_minimized.push_back( rebuild_suite );
	for ( Size i = 1; i <= cutpoint_closed.size(); i++ ) {
		suites_that_must_be_minimized.push_back( cutpoint_closed[i] );
		runtime_assert ( pose.fold_tree().is_cutpoint( cutpoint_closed[i] ) );
	}
	std::string full_sequence = full_model_info.full_sequence();
	Size nres = full_sequence.size();
	if ( full_sequence[ nres - 1] == 'X' ){
		full_sequence = full_sequence.substr( 0, nres - 1 );
		nres -= 1;
	}
	utility::vector1< bool > is_working_res;
	for ( Size n = 1; n <= nres; n++ ) is_working_res.push_back( full_model_info.res_list().has_value( n ) );
	bool const is_prepend = rebuild_res < reference_res;
	std::map< Size, bool > is_prepend_map; // used for rmsd calculations -- include phosphate. May not be necessary.
	utility::vector1< Size > rmsd_res_list_ = rmsd_res_list;
	if ( rmsd_res_list_.size() == 0 ) rmsd_res_list_.push_back( full_model_info.sub_to_full( rebuild_res ) );
	for ( Size i = 1; i <= rmsd_res_list_.size(); i++ ) {
		Size const & n = rmsd_res_list_[ i ];
		if ( is_working_res[ n ] ) {
			is_prepend_map[ n ] = is_prepend; //pose.residue_type( full_model_info.full_to_sub( n ) ).has_variant_type( "VIRTUAL_PHOSPHATE" );
		}
	}

	StepWiseRNA_JobParametersOP job_parameters = new StepWiseRNA_JobParameters;
	job_parameters->set_output_extra_RMSDs( false );
	job_parameters->set_is_simple_full_length_job_params( false );
	job_parameters->set_full_to_sub( full_model_info.full_to_sub() ); // will update sub_to_full
	job_parameters->set_full_sequence( full_sequence ); // will update working_sequence.
	job_parameters->set_moving_res( full_model_info.sub_to_full( moving_res ) );
	job_parameters->set_is_prepend( is_prepend );
	job_parameters->set_is_prepend_map( is_prepend_map );
	job_parameters->set_is_internal( is_internal );
	job_parameters->set_floating_base( floating_base );
	job_parameters->set_floating_base_anchor_res( full_model_info.sub_to_full( floating_base_anchor_res ) );
	job_parameters->set_rebuild_bulge_mode( false );
	job_parameters->set_rebuild_bulge_mode( figure_out_rebuild_bulge_mode( pose, rebuild_res ) );
	job_parameters->set_sample_both_sugar_base_rotamer( figure_out_sample_both_sugar_base_rotamer( pose, floating_base, rebuild_suite ) );
	job_parameters->set_partition_definition( partition_definition );
	job_parameters->set_is_working_res( is_working_res ); // may not be in use later. anyway, set it, just in case.
	job_parameters->set_working_moving_res_list( make_vector1( moving_res ) ); // also sets up working_moving_suite and working_moving_suite_list
	//job_parameters->set_chain_boundaries(); // Not in use later.
	job_parameters->set_five_prime_chain_break_res( five_prime_chain_break_res );
	job_parameters->set_gap_size( gap_size );
	job_parameters->set_working_native_pose( native_pose );
	job_parameters->set_working_fixed_res( fixed_res_guess ); // will get updated later, I think.
	job_parameters->set_rmsd_res_list( rmsd_res_list_ );
	job_parameters->set_terminal_res( terminal_res );
	job_parameters->set_working_terminal_res( full_model_info.full_to_sub( terminal_res ) );
	job_parameters->set_working_moving_partition_pos( moving_partition_res );
	job_parameters->set_input_res_vectors( make_vector1( full_model_info.sub_to_full( input_res1 ),
																											 full_model_info.sub_to_full( input_res2 ) ) );
	job_parameters->set_cutpoint_closed_list( full_model_info.sub_to_full( cutpoint_closed ) );
	job_parameters->set_cutpoint_open_list( full_model_info.sub_to_full( cutpoint_open ) );
	job_parameters->set_working_best_alignment( fixed_res_guess );
	job_parameters->set_working_native_alignment( fixed_res_guess );
	job_parameters->set_native_alignment( full_model_info.sub_to_full( fixed_res_guess ) ); // is this in use?
	job_parameters->set_global_sample_res_list( full_model_info.sub_to_full( make_vector1( moving_res ) ) );
	job_parameters->set_force_syn_chi_res_list( syn_chi_res_list ); // will get automatically converted to working numbering.
	job_parameters->set_fold_tree( pose.fold_tree() );

	TR.Debug << "past job_parameters initialization " << std::endl;

	return job_parameters;
}



/////////////////////////////////////////////////////////////////////////////////////////////////
StepWiseRNA_JobParametersOP
setup_job_parameters_legacy( Size const rebuild_res,
														 pose::Pose const & pose,
														 pose::PoseCOP native_pose,
														 utility::vector1< Size > const & rmsd_res_list,
														 utility::vector1< Size > const & terminal_res,
														 utility::vector1< Size > const & syn_chi_res_list,
														 utility::vector1< Size > & fixed_res_guess,
														 utility::vector1< Size > & suites_that_must_be_minimized ){
	using namespace chemical;
	using namespace pose::full_model_info;

	FullModelInfo const & full_model_info = const_full_model_info( pose );

	TR.Debug << pose.fold_tree() << std::endl;
	TR << "Rebuild residue: " << rebuild_res << std::endl;

	utility::vector1< Size > input_res1, input_res2 /*blank*/, cutpoint_open, cutpoint_closed;
	utility::vector1< Size > moving_res_list = make_vector1( rebuild_res ); // legacy
	Size cutpoint_closed_distal( 0 );
	input_res1 = fixed_res_guess; // not_rebuild_res

	Size rebuild_suite( 0 );
	bool floating_base( false ), is_internal( false );
	bool cut_at_previous =  (rebuild_res == 1) || pose.fold_tree().is_cutpoint( rebuild_res - 1 );
	bool cut_at_next     =  (rebuild_res == pose.total_residue() ) || pose.fold_tree().is_cutpoint( rebuild_res );
	if ( cut_at_next && cut_at_previous ){
		floating_base = true;
	} else {
		cut_at_previous = cut_at_previous && !check_jump_to_previous_residue_in_chain( pose, rebuild_res, moving_res_list );
		cut_at_next     = cut_at_next     && !check_jump_to_next_residue_in_chain( pose, rebuild_res, moving_res_list );
		if ( !cut_at_previous  && cut_at_next ){
			rebuild_suite = rebuild_res - 1;
		} else if ( cut_at_previous && !cut_at_next ){
			rebuild_suite = rebuild_res;
		} else if ( !cut_at_previous && !cut_at_next ){ // internal
			rebuild_suite = rebuild_res;
			is_internal = true;
		}
	}

	// for internal moves, need to be smart about input_res definitions -- 'domains' that are separated by moving residue.
	// this loop will also determine any chainbreak that requires closure
	utility::vector1< bool > partition_definition;
	Size floating_base_anchor_res( 0 );
	if ( floating_base ){
		floating_base_anchor_res = get_anchor_res( rebuild_res, pose );
		partition_definition = get_partition_definition_floating_base( pose, rebuild_res );
	} else {
		partition_definition = get_partition_definition( pose, rebuild_suite );
	}

	utility::vector1< Size > const chains = figure_out_chains_from_full_model_info_const( pose );
	//bool found_moving_cutpoint( false );
	input_res1.clear();
	input_res2.clear();
	for ( Size n = 1; n <= pose.total_residue(); n++ ){
		if ( !partition_definition[ n ] ) {
			input_res1.push_back( n );
		}	else {
			input_res2.push_back( n );
		}
		// look for cutpoints
		if ( n == pose.total_residue() ) continue;
		if ( pose.fold_tree().is_cutpoint( n ) && ( partition_definition[ n ] != partition_definition[ n+1 ] ) ){
			if ( n == rebuild_res && n+1 == floating_base_anchor_res ) continue;
			if ( n == floating_base_anchor_res && ( n + 1 == rebuild_res) ) continue;
			//				if ( !floating_base ) runtime_assert( !found_moving_cutpoint );
			//found_moving_cutpoint = true;
			if ( pose.residue_type( n ).has_variant_type( CUTPOINT_LOWER ) ){
				runtime_assert( pose.residue_type( n + 1 ).has_variant_type( CUTPOINT_UPPER  ) );
				cutpoint_closed.push_back( n );
				if ( !floating_base ||
						 ( floating_base_anchor_res > rebuild_res &&  rebuild_res >  n ) ||
						 ( floating_base_anchor_res < rebuild_res &&  rebuild_res <= n ) ){
					if ( cutpoint_closed_distal == 0 ) cutpoint_closed_distal = n;
				}
			} else if ( chains[ n ] != chains[ n+1 ]  ){
				cutpoint_open.push_back( n );
			}
		}
	}

	TR << "INPUT_RES1 " << make_tag_with_dashes(input_res1) << " INPUT_RES2 " << make_tag_with_dashes(input_res2) << " REBUILD_RES " << rebuild_res << " REBUILD_SUITE " <<  rebuild_suite << " CUTPOINT_CLOSED " << cutpoint_closed << " CUTPOINT_CLOSED_DISTAL " << cutpoint_closed_distal << " FLOATING_BASE " << floating_base << std::endl;

	// if there's really just a single nucleotide being rebuilt, its nucleoside is not fixed.
	// but if there's a whole chunk of stuff, its sugar & base are assumed fixed.
	//note that fixed_res_guess, which is really a list of fixed nucleosides,
	// should now include the 'moving res' [unless re-specified by user down below].

	// check if this is an 'internal' move.
	if ( input_res1.size() > 1 && input_res2.size() > 1 ) {
		fixed_res_guess.push_back( rebuild_res );
		is_internal = true;
	}

	// To specify that the suite moves, we actually need to directly address the movemap... see below.
	suites_that_must_be_minimized.clear();
	if ( rebuild_suite > 0 ) suites_that_must_be_minimized.push_back( rebuild_suite );
	if ( cutpoint_closed.size() > 0 ) {
		for ( Size i = 1; i <= cutpoint_closed.size(); i++ ) {
			suites_that_must_be_minimized.push_back( cutpoint_closed[i] );
			if ( !pose.fold_tree().is_cutpoint( cutpoint_closed[i] ) ) utility_exit_with_message( "StepWiseRNA requires a chainbreak right at sampled residue" );
		}
	}

	std::string full_sequence = full_model_info.full_sequence();
	Size nres = pose.total_residue();
	if ( full_sequence[ nres - 1] == 'X' ){
		full_sequence = full_sequence.substr( 0, nres - 1 );
		nres -= 1;
	}
	legacy::StepWiseRNA_JobParametersSetup stepwise_rna_job_parameters_setup( full_model_info.sub_to_full( moving_res_list ),
																																						full_sequence,
																																						full_model_info.sub_to_full( input_res1 ),
																																						full_model_info.sub_to_full( input_res2 ),
																																						full_model_info.sub_to_full( cutpoint_open ),
																																						full_model_info.sub_to_full( cutpoint_closed_distal ) );
	// following still might be worth doing -- just didn't seem necessary at this point.
	//		stepwise_rna_job_parameters_setup.set_cutpoint_closed_list( full_model_info.sub_to_full( cutpoint_closed ) );
	stepwise_rna_job_parameters_setup.set_fixed_res( full_model_info.sub_to_full( fixed_res_guess ) );
	stepwise_rna_job_parameters_setup.set_floating_base( floating_base );
	if ( is_internal ) stepwise_rna_job_parameters_setup.set_force_internal( true );
	if ( floating_base ){
		stepwise_rna_job_parameters_setup.set_assert_jump_point_in_fixed_res( false );
		stepwise_rna_job_parameters_setup.set_floating_base_anchor_res( full_model_info.sub_to_full( floating_base_anchor_res ) );
	}
	if ( rmsd_res_list.size() > 0 ) stepwise_rna_job_parameters_setup.set_rmsd_res_list( rmsd_res_list /*global numbering*/ );
	else  stepwise_rna_job_parameters_setup.set_rmsd_res_list( full_model_info.sub_to_full( make_vector1( rebuild_res ) ) );

	stepwise_rna_job_parameters_setup.set_rebuild_bulge_mode( figure_out_rebuild_bulge_mode( pose, rebuild_res ) );
	stepwise_rna_job_parameters_setup.set_sample_both_sugar_base_rotamer( figure_out_sample_both_sugar_base_rotamer( pose, floating_base, rebuild_suite ) );

	utility::vector1< std::string > jump_point_pair_list;
	kinematics::FoldTree const & f = pose.fold_tree();
	for ( Size n = 1; n <= f.num_jump(); n++ ){
		jump_point_pair_list.push_back( ObjexxFCL::string_of( full_model_info.sub_to_full( f.upstream_jump_residue( n ) ) ) + "-" +
																		ObjexxFCL::string_of( full_model_info.sub_to_full( f.downstream_jump_residue( n ) ) ) );
	}
	stepwise_rna_job_parameters_setup.set_jump_point_pair_list( jump_point_pair_list ); //Important!: Needs to be called after set_fixed_res

	utility::vector1< std::string > alignment_res; //why is this a string vector?????
	for ( Size n = 1; n <= fixed_res_guess.size(); n++ ) alignment_res.push_back( ObjexxFCL::string_of( full_model_info.sub_to_full( fixed_res_guess[ n ] ) ) );
	stepwise_rna_job_parameters_setup.set_alignment_res( alignment_res );
	stepwise_rna_job_parameters_setup.set_native_alignment_res( full_model_info.sub_to_full( fixed_res_guess ) );
	stepwise_rna_job_parameters_setup.set_terminal_res( terminal_res );

	// could use this later to minimize more residues...
	//stepwise_rna_job_parameters_setup.set_global_sample_res_list( option[ global_sample_res_list ]() ); //March 20, 2011

	// NOT SURE ABOUT THIS. false by default, but shows up as true in 'normal' erraser runs.
	stepwise_rna_job_parameters_setup.set_allow_chain_boundary_jump_partner_right_at_fixed_BP( true );

	// NOT SURE ABOUT THIS...
	stepwise_rna_job_parameters_setup.set_add_virt_res_as_root( true );

	// ignore fold tree setup which is hopelessly complicated in stepwise_rna_job_parameters_setup.
	stepwise_rna_job_parameters_setup.force_fold_tree( f );
	stepwise_rna_job_parameters_setup.apply();
	StepWiseRNA_JobParametersOP job_parameters = stepwise_rna_job_parameters_setup.job_parameters();
	job_parameters->set_working_native_pose( native_pose );
	job_parameters->set_force_syn_chi_res_list( syn_chi_res_list ); // will get automatically converted to working numbering.

	if ( is_internal ) runtime_assert( job_parameters->is_internal() );
	TR.Debug << "past job_parameters initialization " << std::endl;

	return job_parameters;
}


/////////////////////////////////////////////////////////////////////////
bool
figure_out_rebuild_bulge_mode( pose::Pose const & pose, Size const rebuild_res ){
	kinematics::FoldTree const & f = pose.fold_tree();
	if ( rebuild_res > 1 &&
			 pose.residue( rebuild_res - 1 ).has_variant_type( "VIRTUAL_RIBOSE" ) &&
			 ( !f.is_cutpoint( rebuild_res - 1 ) ||
				 is_cutpoint_closed( pose, rebuild_res - 1 ) ) &&
			 !f.is_cutpoint( rebuild_res ) &&
			 f.jump_nr( rebuild_res - 1, rebuild_res + 1) > 0 ) return true;
	if ( rebuild_res < pose.total_residue() &&
			 pose.residue( rebuild_res + 1 ).has_variant_type( "VIRTUAL_RIBOSE" ) &&
			 ( !f.is_cutpoint( rebuild_res ) ||
				 is_cutpoint_closed( pose, rebuild_res ) ) &&
			 !f.is_cutpoint( rebuild_res - 1 ) &&
			 f.jump_nr( rebuild_res - 1, rebuild_res + 1 ) > 0 ) return true;
	return false;
}

/////////////////////////////////////////////////////////////////////////
bool
figure_out_sample_both_sugar_base_rotamer( pose::Pose const & pose, bool const floating_base, Size const rebuild_suite ){
	if ( !floating_base &&
			 rotamer_sampler::rna::sampling_sugar_at_five_prime(  pose, rebuild_suite ) &&
			 rotamer_sampler::rna::sampling_sugar_at_three_prime( pose, rebuild_suite ) &&
			 !pose.residue( rebuild_suite ).has_variant_type( "VIRTUAL_RIBOSE" ) &&
			 !pose.residue( rebuild_suite + 1 ).has_variant_type( "VIRTUAL_RIBOSE" ) ) {
		TR << "SAMPLE_BOTH_SUGAR_BASE_ROTAMER!" << std::endl;
		return true;
	}
	return false;
}

} //rna
} //sampling
} //stepwise
} //protocols
