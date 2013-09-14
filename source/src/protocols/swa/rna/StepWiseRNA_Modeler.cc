// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/swa/rna/SWA_Modeler.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu

#include <protocols/swa/rna/StepWiseRNA_Modeler.hh>
#include <protocols/swa/rna/StepWiseRNA_ResidueSampler.hh>
#include <protocols/swa/rna/StepWiseRNA_BaseCentroidScreener.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParametersSetup.hh>
#include <protocols/swa/rna/StepWiseRNA_Minimizer.hh>
#include <protocols/swa/rna/StepWiseRNA_BaseCentroidScreener.hh>
#include <protocols/swa/StepWiseUtil.hh>

#include <core/chemical/VariantType.hh>
#include <core/id/TorsionID.hh>

#include <core/types.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh> // should remove this later.

#include <utility/tools/make_vector1.hh>
#include <basic/Tracer.hh>

#include <ObjexxFCL/string.functions.hh>

static basic::Tracer TR( "protocols.swa.rna.StepWiseRNA_Modeler" );

using utility::tools::make_vector1;

namespace protocols {
namespace swa {
namespace rna {

//Constructor
StepWiseRNA_Modeler::StepWiseRNA_Modeler( core::Size const sample_res,
																	core::scoring::ScoreFunctionOP scorefxn ) :
	moving_res_( make_vector1( sample_res ) ),
	scorefxn_( scorefxn ),
	silent_file_( "" ),
	sampler_num_pose_kept_( 108 ),
	num_pose_minimize_( 99999 ),
	num_sampled_( 0 ),
	sampler_native_screen_rmsd_cutoff_( 2.0 ),
	cluster_rmsd_( 0.5 ),
	native_edensity_score_cutoff_( -1 ),
	sampler_native_rmsd_screen_( false ),
	o2star_screen_( true ),
	verbose_( false ),
	distinguish_pucker_( true ),
	finer_sampling_at_chain_closure_( false ),
	PBP_clustering_at_chain_closure_( false ),
	allow_syn_pyrimidine_( false ),
	extra_chi_( false ),
	use_phenix_geo_( false ),
	kic_sampling_( false ),
	kic_sampling_if_relevant_( false ),
	centroid_screen_( true ),
	VDW_atr_rep_screen_( true ),
	force_centroid_interaction_( false ),
	choose_random_( false ),
	num_random_samples_( 1 ),
	skip_sampling_( false ),
	perform_minimize_( true ),
	minimize_and_score_sugar_( true ),
	minimize_and_score_native_pose_( false ),
	rm_virt_phosphate_( false ),
	VDW_rep_alignment_RMSD_CUTOFF_( 0.001 ),
	output_pdb_( false ),
	output_minimized_pose_data_list_( false ),
	// following are options that are not in SWA_Modeler (yet) and need to be treated carefully during integration.
	VDW_rep_screen_physical_pose_clash_dist_cutoff_( false ),
	integration_test_mode_( false ),
	allow_bulge_at_chainbreak_( false ),
	parin_favorite_output_( false ),
	floating_base_( false ),
	reinitialize_CCD_torsions_( false ),
	sampler_extra_epsilon_rotamer_( false ),
	sampler_extra_beta_rotamer_( false ),
	sample_both_sugar_base_rotamer_( false ),
	sampler_include_torsion_value_in_tag_( false ),
	rebuild_bulge_mode_( false ),
	combine_long_loop_mode_( false ),
	do_not_sample_multiple_virtual_sugar_( false ),
	sample_ONLY_multiple_virtual_sugar_( false ),
	sampler_assert_no_virt_ribose_sampling_( false ),
	allow_base_pair_only_centroid_screen_( false ),
	minimizer_perform_o2star_pack_( false ),
	minimizer_output_before_o2star_pack_( false ),
	minimizer_rename_tag_( false )
{}

//Destructor
StepWiseRNA_Modeler::~StepWiseRNA_Modeler()
{}

/////////////////////
std::string
StepWiseRNA_Modeler::get_name() const {
	return "StepWiseRNA_Modeler";
}
//////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Modeler::apply( core::pose::Pose & pose ){

	using namespace core::pose;
	using namespace core::chemical;
	using namespace core::kinematics;
	using namespace core::scoring;
	using namespace protocols::swa::rna;

	if ( !job_parameters_ ) job_parameters_ = setup_job_parameters_for_swa( moving_res_, pose );
	if ( !job_parameters_ ) utility_exit_with_message( "You must supply job parameters!" );

	StepWiseRNA_ResidueSampler stepwise_rna_residue_sampler( job_parameters_ );
	stepwise_rna_residue_sampler.set_silent_file ( silent_file_ + "_sampling" );
	stepwise_rna_residue_sampler.set_scorefxn ( scorefxn_ );
	stepwise_rna_residue_sampler.set_num_pose_kept ( sampler_num_pose_kept_ );
	stepwise_rna_residue_sampler.set_native_rmsd_screen ( sampler_native_rmsd_screen_ );
	stepwise_rna_residue_sampler.set_native_screen_rmsd_cutoff ( sampler_native_screen_rmsd_cutoff_ );
	stepwise_rna_residue_sampler.set_perform_o2star_pack ( o2star_screen_ );
	stepwise_rna_residue_sampler.set_verbose ( verbose_ );
	stepwise_rna_residue_sampler.set_cluster_rmsd (	cluster_rmsd_	);
	stepwise_rna_residue_sampler.set_distinguish_pucker ( distinguish_pucker_ );
	stepwise_rna_residue_sampler.set_finer_sampling_at_chain_closure ( finer_sampling_at_chain_closure_ );
	stepwise_rna_residue_sampler.set_PBP_clustering_at_chain_closure ( PBP_clustering_at_chain_closure_ );
	stepwise_rna_residue_sampler.set_allow_syn_pyrimidine( allow_syn_pyrimidine_ );
	stepwise_rna_residue_sampler.set_extra_chi( extra_chi_ );
	stepwise_rna_residue_sampler.set_use_phenix_geo ( use_phenix_geo_  );
	stepwise_rna_residue_sampler.set_kic_sampling( kic_sampling_ );
	stepwise_rna_residue_sampler.set_centroid_screen ( centroid_screen_ );
	stepwise_rna_residue_sampler.set_VDW_atr_rep_screen ( VDW_atr_rep_screen_ );
	stepwise_rna_residue_sampler.set_force_centroid_interaction ( force_centroid_interaction_ );

	stepwise_rna_residue_sampler.set_integration_test_mode( integration_test_mode_ ); //Should set after setting sampler_native_screen_rmsd_cutoff, fast, medium_fast options.
	stepwise_rna_residue_sampler.set_allow_bulge_at_chainbreak( allow_bulge_at_chainbreak_ );
	stepwise_rna_residue_sampler.set_parin_favorite_output( parin_favorite_output_ );
	stepwise_rna_residue_sampler.set_floating_base( floating_base_ );
	stepwise_rna_residue_sampler.set_reinitialize_CCD_torsions( reinitialize_CCD_torsions_ );
	stepwise_rna_residue_sampler.set_extra_epsilon_rotamer( sampler_extra_epsilon_rotamer_ );
	stepwise_rna_residue_sampler.set_extra_beta_rotamer( sampler_extra_beta_rotamer_ );
	stepwise_rna_residue_sampler.set_sample_both_sugar_base_rotamer( sample_both_sugar_base_rotamer_ ); //Nov 12, 2010
	stepwise_rna_residue_sampler.set_include_torsion_value_in_tag( sampler_include_torsion_value_in_tag_ );
	stepwise_rna_residue_sampler.set_rebuild_bulge_mode( rebuild_bulge_mode_ );
	stepwise_rna_residue_sampler.set_combine_long_loop_mode( combine_long_loop_mode_ );
	stepwise_rna_residue_sampler.set_do_not_sample_multiple_virtual_sugar( do_not_sample_multiple_virtual_sugar_ );
	stepwise_rna_residue_sampler.set_sample_ONLY_multiple_virtual_sugar( sample_ONLY_multiple_virtual_sugar_ );
	stepwise_rna_residue_sampler.set_assert_no_virt_ribose_sampling( sampler_assert_no_virt_ribose_sampling_ );
	stepwise_rna_residue_sampler.set_allow_base_pair_only_centroid_screen( allow_base_pair_only_centroid_screen_ );

	StepWiseRNA_BaseCentroidScreenerOP base_centroid_screener = new StepWiseRNA_BaseCentroidScreener ( pose, job_parameters_ );
	stepwise_rna_residue_sampler.set_base_centroid_screener ( base_centroid_screener );

	StepWiseRNA_VDW_BinScreenerOP user_input_VDW_bin_screener = new StepWiseRNA_VDW_BinScreener();
	if ( VDW_rep_screen_info_.size() > 0 ) {
		user_input_VDW_bin_screener->set_VDW_rep_alignment_RMSD_CUTOFF ( VDW_rep_alignment_RMSD_CUTOFF_ );
		user_input_VDW_bin_screener->set_VDW_rep_delete_matching_res( VDW_rep_delete_matching_res_ );
		user_input_VDW_bin_screener->set_physical_pose_clash_dist_cutoff( VDW_rep_screen_physical_pose_clash_dist_cutoff_ );
		user_input_VDW_bin_screener->setup_using_user_input_VDW_pose( VDW_rep_screen_info_, pose, job_parameters_ );
		user_input_VDW_bin_screener->set_output_pdb( output_pdb_ );
	}
	stepwise_rna_residue_sampler.set_user_input_VDW_bin_screener ( user_input_VDW_bin_screener );

	if ( choose_random_ ){
		stepwise_rna_residue_sampler.set_choose_random( true );
		stepwise_rna_residue_sampler.set_num_random_samples( num_random_samples_ );
		stepwise_rna_residue_sampler.set_cluster_rmsd( 0.0 ); // don't cluster.
	}
	print_JobParameters_info( job_parameters_, "job_parameters_COP", TR.Debug );

	if ( ! skip_sampling_ ) stepwise_rna_residue_sampler.apply( pose );

	utility::vector1< pose_data_struct2 > & pose_data_list = stepwise_rna_residue_sampler.get_pose_data_list();
	num_sampled_ = pose_data_list.size();
	if ( num_sampled_ == 0 && !skip_sampling_ ){
		TR << "WARNING! WARNING! WARNING! pose_data_list.size() == 0! " << std::endl;
		if ( !output_minimized_pose_data_list_ ) return; // don't do a minimize...
	}

	if ( skip_sampling_ /* still want to minimize*/ ){
		pose_data_struct2 data_struct;
		data_struct.pose_OP = new Pose( pose );
		data_struct.score = 0.0;
		data_struct.tag = "input_pose";
		pose_data_list.push_back( data_struct );
	}

	if ( minimize_and_score_native_pose_ ) {
		if ( !get_native_pose() ) utility_exit_with_message ( "minimize_and_score_native_pose == True but user did not pass in native pose" );
		pose_data_struct2 native_data_struct;
		native_data_struct.pose_OP = new Pose( *get_native_pose() ); //hopefully this clones...
		native_data_struct.score = 0.0;
		native_data_struct.tag = "working_native_pose";
		pose_data_list.push_back ( native_data_struct );
	}

	// let's output the final pose_data_list, just to have a look
	if ( verbose_ ) stepwise_rna_residue_sampler.output_pose_data_list ( silent_file_ + "_final_sample" );

	////////////////////////////////////////////////////////////////
	stepwise_rna_minimizer_ = new StepWiseRNA_Minimizer( stepwise_rna_residue_sampler.get_pose_data_list(), job_parameters_ );
	stepwise_rna_minimizer_->set_silent_file ( silent_file_ );
	stepwise_rna_minimizer_->set_verbose (  verbose_ );
	stepwise_rna_minimizer_->set_scorefxn ( scorefxn_ );
	stepwise_rna_minimizer_->set_centroid_screen (  centroid_screen_ );
	stepwise_rna_minimizer_->set_base_centroid_screener ( base_centroid_screener );
	stepwise_rna_minimizer_->set_perform_minimize(  perform_minimize_ );
	stepwise_rna_minimizer_->set_native_rmsd_screen (  sampler_native_rmsd_screen_ );
	stepwise_rna_minimizer_->set_native_edensity_score_cutoff ( native_edensity_score_cutoff_ );
	stepwise_rna_minimizer_->set_rm_virt_phosphate (  rm_virt_phosphate_ );
	stepwise_rna_minimizer_->set_native_screen_rmsd_cutoff (  sampler_native_screen_rmsd_cutoff_ + 1 ); //+1 for leniency Sept 20, 2010

	if ( integration_test_mode_ ) num_pose_minimize_ = 1;
	if ( num_pose_minimize_ > 0 ) stepwise_rna_minimizer_->set_num_pose_minimize ( num_pose_minimize_ );
	stepwise_rna_minimizer_->set_minimize_and_score_sugar ( minimize_and_score_sugar_ );
	stepwise_rna_minimizer_->set_user_input_VDW_bin_screener ( user_input_VDW_bin_screener );
	stepwise_rna_minimizer_->set_output_minimized_pose_data_list( output_minimized_pose_data_list_ );
	if ( minimize_move_map_ ) stepwise_rna_minimizer_->set_move_map_list( make_vector1( *minimize_move_map_ ) );
	stepwise_rna_minimizer_->set_perform_o2star_pack(  minimizer_perform_o2star_pack_ );
	stepwise_rna_minimizer_->set_output_before_o2star_pack( minimizer_output_before_o2star_pack_ );
	stepwise_rna_minimizer_->set_rename_tag( minimizer_rename_tag_ );

	stepwise_rna_minimizer_->apply ( pose );

	// Need to make sure that final pose output is the lowest scoring one. Is it?

}


///////////////////////////////////////////////////////////////////////////////////////////////////
// This could go into a Util.hh if it ends up being more useful.
// Briefly, we need to make a StepWiseRNA_JobParameters object that will
// be fed into various StepWiseRNA movers. Setting this up can be quite complicated...
// it has become a grab bag of residue lists referring to the global pose, the working pose,
// sequence mappings, "Is_Prepend_map", etc.
//
StepWiseRNA_JobParametersOP
StepWiseRNA_Modeler::setup_job_parameters_for_swa( utility::vector1< Size > moving_res, core::pose::Pose const & pose ){

	using namespace core::pose;
	using namespace core::chemical;
	using namespace core::id;
	using namespace protocols::swa;

	if ( moving_res.size() != 1 ) utility_exit_with_message( "For now, StepWiseRNA_Modeler requires exactly 1 number in -moving_res unless you feed it a JobParameters object." );
	Size const rebuild_res =  moving_res[1];

	std::string full_sequence = pose.sequence();
	Size nres = pose.total_residue();

	// what if there is a virtual residue? need to remove it, actually, before running stepwise_rna_job_parameters_setup.
	if ( full_sequence[nres - 1] == 'X' ){
		full_sequence = full_sequence.substr( 0, nres - 1 );
		nres -= 1;
	}

	TR.Debug << pose.fold_tree() << std::endl;
	TR.Debug << "Rebuild residue: " << rebuild_res << std::endl;

	utility::vector1< Size > not_rebuild_res;
	for ( Size n = 1; n <= nres; n++ ) if ( n != rebuild_res ) not_rebuild_res.push_back( n );
	utility::vector1< Size > input_res1, input_res2 /*blank*/, cutpoint_open;
	input_res1 = not_rebuild_res;
	utility::vector1< Size > fixed_res_guess = not_rebuild_res; // may be revised below.

	Size cutpoint_closed( 0 );
	utility::vector1< Size > suites_that_must_be_minimized;
	// check for cutpoint variant.
	if ( pose.residue_type( rebuild_res ).has_variant_type( CUTPOINT_UPPER ) ){
		runtime_assert( pose.residue_type( rebuild_res - 1 ).has_variant_type( CUTPOINT_LOWER  ) );
		cutpoint_closed = rebuild_res - 1;
	} else if (  pose.residue_type( rebuild_res ).has_variant_type( CUTPOINT_LOWER ) ){
		runtime_assert( pose.residue_type( rebuild_res + 1 ).has_variant_type( CUTPOINT_UPPER  ) );
		cutpoint_closed = rebuild_res;
	} else if ( rebuild_res == 1  ){
		/*must be a prepend!*/;
	} else if ( rebuild_res == nres  ){
		/*must be an append!*/;
	} else if ( pose.fold_tree().is_cutpoint( rebuild_res - 1 ) ){
		cutpoint_open.push_back( rebuild_res - 1 );
	} else if ( pose.fold_tree().is_cutpoint( rebuild_res ) ){
		cutpoint_open.push_back( rebuild_res );
	} else {
		// internal. need to be smart about input_res definitions -- 'domains' that are separated by moving residue.
		input_res1.clear();
		input_res2.clear();
		utility::vector1< bool > partition_definition =	get_partition_definition( pose, rebuild_res );
		bool found_moving_cutpoint( false );
		for ( Size n = 1; n <= pose.total_residue(); n++ ){
			if ( !partition_definition[ n ] ) {
				input_res1.push_back( n );
			}	else {
				input_res2.push_back( n );
			}
			// look for cutpoints
			if ( n == pose.total_residue() ) continue;
			if ( pose.fold_tree().is_cutpoint( n ) && ( partition_definition[ n ] != partition_definition[ n+1 ] ) ){
				runtime_assert( !found_moving_cutpoint );
				found_moving_cutpoint = true;
				if ( pose.residue_type( n ).has_variant_type( CUTPOINT_LOWER ) ){
					runtime_assert( pose.residue_type( n + 1 ).has_variant_type( CUTPOINT_UPPER  ) );
					if ( cutpoint_closed > 0 ) utility_exit_with_message( "Right now, StepWiseRNA Modeler can only handle poses with single movable closable cutpoints!" );
					cutpoint_closed = n;
				} else {
					cutpoint_open.push_back( n );
				}
			}
		}
		//note that fixed_res_guess, which is really a list of fixed nucleosides,
		// should now include the 'moving res' [unless re-specified by user down below].
		fixed_res_guess.push_back( rebuild_res );

		// To specify that the suite move, we actually need to directly address the movemap... see below.
		suites_that_must_be_minimized.push_back( rebuild_res );
		if ( cutpoint_closed > 0 ) suites_that_must_be_minimized.push_back( cutpoint_closed );

		// last, but not least, there might be some information in the domain map. Note
		// that generally we could instead replace fixed_res with an inputted domain map.
		// that is, get rid of fixed_res_ & minimize_res_ and instead have a local fixed_domain_map,
		// which can instead be updated by set_fixed_res.
		using namespace core::pose::full_model_info;
		// following is dangerous -- what if full_model_info is not set properly?
		FullModelInfo const & full_model_info = const_full_model_info_from_pose( pose );
		utility::vector1< Size > const & fixed_domain_map = full_model_info.fixed_domain_map();
		utility::vector1< Size > const & res_list = full_model_info.res_list();
		for ( Size n = 1; n < pose.total_residue(); n++ ){
			if ( !pose.fold_tree().is_cutpoint(n) &&
					 ( res_list[ n + 1 ]  == res_list[ n ] + 1 ) &&
					 fixed_res_guess.has_value( n )  && fixed_res_guess.has_value( n+1 )  &&
					 ( fixed_domain_map[ res_list[ n + 1 ] ] !=  fixed_domain_map[ res_list[ n ] ] ) &&
					 !suites_that_must_be_minimized.has_value( n ) ){
				TR.Debug << "ADDING NEW SUITE TO BE MINIMIZED BASED ON LOCATION AT DOMAIN BOUNDARY: " << n << std::endl;
				suites_that_must_be_minimized.push_back( n );
			}

		}

	}


	if ( cutpoint_closed > 0 ) {
		if ( kic_sampling_if_relevant_ ) kic_sampling_ = true;
		if ( !pose.fold_tree().is_cutpoint( cutpoint_closed ) ) utility_exit_with_message( "StepWiseRNA requires a chainbreak right at sampled residue" );
		if ( rebuild_res == 1 || rebuild_res == pose.total_residue() ) utility_exit_with_message( "StepWiseRNA requires that residue is not at terminus!");
	}

	StepWiseRNA_JobParametersSetup stepwise_rna_job_parameters_setup( moving_res,
																																		 full_sequence,
																																		 input_res1,
																																		 input_res2,
																																		 cutpoint_open,
																																		 cutpoint_closed );

	stepwise_rna_job_parameters_setup.set_fixed_res( fixed_res_guess );

	utility::vector1< Size > rmsd_res_list;
	rmsd_res_list.push_back( rebuild_res );
	stepwise_rna_job_parameters_setup.set_rmsd_res_list( rmsd_res_list );

	// not sure about the following -- instead, how about reading jump residues from within pose itself?
	if ( cutpoint_closed > 0 || cutpoint_open.size() > 0 ){
		utility::vector1< std::string > jump_point_pair_list;
		jump_point_pair_list.push_back( ObjexxFCL::string_of( rebuild_res - 1 ) + "-" + ObjexxFCL::string_of( rebuild_res + 1 ) );
		stepwise_rna_job_parameters_setup.set_jump_point_pair_list( jump_point_pair_list ); //Important!: Needs to be called after set_fixed_res
	}

	utility::vector1< std::string > alignment_res; //why is this a string vector?????
	for ( Size n = 1; n <= fixed_res_guess.size(); n++ ) alignment_res.push_back( ObjexxFCL::string_of( fixed_res_guess[ n ] ) );
	stepwise_rna_job_parameters_setup.set_alignment_res( alignment_res );
	stepwise_rna_job_parameters_setup.set_native_alignment_res( fixed_res_guess );

	// could use this later to minimize more residues...
	//stepwise_rna_job_parameters_setup.set_global_sample_res_list( option[ global_sample_res_list ]() ); //March 20, 2011

	// NOT SURE ABOUT THIS. false by default, but shows up as true in 'normal' erraser runs.
	stepwise_rna_job_parameters_setup.set_allow_chain_boundary_jump_partner_right_at_fixed_BP ( true );

	// NOT SURE ABOUT THIS...
	stepwise_rna_job_parameters_setup.set_add_virt_res_as_root( true );

	// ignore fold tree setup which is hopelessly complicated in stepwise_rna_job_parameters_setup.
	stepwise_rna_job_parameters_setup.force_fold_tree( pose.fold_tree() );

	stepwise_rna_job_parameters_setup.apply();

	StepWiseRNA_JobParametersOP job_parameters = stepwise_rna_job_parameters_setup.job_parameters();

	job_parameters->set_working_native_pose( get_native_pose() );

	// user input minimize_res...
	if ( minimize_res_.size() > 0 ) { // specifying more residues which could move during the minimize step.
		fixed_res_.clear();
		for ( Size n = 1; n <= nres; n++ ) {
			if ( !minimize_res_.has_value( n ) )	fixed_res_.push_back( n );
		}
	} else {
		for ( Size n = 1; n <= nres; n++ ) {
			if ( !fixed_res_guess.has_value( n ) )	minimize_res_.push_back( n );
		}
	}

	if ( fixed_res_.size() > 0 ) 	job_parameters->set_working_fixed_res( fixed_res_ ); // is this necessary if we just supply movemap?

	minimize_move_map_ = new core::kinematics::MoveMap;
	figure_out_swa_rna_movemap( *minimize_move_map_, pose, minimize_res_ );
	for ( Size n = 1; n <= suites_that_must_be_minimized.size(); n++ ){
		Size const suite_num = suites_that_must_be_minimized[ n ];
		minimize_move_map_->set( TorsionID( suite_num,   id::BB, 5 ), true ); // epsilon
		minimize_move_map_->set( TorsionID( suite_num,   id::BB, 6 ), true ); // zeta
		minimize_move_map_->set( TorsionID( suite_num+1, id::BB, 1 ), true ); // alpha
		minimize_move_map_->set( TorsionID( suite_num+1, id::BB, 2 ), true ); // beta
		minimize_move_map_->set( TorsionID( suite_num+1, id::BB, 3 ), true ); // gamma
	}


	return job_parameters;
}


///////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Modeler::set_job_parameters( StepWiseRNA_JobParametersCOP job_parameters ){ job_parameters_ = job_parameters;	}

///////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Modeler::set_native_pose( core::pose::PoseCOP native_pose ){ native_pose_ = native_pose; }

///////////////////////////////////////////////////////////////////////////////
core::pose::PoseCOP
StepWiseRNA_Modeler::get_native_pose(){ return native_pose_; }


///////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Modeler::output_pose(
																 pose::Pose & pose,
																 std::string const & out_tag,
																 std::string const out_silent_file ) const {
	runtime_assert(  stepwise_rna_minimizer_ != 0 );
	stepwise_rna_minimizer_->output_pose_data_wrapper( out_tag, pose, out_silent_file );
}


} //rna
} //swa
} //protocols
