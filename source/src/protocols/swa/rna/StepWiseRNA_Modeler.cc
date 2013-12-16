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
#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParametersSetup.hh>
#include <protocols/swa/rna/StepWiseRNA_Minimizer.hh>
#include <protocols/swa/rna/screener/StepWiseRNA_BaseCentroidScreener.hh>
#include <protocols/swa/rna/screener/StepWiseRNA_VDW_BinScreener.hh>
#include <protocols/swa/StepWiseUtil.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>

#include <core/chemical/VariantType.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/util.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelInfoUtil.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>

#include <ObjexxFCL/string.functions.hh>

static basic::Tracer TR( "protocols.swa.rna.StepWiseRNA_Modeler" );

using utility::tools::make_vector1;

////////////////////////////////////////////////////////////////////////////////
//
// This is meant to be a simple 'master interface' to StepWiseAssembly and
//  StepWiseMonteCarlo functions, requiring only a pose and the residue to be
//  sampled.
//
// All of the complexities of setting up the JobParameters, etc. are hidden
//  in a single wrapper function below.
//
//  -- Rhiju
////////////////////////////////////////////////////////////////////////////////


namespace protocols {
namespace swa {
namespace rna {

//Constructor
StepWiseRNA_Modeler::StepWiseRNA_Modeler( core::scoring::ScoreFunctionOP scorefxn ) :
	scorefxn_( scorefxn )
{
	initialize_variables();
}


StepWiseRNA_Modeler::StepWiseRNA_Modeler( core::Size const sample_res,
																	core::scoring::ScoreFunctionOP scorefxn ) :
	moving_res_list_( make_vector1( sample_res ) ),
	scorefxn_( scorefxn )
{
	initialize_variables();
}

// If you add a variable, initialize it here, and include in operator= definition below!
void
StepWiseRNA_Modeler::initialize_variables(){
	silent_file_ = "";
	sampler_num_pose_kept_ = 108;
	num_pose_minimize_ = 99999;
	num_sampled_ = 0;
	sampler_native_screen_rmsd_cutoff_ = 2.0;
	cluster_rmsd_ = 0.5;
	native_edensity_score_cutoff_ = -1;
	sampler_native_rmsd_screen_ = false;
	o2prime_screen_ = true;
	verbose_ = false;
	distinguish_pucker_ = true;
	finer_sampling_at_chain_closure_ = false;
	PBP_clustering_at_chain_closure_ = false;
	allow_syn_pyrimidine_ = false;
	extra_chi_ = false;
	use_phenix_geo_ = false;
	virtual_sugar_legacy_mode_ = false;
	kic_sampling_ = false;
	kic_sampling_if_relevant_ = false;
	centroid_screen_ = true;
	VDW_atr_rep_screen_ = true;
	force_centroid_interaction_ = false;
	choose_random_ = false;
	num_random_samples_ = 1;
	skip_sampling_ = false;
	perform_minimize_ = true;
	minimize_and_score_sugar_ = true;
	minimize_and_score_native_pose_ = false;
	rm_virt_phosphate_ = false;
	VDW_rep_alignment_RMSD_CUTOFF_ = 0.001;
	output_pdb_ = false;
	output_minimized_pose_list_ = false;
	VDW_rep_screen_physical_pose_clash_dist_cutoff_ = false;
	integration_test_mode_ = false;
	allow_bulge_at_chainbreak_ = false;
	parin_favorite_output_ = false;
	reinitialize_CCD_torsions_ = false;
	sampler_extra_epsilon_rotamer_ = false;
	sampler_extra_beta_rotamer_ = false;
	sample_both_sugar_base_rotamer_ = false;
	sampler_include_torsion_value_in_tag_ = false;
	combine_long_loop_mode_ = false;
	do_not_sample_multiple_virtual_sugar_ = false;
	sample_ONLY_multiple_virtual_sugar_ = false;
	sampler_assert_no_virt_sugar_sampling_ = false;
	allow_base_pair_only_centroid_screen_ = false;
	minimizer_perform_o2prime_pack_ = false;
	minimizer_output_before_o2prime_pack_ = false;
	minimizer_rename_tag_ = false;
	minimizer_allow_variable_bond_geometry_ = false;
	minimizer_vary_bond_geometry_frequency_ = 0.0;
}

//Destructor
StepWiseRNA_Modeler::~StepWiseRNA_Modeler()
{}

StepWiseRNA_Modeler::StepWiseRNA_Modeler( StepWiseRNA_Modeler const & src ):
	Mover( src )
{
	*this = src;
}

/// @brief clone the conformation
StepWiseRNA_ModelerOP
StepWiseRNA_Modeler::clone_modeler() const
{
	return new StepWiseRNA_Modeler( *this );
}

StepWiseRNA_Modeler &
StepWiseRNA_Modeler::operator=( StepWiseRNA_Modeler const & src )
{
	job_parameters_ = src.job_parameters_;
	native_pose_ = src.native_pose_;
	moving_res_list_ = src.moving_res_list_;
	fixed_res_ = src.fixed_res_;
	minimize_res_ = src.minimize_res_;
	scorefxn_ = src.scorefxn_;
	silent_file_ = src.silent_file_;
	sampler_num_pose_kept_ = src.sampler_num_pose_kept_;
	num_pose_minimize_ = src.num_pose_minimize_;
	num_sampled_ = src.num_sampled_;
	sampler_native_screen_rmsd_cutoff_ = src.sampler_native_screen_rmsd_cutoff_;
	cluster_rmsd_ = src.cluster_rmsd_;
	native_edensity_score_cutoff_ = src.native_edensity_score_cutoff_;
	sampler_native_rmsd_screen_ = src.sampler_native_rmsd_screen_;
	o2prime_screen_ = src.o2prime_screen_;
	verbose_ = src.verbose_;
	distinguish_pucker_ = src.distinguish_pucker_;
	finer_sampling_at_chain_closure_ = src.finer_sampling_at_chain_closure_;
	PBP_clustering_at_chain_closure_ = src.PBP_clustering_at_chain_closure_;
	allow_syn_pyrimidine_ = src.allow_syn_pyrimidine_;
	extra_chi_ = src.extra_chi_;
	use_phenix_geo_ = src.use_phenix_geo_;
	virtual_sugar_legacy_mode_ = src.virtual_sugar_legacy_mode_;
	kic_sampling_ = src.kic_sampling_;
	kic_sampling_if_relevant_ = src.kic_sampling_if_relevant_;
	centroid_screen_ = src.centroid_screen_;
	VDW_atr_rep_screen_ = src.VDW_atr_rep_screen_;
	force_centroid_interaction_ = src.force_centroid_interaction_;
	choose_random_ = src.choose_random_;
	num_random_samples_ = src.num_random_samples_;
	skip_sampling_ = src.skip_sampling_;
	perform_minimize_ = src.perform_minimize_;
	minimize_and_score_sugar_ = src.minimize_and_score_sugar_;
	minimize_and_score_native_pose_ = src.minimize_and_score_native_pose_;
	rm_virt_phosphate_ = src.rm_virt_phosphate_;
	VDW_rep_alignment_RMSD_CUTOFF_ = src.VDW_rep_alignment_RMSD_CUTOFF_;
	output_pdb_ = src.output_pdb_;
	output_minimized_pose_list_ = src.output_minimized_pose_list_;
	VDW_rep_screen_physical_pose_clash_dist_cutoff_ = src.VDW_rep_screen_physical_pose_clash_dist_cutoff_;
	integration_test_mode_ = src.integration_test_mode_;
	allow_bulge_at_chainbreak_ = src.allow_bulge_at_chainbreak_;
	parin_favorite_output_ = src.parin_favorite_output_;
	reinitialize_CCD_torsions_ = src.reinitialize_CCD_torsions_;
	sampler_extra_epsilon_rotamer_ = src.sampler_extra_epsilon_rotamer_;
	sampler_extra_beta_rotamer_ = src.sampler_extra_beta_rotamer_;
	sample_both_sugar_base_rotamer_ = src.sample_both_sugar_base_rotamer_;
	sampler_include_torsion_value_in_tag_ = src.sampler_include_torsion_value_in_tag_;
	combine_long_loop_mode_ = src.combine_long_loop_mode_;
	do_not_sample_multiple_virtual_sugar_ = src.do_not_sample_multiple_virtual_sugar_;
	sample_ONLY_multiple_virtual_sugar_ = src.sample_ONLY_multiple_virtual_sugar_;
	sampler_assert_no_virt_sugar_sampling_ = src.sampler_assert_no_virt_sugar_sampling_;
	sampler_try_sugar_instantiation_ = src.sampler_try_sugar_instantiation_;
	allow_base_pair_only_centroid_screen_ = src.allow_base_pair_only_centroid_screen_;
	minimizer_perform_o2prime_pack_ = src.minimizer_perform_o2prime_pack_;
	minimizer_output_before_o2prime_pack_ = src.minimizer_output_before_o2prime_pack_;
	minimizer_rename_tag_ = src.minimizer_rename_tag_;
	stepwise_rna_minimizer_ = src.stepwise_rna_minimizer_;
	minimize_move_map_ = src.minimize_move_map_;
	minimizer_vary_bond_geometry_frequency_ = src.minimizer_vary_bond_geometry_frequency_;
	minimizer_allow_variable_bond_geometry_ = src.minimizer_allow_variable_bond_geometry_;
	minimizer_extra_minimize_res_ = src.minimizer_extra_minimize_res_;
	return *this;
}

/////////////////////
std::string
StepWiseRNA_Modeler::get_name() const {
	return "StepWiseRNA_Modeler";
}

/////////////////////
void
StepWiseRNA_Modeler::set_moving_res_and_reset( core::Size const moving_res ){
	moving_res_list_.clear();
	if ( moving_res > 0 ) moving_res_list_ = utility::tools::make_vector1( moving_res );
	job_parameters_ = 0; // Important: will trigger reset of job parameters when we get pose.
}

//////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Modeler::apply( core::pose::Pose & pose ){

	using namespace core::pose;
	using namespace core::chemical;
	using namespace core::kinematics;
	using namespace core::scoring;
	using namespace protocols::swa::rna;

	initialize_job_parameters( pose );

	utility::vector1< PoseOP > pose_list;
	do_residue_sampling( pose, pose_list );
	if ( sampling_successful( pose_list ) ) do_minimizing( pose, pose_list );

	job_parameters_ = 0; // Important: make sure that the next time this is used, job parameters is set explicitly -- or it will be reset.
}


//////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Modeler::initialize_job_parameters( pose::Pose & pose ){
	if ( job_parameters_ ) return;
	job_parameters_ = setup_job_parameters_for_swa_with_full_model_info( moving_res_list_, pose );
}

//////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Modeler::do_residue_sampling( pose::Pose & pose,
																					utility::vector1< PoseOP > & pose_list ){

	if ( ! skip_sampling_ && moving_res_list_.size() > 0 ) {
		// let's actually sample.
		StepWiseRNA_ResidueSampler stepwise_rna_residue_sampler( job_parameters_ );
		stepwise_rna_residue_sampler.set_silent_file ( silent_file_ + "_sampling" );
		stepwise_rna_residue_sampler.set_scorefxn ( scorefxn_ );
		stepwise_rna_residue_sampler.set_num_pose_kept ( sampler_num_pose_kept_ );
		stepwise_rna_residue_sampler.set_native_rmsd_screen ( sampler_native_rmsd_screen_ );
		stepwise_rna_residue_sampler.set_native_screen_rmsd_cutoff ( sampler_native_screen_rmsd_cutoff_ );
		stepwise_rna_residue_sampler.set_perform_o2prime_pack ( o2prime_screen_ );
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
		stepwise_rna_residue_sampler.set_reinitialize_CCD_torsions( reinitialize_CCD_torsions_ );
		stepwise_rna_residue_sampler.set_extra_epsilon_rotamer( sampler_extra_epsilon_rotamer_ );
		stepwise_rna_residue_sampler.set_extra_beta_rotamer( sampler_extra_beta_rotamer_ );
		stepwise_rna_residue_sampler.set_sample_both_sugar_base_rotamer( sample_both_sugar_base_rotamer_ ); //Nov 12, 2010
		stepwise_rna_residue_sampler.set_include_torsion_value_in_tag( sampler_include_torsion_value_in_tag_ );
		stepwise_rna_residue_sampler.set_combine_long_loop_mode( combine_long_loop_mode_ );
		stepwise_rna_residue_sampler.set_do_not_sample_multiple_virtual_sugar( do_not_sample_multiple_virtual_sugar_ );
		stepwise_rna_residue_sampler.set_sample_ONLY_multiple_virtual_sugar( sample_ONLY_multiple_virtual_sugar_ );
		stepwise_rna_residue_sampler.set_assert_no_virt_sugar_sampling( sampler_assert_no_virt_sugar_sampling_ );
		stepwise_rna_residue_sampler.set_try_sugar_instantiation( sampler_try_sugar_instantiation_ );

		base_centroid_screener_ = new screener::StepWiseRNA_BaseCentroidScreener ( pose, job_parameters_ );
		base_centroid_screener_->set_floating_base( job_parameters_->floating_base() );
		base_centroid_screener_->set_allow_base_pair_only_screen( allow_base_pair_only_centroid_screen_ );
		stepwise_rna_residue_sampler.set_base_centroid_screener( base_centroid_screener_ );

		user_input_VDW_bin_screener_ = new screener::StepWiseRNA_VDW_BinScreener();
		if ( VDW_rep_screen_info_.size() > 0 ) {
			user_input_VDW_bin_screener_->set_VDW_rep_alignment_RMSD_CUTOFF ( VDW_rep_alignment_RMSD_CUTOFF_ );
			user_input_VDW_bin_screener_->set_VDW_rep_delete_matching_res( VDW_rep_delete_matching_res_ );
			user_input_VDW_bin_screener_->set_physical_pose_clash_dist_cutoff( VDW_rep_screen_physical_pose_clash_dist_cutoff_ );
			user_input_VDW_bin_screener_->setup_using_user_input_VDW_pose( VDW_rep_screen_info_, pose, job_parameters_ );
			user_input_VDW_bin_screener_->set_output_pdb( output_pdb_ );
		}
		stepwise_rna_residue_sampler.set_user_input_VDW_bin_screener ( user_input_VDW_bin_screener_ );

		if ( choose_random_ ){
			stepwise_rna_residue_sampler.set_choose_random( true );
			stepwise_rna_residue_sampler.set_num_random_samples( num_random_samples_ );
			stepwise_rna_residue_sampler.set_cluster_rmsd( 0.0 ); // don't cluster.
		}
		print_JobParameters_info( job_parameters_, "job_parameters_COP", TR.Debug );

		stepwise_rna_residue_sampler.apply( pose );

		pose_list = stepwise_rna_residue_sampler.get_pose_list();
		if ( verbose_ ) stepwise_rna_residue_sampler.output_pose_list( silent_file_ + "_final_sample" );

	} else {
		add_to_pose_list( pose_list, pose, "input_pose" );
	}
}

//////////////////////////////////////////////////////////////////////////////
bool
StepWiseRNA_Modeler::sampling_successful( utility::vector1< PoseOP > & pose_list ){
	num_sampled_ = pose_list.size();
	if ( num_sampled_ == 0 ){
		TR << "WARNING! WARNING! WARNING! pose_list_.size() == 0! " << std::endl;
		if ( !output_minimized_pose_list_ ) return false; // don't do a minimize...
	}
	return true;
}

//////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Modeler::do_minimizing( pose::Pose & pose, utility::vector1< PoseOP > & pose_list ){

	if ( minimize_and_score_native_pose_ ) {
		runtime_assert ( get_native_pose() );
		add_to_pose_list( pose_list, *get_native_pose(), "working_native_pose" );
	}
	////////////////////////////////////////////////////////////////
	stepwise_rna_minimizer_ = new StepWiseRNA_Minimizer( pose_list, job_parameters_ );
	stepwise_rna_minimizer_->set_silent_file ( silent_file_ );
	stepwise_rna_minimizer_->set_verbose (  verbose_ );
	stepwise_rna_minimizer_->set_scorefxn ( scorefxn_ );
	stepwise_rna_minimizer_->set_base_centroid_screener ( base_centroid_screener_ );
	stepwise_rna_minimizer_->set_centroid_screen ( ( base_centroid_screener_ != 0 ) );
	stepwise_rna_minimizer_->set_perform_minimize(  perform_minimize_ );
	stepwise_rna_minimizer_->set_native_rmsd_screen (  sampler_native_rmsd_screen_ );
	stepwise_rna_minimizer_->set_native_edensity_score_cutoff ( native_edensity_score_cutoff_ );
	stepwise_rna_minimizer_->set_rm_virt_phosphate (  rm_virt_phosphate_ );
	stepwise_rna_minimizer_->set_native_screen_rmsd_cutoff (  sampler_native_screen_rmsd_cutoff_ + 1 ); //+1 for leniency Sept 20, 2010

	if ( integration_test_mode_ ) num_pose_minimize_ = 1;
	if ( num_pose_minimize_ > 0 ) stepwise_rna_minimizer_->set_num_pose_minimize ( num_pose_minimize_ );
	stepwise_rna_minimizer_->set_minimize_and_score_sugar ( minimize_and_score_sugar_ );
	stepwise_rna_minimizer_->set_user_input_VDW_bin_screener ( user_input_VDW_bin_screener_ );
	stepwise_rna_minimizer_->set_output_minimized_pose_list( output_minimized_pose_list_ );
	
	if ( minimize_move_map_ ) {
		stepwise_rna_minimizer_->set_move_map_list( make_vector1( *minimize_move_map_ ) );
		stepwise_rna_minimizer_->set_allow_insert( allow_insert_ );
	}
	
	stepwise_rna_minimizer_->set_perform_o2prime_pack(  minimizer_perform_o2prime_pack_ );
	stepwise_rna_minimizer_->set_output_before_o2prime_pack( minimizer_output_before_o2prime_pack_ );
	stepwise_rna_minimizer_->set_rename_tag( minimizer_rename_tag_ );
	stepwise_rna_minimizer_->set_extra_minimize_res( minimizer_extra_minimize_res_ );
	stepwise_rna_minimizer_->set_allow_variable_bond_geometry( minimizer_allow_variable_bond_geometry_ );
	stepwise_rna_minimizer_->set_vary_bond_geometry_frequency( minimizer_vary_bond_geometry_frequency_ );

	stepwise_rna_minimizer_->apply ( pose );

}

// Wrapper for routine below (which requires a const pose). If you call this,
// note that the pose's full_model_info object will be initialized based on its
// current fold_tree, cutpoint_variants, and any chain/residue-numbering information in
// PDBInfo.
StepWiseRNA_JobParametersOP
StepWiseRNA_Modeler::setup_job_parameters_for_swa_with_full_model_info( utility::vector1< Size > moving_res,
																																				core::pose::Pose & pose ){
	pose::full_model_info::make_sure_full_model_info_is_setup( pose );
	return setup_job_parameters_for_swa( moving_res_list_, pose );
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// This could go into a Util.hh if it ends up being more useful.
// Briefly, we need to make a StepWiseRNA_JobParameters object that will
// be fed into various StepWiseRNA movers. Setting this up can be quite complicated...
// it has become a grab bag of residue lists referring to the global pose, the working pose,
// sequence mappings, "is_Prepend_map", etc.
//
StepWiseRNA_JobParametersOP
StepWiseRNA_Modeler::setup_job_parameters_for_swa( utility::vector1< Size > moving_res,
																									 core::pose::Pose const & pose ){

	using namespace core::pose;
	using namespace core::chemical;
	using namespace core::id;
	using namespace protocols::swa;
	using namespace core::pose::full_model_info;

	runtime_assert( moving_res.size() <= 1 );

	// FullModelInfo can be setup with the command make_sure_full_model_info_is_setup( pose ) before calling this function.
	// It is used to figure out which suites need to be minimized (using info in domain map), and
	// is also necessary if there are jumps to floating bases -- need to know how many intervening residues are skipped.
	FullModelInfo const & full_model_info = const_full_model_info( pose );

	StepWiseRNA_JobParametersOP job_parameters = new StepWiseRNA_JobParameters;

	// not actually sure if we need this to be filled... should be taken care of by fixed_domain testing below.
	utility::vector1< Size > suites_that_must_be_minimized;
	std::string full_sequence = pose.sequence();

	// what if there is a virtual residue? need to remove it, actually, before running stepwise_rna_job_parameters_setup.
	Size nres = pose.total_residue();
	Size const rebuild_res =  ( moving_res.size() == 1 ) ? moving_res[1] : 0;
	utility::vector1< Size > not_rebuild_res;
	for ( Size n = 1; n <= nres; n++ ) if ( n != rebuild_res ) not_rebuild_res.push_back( n );
	utility::vector1< Size > fixed_res_guess = not_rebuild_res; // may be revised below.
	kic_sampling_ = false;

	if ( moving_res.size() == 1 ){

		TR.Debug << pose.fold_tree() << std::endl;
		TR.Debug << "Rebuild residue: " << rebuild_res << std::endl;

		utility::vector1< Size > input_res1, input_res2 /*blank*/, cutpoint_open, cutpoint_closed;
		Size cutpoint_closed_distal( 0 );
		input_res1 = not_rebuild_res;

		Size rebuild_suite( 0 );
		bool floating_base( false );
		bool const cut_at_previous =  (rebuild_res == 1) || pose.fold_tree().is_cutpoint( rebuild_res - 1 );
		bool const cut_at_next     =  (rebuild_res == pose.total_residue() ) || pose.fold_tree().is_cutpoint( rebuild_res );
		if ( cut_at_next && !cut_at_previous ){
			rebuild_suite = rebuild_res - 1;
		} else if ( cut_at_previous && !cut_at_next ){
			rebuild_suite = rebuild_res;
		} else if ( !cut_at_previous && !cut_at_next ){ // internal
			rebuild_suite = rebuild_res;
		} else {
			floating_base = true;
		}

		// for internal moves, need to be smart about input_res definitions -- 'domains' that are separated by moving residue.
		// this loop will also determine any chainbreak that requires closure
		utility::vector1< bool > partition_definition;
		Size floating_base_anchor_res( 0 );
		if ( floating_base ){
			floating_base_anchor_res = get_anchor_res( rebuild_res, pose );
			//if ( floating_base_anchor_res == rebuild_res + 2 ) floating_base_bulge_res = rebuild_res + 1;
			//			if ( floating_base_anchor_res == rebuild_res - 2 ) floating_base_bulge_res = rebuild_res - 1;
			//TR << TR.Red << "Floating_base bulge_res " << floating_base_bulge_res << TR.Reset << std::endl;
			partition_definition = get_partition_definition_floating_base( pose, rebuild_res );
		} else {
			partition_definition = get_partition_definition( pose, rebuild_suite );
		}

		utility::vector1< Size > const chains = figure_out_chains_from_full_model_info_const( pose );
		bool found_moving_cutpoint( false );
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
				if ( !floating_base ) runtime_assert( !found_moving_cutpoint );
				found_moving_cutpoint = true;
				if ( pose.residue_type( n ).has_variant_type( CUTPOINT_LOWER ) ){
					runtime_assert( pose.residue_type( n + 1 ).has_variant_type( CUTPOINT_UPPER  ) );
					cutpoint_closed.push_back( n );
					if ( !floating_base ||
							 ( floating_base_anchor_res > rebuild_res &&  n < rebuild_res ) ||
							 ( floating_base_anchor_res < rebuild_res &&  n > rebuild_res ) ){
						runtime_assert( cutpoint_closed_distal == 0 ); // can only handle one such distal cutpoint, at present.
						cutpoint_closed_distal = n;
					}
				} else if ( chains[ n ] != chains[ n+1 ]  ){
					cutpoint_open.push_back( n );
				}
			}
		}

		TR.Debug << "INPUT_RES1 " << make_tag_with_dashes(input_res1) << " INPUT_RES2 " << make_tag_with_dashes(input_res2) << " REBUILD_RES " << rebuild_res << " REBUILD_SUITE " <<  rebuild_suite << " CUTPOINT_CLOSED " << cutpoint_closed << std::endl;

		// if there's really just a single nucleotide being rebuilt, its nucleoside is not fixed.
		// but if there's a whole chunk of stuff, its sugar & base are assumed fixed.
		//note that fixed_res_guess, which is really a list of fixed nucleosides,
		// should now include the 'moving res' [unless re-specified by user down below].

		// check if this is an 'internal' move.
		if ( input_res1.size() > 1 && input_res2.size() > 1 ) fixed_res_guess.push_back( rebuild_res );

		// To specify that the suite moves, we actually need to directly address the movemap... see below.
		if ( rebuild_suite > 0 ) suites_that_must_be_minimized.push_back( rebuild_suite );
		if ( cutpoint_closed.size() > 0 ) {
			if ( kic_sampling_if_relevant_ ) kic_sampling_ = true;
			for ( Size i = 1; i <= cutpoint_closed.size(); i++ ) {
				suites_that_must_be_minimized.push_back( cutpoint_closed[i] );
				if ( !pose.fold_tree().is_cutpoint( cutpoint_closed[i] ) ) utility_exit_with_message( "StepWiseRNA requires a chainbreak right at sampled residue" );
			}
			if ( rebuild_res == 1 || rebuild_res == pose.total_residue() ) utility_exit_with_message( "StepWiseRNA requires that residue is not at terminus!");
		}

		std::string full_sequence = full_model_info.full_sequence();
		if ( full_sequence[nres - 1] == 'X' ){
			full_sequence = full_sequence.substr( 0, nres - 1 );
			nres -= 1;
		}
		StepWiseRNA_JobParametersSetup stepwise_rna_job_parameters_setup( full_model_info.sub_to_full( moving_res ),
																																			full_sequence,
																																			full_model_info.sub_to_full( input_res1 ),
																																			full_model_info.sub_to_full( input_res2 ),
																																			full_model_info.sub_to_full( cutpoint_open ),
																																			full_model_info.sub_to_full( cutpoint_closed_distal ) );
		// following still might be worth doing -- just didn't seem necessary at this point.
		//		stepwise_rna_job_parameters_setup.set_cutpoint_closed_list( full_model_info.sub_to_full( cutpoint_closed ) );
		stepwise_rna_job_parameters_setup.set_fixed_res( full_model_info.sub_to_full( fixed_res_guess ) );
		stepwise_rna_job_parameters_setup.set_floating_base( floating_base );
		if ( floating_base ){
			stepwise_rna_job_parameters_setup.set_assert_jump_point_in_fixed_res( false );
			stepwise_rna_job_parameters_setup.set_floating_base_anchor_res( full_model_info.sub_to_full( floating_base_anchor_res ) );
		}
		if ( rmsd_res_list_.size() > 0 ) stepwise_rna_job_parameters_setup.set_rmsd_res_list( rmsd_res_list_ /*global numbering*/ );
		else  stepwise_rna_job_parameters_setup.set_rmsd_res_list( full_model_info.sub_to_full( make_vector1( rebuild_res ) ) );

		if ( rebuild_res > 1 &&
				 pose.residue( rebuild_res - 1 ).has_variant_type( "VIRTUAL_RIBOSE" ) &&
				 ( !pose.fold_tree().is_cutpoint( rebuild_res - 1 ) ||
					 is_cutpoint_closed( pose, rebuild_res - 1 ) ) &&
				 !pose.fold_tree().is_cutpoint( rebuild_res ) &&
				 pose.fold_tree().jump_nr( rebuild_res - 1, rebuild_res + 1) > 0 ) stepwise_rna_job_parameters_setup.set_rebuild_bulge_mode( true );
		if ( rebuild_res < pose.total_residue() &&
				 pose.residue( rebuild_res + 1 ).has_variant_type( "VIRTUAL_RIBOSE" ) &&
				 ( !pose.fold_tree().is_cutpoint( rebuild_res ) ||
					 is_cutpoint_closed( pose, rebuild_res ) ) &&
				 !pose.fold_tree().is_cutpoint( rebuild_res - 1 ) &&
				 pose.fold_tree().jump_nr( rebuild_res - 1, rebuild_res + 1 ) > 0 ) stepwise_rna_job_parameters_setup.set_rebuild_bulge_mode( true );

		utility::vector1< std::string > jump_point_pair_list;
		core::kinematics::FoldTree const & f = pose.fold_tree();
		for ( Size n = 1; n <= f.num_jump(); n++ ){
			jump_point_pair_list.push_back( ObjexxFCL::string_of( full_model_info.sub_to_full( f.upstream_jump_residue( n ) ) ) + "-" +
																			ObjexxFCL::string_of( full_model_info.sub_to_full( f.downstream_jump_residue( n ) ) ) );
		}
		stepwise_rna_job_parameters_setup.set_jump_point_pair_list( jump_point_pair_list ); //Important!: Needs to be called after set_fixed_res

		utility::vector1< std::string > alignment_res; //why is this a string vector?????
		for ( Size n = 1; n <= fixed_res_guess.size(); n++ ) alignment_res.push_back( ObjexxFCL::string_of( full_model_info.sub_to_full( fixed_res_guess[ n ] ) ) );
		stepwise_rna_job_parameters_setup.set_alignment_res( alignment_res );
		stepwise_rna_job_parameters_setup.set_native_alignment_res( full_model_info.sub_to_full( fixed_res_guess ) );

		// could use this later to minimize more residues...
		//stepwise_rna_job_parameters_setup.set_global_sample_res_list( option[ global_sample_res_list ]() ); //March 20, 2011

		// NOT SURE ABOUT THIS. false by default, but shows up as true in 'normal' erraser runs.
		stepwise_rna_job_parameters_setup.set_allow_chain_boundary_jump_partner_right_at_fixed_BP ( true );

		// NOT SURE ABOUT THIS...
		stepwise_rna_job_parameters_setup.set_add_virt_res_as_root( true );

		// ignore fold tree setup which is hopelessly complicated in stepwise_rna_job_parameters_setup.
		stepwise_rna_job_parameters_setup.force_fold_tree( pose.fold_tree() );
		stepwise_rna_job_parameters_setup.apply();
		job_parameters = stepwise_rna_job_parameters_setup.job_parameters();
		job_parameters->set_working_native_pose( get_native_pose() );
		TR.Debug << "past job_parameters initialization " << std::endl;
	}

	// If setup_job_parameters_for_swa() is called, then the user has not supplied their own StepWiseRNA_JobParameters object to the modeler. This means that StepWiseRNAMinimizer will not be generating a move map on its own, so we need to supply it with an AllowInsert object in order to handle the possibility of variable bond geometries.
	allow_insert_ = new toolbox::AllowInsert(pose); // Default constructor that allows everything to move
	
	// user input minimize_res...
	if ( minimize_res_.size() > 0 ) { // specifying more residues which could move during the minimize step -- standard for SWM.
		fixed_res_.clear();
		for ( Size n = 1; n <= nres; n++ ) {
			if ( !minimize_res_.has_value( n ) )	{
				fixed_res_.push_back( n );
				allow_insert_->set( n, false );
			}
		}
	} else if ( fixed_res_.size() > 0 ){ // how 'standard' SWA specifies moving information.
		runtime_assert( minimize_res_.size() == 0 );
		for ( Size n = 1; n <= nres; n++ ) {
			if ( !fixed_res_.has_value( n ) )	{
				minimize_res_.push_back( n );
			} else {
				allow_insert_->set( n, false );
			}
		}
	} else { // 'reasonable' default behavior, inferred above.
		for ( Size n = 1; n <= nres; n++ ) {
			if ( !fixed_res_guess.has_value( n ) )	{
				minimize_res_.push_back( n );
			} else {
				allow_insert_->set( n, false );
			}
		}
		fixed_res_ = fixed_res_guess;
	}

	if ( fixed_res_.size() > 0 ) 	job_parameters->set_working_fixed_res( fixed_res_ ); // is this necessary if we just supply movemap?
	job_parameters->set_working_native_alignment( fixed_res_ );
	
	//Now we perform the additional task of updating the AllowInsert object based on any optional additional residues that the user wants minimized, as specified in minimizer_extra_minimize_res. The intended mode of operation is that the user decides on extra_minimize_res at the high level for an entire SWM run, and then changes minimize_res_ to control a specific minimization event.
	update_allow_insert_with_extra_minimize_res( pose, allow_insert_, minimizer_extra_minimize_res_ );

	minimize_move_map_ = new core::kinematics::MoveMap;
	//figure_out_swa_rna_movemap( *minimize_move_map_, pose, minimize_res_ );
	figure_out_swa_rna_movemap( *minimize_move_map_, pose, allow_insert_ );

	// last, but not least, there might be some information in the domain map. Note
	// that generally we could instead replace fixed_res with an inputted domain map.
	// that is, get rid of fixed_res_ & minimize_res_ and instead have a local fixed_domain_map,
	// which can instead be updated by set_fixed_res.
	// how to tell Modeler to *not* minimize additional suites?
	utility::vector1< Size > const & res_list = full_model_info.res_list();
	utility::vector1< Size > const & fixed_domain_map = full_model_info.fixed_domain_map();
	utility::vector1< Size > const & cutpoint_open_in_full_model = full_model_info.cutpoint_open_in_full_model();
	for ( Size n = 1; n < pose.total_residue(); n++ ){
		if ( !cutpoint_open_in_full_model.has_value( res_list[ n ] ) &&
				( res_list[ n + 1 ]  == res_list[ n ] + 1 ) &&
				 !minimize_res_.has_value( n )  && !minimize_res_.has_value( n+1 )  &&
				 ( fixed_domain_map[ res_list[ n + 1 ] ] !=  fixed_domain_map[ res_list[ n ] ] ) &&
				 !suites_that_must_be_minimized.has_value( n ) ){
			TR.Debug << "ADDING NEW SUITE TO BE MINIMIZED BASED ON LOCATION AT DOMAIN BOUNDARY: " << n << std::endl;
			suites_that_must_be_minimized.push_back( n );
		}
	}
	TR.Debug << "SUITES_THAT_MUST_BE_MINIMIZED " << suites_that_must_be_minimized << std::endl;

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
	if ( !stepwise_rna_minimizer_ ) return;
	stepwise_rna_minimizer_->output_pose_wrapper( out_tag, pose, out_silent_file );
}

///////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Modeler::add_to_pose_list( utility::vector1< core::pose::PoseOP > & pose_list, pose::Pose const & pose, std::string const pose_tag ) const {
	core::pose::PoseOP pose_op = pose.clone();
	tag_into_pose( *pose_op, pose_tag );
	pose_list.push_back( pose_op );
}


} //rna
} //swa
} //protocols
