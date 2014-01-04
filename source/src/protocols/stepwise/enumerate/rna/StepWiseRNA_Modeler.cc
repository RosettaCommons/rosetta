// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/enumerate/rna/SWA_Modeler.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu

#include <protocols/stepwise/enumerate/rna/StepWiseRNA_Modeler.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_ModelerOptions.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_ResidueSampler.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_JobParametersSetup.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_Minimizer.hh>
#include <protocols/stepwise/enumerate/rna/screener/StepWiseRNA_BaseCentroidScreener.hh>
#include <protocols/stepwise/enumerate/rna/screener/StepWiseRNA_VDW_BinScreener.hh>
#include <protocols/stepwise/StepWiseUtil.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_Util.hh>

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

static basic::Tracer TR( "protocols.stepwise.rna.StepWiseRNA_Modeler" );

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
namespace stepwise {
namespace enumerate {
namespace rna {

//Constructor
StepWiseRNA_Modeler::StepWiseRNA_Modeler( core::scoring::ScoreFunctionOP scorefxn ) :
	scorefxn_( scorefxn ),
	options_( new StepWiseRNA_ModelerOptions ),
	num_sampled_( 0 )
{
}


StepWiseRNA_Modeler::StepWiseRNA_Modeler( core::Size const sample_res,
																	core::scoring::ScoreFunctionOP scorefxn ) :
	moving_res_list_( make_vector1( sample_res ) ),
	scorefxn_( scorefxn ),
	options_( new StepWiseRNA_ModelerOptions ),
	num_sampled_( 0 )
{
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
	num_sampled_ = src.num_sampled_;
	native_pose_ = src.native_pose_;
	moving_res_list_ = src.moving_res_list_;
	fixed_res_ = src.fixed_res_;
	minimize_res_ = src.minimize_res_;
	terminal_res_ = src.terminal_res_;
	scorefxn_ = src.scorefxn_;
	options_ = src.options_;
	stepwise_rna_minimizer_ = src.stepwise_rna_minimizer_;
	minimize_move_map_ = src.minimize_move_map_;
	minimizer_extra_minimize_res_ = src.minimizer_extra_minimize_res_;
	syn_chi_res_list_ = src.syn_chi_res_list_;

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
	using namespace protocols::stepwise::enumerate::rna;

	initialize_job_parameters_and_root( pose );

	utility::vector1< PoseOP > pose_list;
	do_residue_sampling( pose, pose_list );

	if ( sampling_successful( pose_list ) ) do_minimizing( pose, pose_list );

	job_parameters_ = 0; // Important: make sure that the next time this is used, job parameters is set explicitly -- or it will be reset.
}


//////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Modeler::initialize_job_parameters_and_root( pose::Pose & pose ){
	if ( job_parameters_ ) return;
	job_parameters_ = setup_job_parameters_for_stepwise_with_full_model_info( moving_res_list_, pose );
	setup_root_based_on_full_model_info( pose, job_parameters_ );
}

//////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Modeler::do_residue_sampling( pose::Pose & pose,
																					utility::vector1< PoseOP > & pose_list ){

	if ( options_->skip_sampling() || moving_res_list_.size() == 0 ) {
		add_to_pose_list( pose_list, pose, "input_pose" );
		return;
	}

	StepWiseRNA_ResidueSampler stepwise_rna_residue_sampler( job_parameters_ );
	stepwise_rna_residue_sampler.set_sampling_silent_file ( silent_file_ + "_sampling" );
	stepwise_rna_residue_sampler.set_scorefxn ( scorefxn_ );
	stepwise_rna_residue_sampler.set_options( options_ );

	base_centroid_screener_ = new screener::StepWiseRNA_BaseCentroidScreener ( pose, job_parameters_ );
	base_centroid_screener_->set_floating_base( job_parameters_->floating_base() );
	base_centroid_screener_->set_allow_base_pair_only_screen( options_->allow_base_pair_only_centroid_screen() );
	stepwise_rna_residue_sampler.set_base_centroid_screener( base_centroid_screener_ );

	user_input_VDW_bin_screener_ = new screener::StepWiseRNA_VDW_BinScreener();
	if ( VDW_rep_screen_info_.size() > 0 ) {
		options_->setup_options_for_VDW_bin_screener( user_input_VDW_bin_screener_ );
		user_input_VDW_bin_screener_->setup_using_user_input_VDW_pose( VDW_rep_screen_info_, pose, job_parameters_ );
	}
	stepwise_rna_residue_sampler.set_user_input_VDW_bin_screener ( user_input_VDW_bin_screener_ );

	print_JobParameters_info( job_parameters_, "job_parameters_COP", TR.Debug );

	stepwise_rna_residue_sampler.apply( pose );

	pose_list = stepwise_rna_residue_sampler.get_pose_list();
	if ( options_->verbose() ) stepwise_rna_residue_sampler.output_pose_list( silent_file_ + "_final_sample" );

}

//////////////////////////////////////////////////////////////////////////////
bool
StepWiseRNA_Modeler::sampling_successful( utility::vector1< PoseOP > & pose_list ){
	num_sampled_ = pose_list.size();
	if ( num_sampled_ == 0 ){
		TR << "WARNING! WARNING! WARNING! pose_list_.size() == 0! " << std::endl;
		if ( !options_->output_minimized_pose_list() ) return false; // don't do a minimize...
	}
	return true;
}

//////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Modeler::do_minimizing( pose::Pose & pose, utility::vector1< PoseOP > & pose_list ){

	if ( options_->minimize_and_score_native_pose() ) {
		runtime_assert ( get_native_pose() );
		add_to_pose_list( pose_list, *get_native_pose(), "working_native_pose" );
	}

	stepwise_rna_minimizer_ = new StepWiseRNA_Minimizer( pose_list, job_parameters_ );
	stepwise_rna_minimizer_->set_scorefxn( scorefxn_ );
	stepwise_rna_minimizer_->set_silent_file( silent_file_ );
	stepwise_rna_minimizer_->set_base_centroid_screener( base_centroid_screener_ );
	stepwise_rna_minimizer_->set_user_input_VDW_bin_screener( user_input_VDW_bin_screener_ );

	StepWiseRNA_ModelerOptionsOP minimizer_options = options_->clone();
	minimizer_options->set_sampler_native_screen_rmsd_cutoff( options_->sampler_native_screen_rmsd_cutoff() + 1.0 );
	stepwise_rna_minimizer_->set_options( minimizer_options );

	if ( minimize_move_map_ ) {
		stepwise_rna_minimizer_->set_move_map_list( make_vector1( *minimize_move_map_ ) );
		stepwise_rna_minimizer_->set_allow_insert( allow_insert_ );
	}
	stepwise_rna_minimizer_->set_extra_minimize_res( minimizer_extra_minimize_res_ );

	stepwise_rna_minimizer_->apply ( pose );

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Wrapper for routine below (which requires a const pose). If you call this,
// note that the pose's full_model_info object will be initialized based on its
// current fold_tree, cutpoint_variants, and any chain/residue-numbering information in
// PDBInfo.
StepWiseRNA_JobParametersOP
StepWiseRNA_Modeler::setup_job_parameters_for_stepwise_with_full_model_info( utility::vector1< Size > moving_res,
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
//    -- rhiju
//
// PS> This whole thing has got to be refactored. Its insane.
//
StepWiseRNA_JobParametersOP
StepWiseRNA_Modeler::setup_job_parameters_for_swa( utility::vector1< Size > moving_res,
																									 core::pose::Pose const & pose ){

	using namespace core::pose;
	using namespace core::chemical;
	using namespace core::id;
	using namespace protocols::stepwise;
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

	if ( moving_res.size() == 1 ){

		TR.Debug << pose.fold_tree() << std::endl;
		TR.Debug << "Rebuild residue: " << rebuild_res << std::endl;

		utility::vector1< Size > input_res1, input_res2 /*blank*/, cutpoint_open, cutpoint_closed;
		Size cutpoint_closed_distal( 0 );
		input_res1 = not_rebuild_res;

		Size rebuild_suite( 0 );
		bool floating_base( false ), is_internal( false );
		bool cut_at_previous =  (rebuild_res == 1) || pose.fold_tree().is_cutpoint( rebuild_res - 1 );
		bool cut_at_next     =  (rebuild_res == pose.total_residue() ) || pose.fold_tree().is_cutpoint( rebuild_res );
		if ( cut_at_next && cut_at_previous ){
			floating_base = true;
		} else {
			cut_at_previous = cut_at_previous && !check_jump_to_previous_residue_in_chain( pose, rebuild_res, moving_res );
			cut_at_next     = cut_at_next     && !check_jump_to_next_residue_in_chain( pose, rebuild_res, moving_res );
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
							 ( floating_base_anchor_res > rebuild_res &&  rebuild_res >  n ) ||
							 ( floating_base_anchor_res < rebuild_res &&  rebuild_res <= n ) ){
						runtime_assert( cutpoint_closed_distal == 0 ); // can only handle one such distal cutpoint, at present.
						cutpoint_closed_distal = n;
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
			runtime_assert( is_internal );
		}

		// To specify that the suite moves, we actually need to directly address the movemap... see below.
		if ( rebuild_suite > 0 ) suites_that_must_be_minimized.push_back( rebuild_suite );
		if ( cutpoint_closed.size() > 0 ) {
			for ( Size i = 1; i <= cutpoint_closed.size(); i++ ) {
				suites_that_must_be_minimized.push_back( cutpoint_closed[i] );
				if ( !pose.fold_tree().is_cutpoint( cutpoint_closed[i] ) ) utility_exit_with_message( "StepWiseRNA requires a chainbreak right at sampled residue" );
			}
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
		if ( is_internal ) stepwise_rna_job_parameters_setup.set_force_internal( true );
		if ( floating_base ){
			stepwise_rna_job_parameters_setup.set_assert_jump_point_in_fixed_res( false );
			stepwise_rna_job_parameters_setup.set_floating_base_anchor_res( full_model_info.sub_to_full( floating_base_anchor_res ) );
		}
		if ( rmsd_res_list_.size() > 0 ) stepwise_rna_job_parameters_setup.set_rmsd_res_list( rmsd_res_list_ /*global numbering*/ );
		else  stepwise_rna_job_parameters_setup.set_rmsd_res_list( full_model_info.sub_to_full( make_vector1( rebuild_res ) ) );

		core::kinematics::FoldTree const & f = pose.fold_tree();

		if ( rebuild_res > 1 &&
				 pose.residue( rebuild_res - 1 ).has_variant_type( "VIRTUAL_RIBOSE" ) &&
				 ( !f.is_cutpoint( rebuild_res - 1 ) ||
					 is_cutpoint_closed( pose, rebuild_res - 1 ) ) &&
				 !f.is_cutpoint( rebuild_res ) &&
				 f.jump_nr( rebuild_res - 1, rebuild_res + 1) > 0 ) stepwise_rna_job_parameters_setup.set_rebuild_bulge_mode( true );
		if ( rebuild_res < pose.total_residue() &&
				 pose.residue( rebuild_res + 1 ).has_variant_type( "VIRTUAL_RIBOSE" ) &&
				 ( !f.is_cutpoint( rebuild_res ) ||
					 is_cutpoint_closed( pose, rebuild_res ) ) &&
				 !f.is_cutpoint( rebuild_res - 1 ) &&
				 f.jump_nr( rebuild_res - 1, rebuild_res + 1 ) > 0 ) stepwise_rna_job_parameters_setup.set_rebuild_bulge_mode( true );

		if ( !floating_base && !is_internal &&
				 ( rebuild_suite == 1 ||
					 ( f.is_cutpoint( rebuild_suite - 1) && !is_cutpoint_closed( pose, rebuild_suite - 1 ) ) ) &&
				 ( rebuild_suite+1 == pose.total_residue() ||
					 ( f.is_cutpoint( rebuild_suite + 1) && !is_cutpoint_closed( pose, rebuild_suite + 1 ) ) ) &&
				 !pose.residue( rebuild_suite ).has_variant_type( "VIRTUAL_RIBOSE" ) &&
				 !pose.residue( rebuild_suite + 1 ).has_variant_type( "VIRTUAL_RIBOSE" ) ) {
			TR << "SAMPLE_BOTH_SUGAR_BASE_ROTAMER!" << std::endl;
			stepwise_rna_job_parameters_setup.set_sample_both_sugar_base_rotamer( true );
		}

		utility::vector1< std::string > jump_point_pair_list;
		for ( Size n = 1; n <= f.num_jump(); n++ ){
			jump_point_pair_list.push_back( ObjexxFCL::string_of( full_model_info.sub_to_full( f.upstream_jump_residue( n ) ) ) + "-" +
																			ObjexxFCL::string_of( full_model_info.sub_to_full( f.downstream_jump_residue( n ) ) ) );
		}
		stepwise_rna_job_parameters_setup.set_jump_point_pair_list( jump_point_pair_list ); //Important!: Needs to be called after set_fixed_res

		utility::vector1< std::string > alignment_res; //why is this a string vector?????
		for ( Size n = 1; n <= fixed_res_guess.size(); n++ ) alignment_res.push_back( ObjexxFCL::string_of( full_model_info.sub_to_full( fixed_res_guess[ n ] ) ) );
		stepwise_rna_job_parameters_setup.set_alignment_res( alignment_res );
		stepwise_rna_job_parameters_setup.set_native_alignment_res( full_model_info.sub_to_full( fixed_res_guess ) );
		stepwise_rna_job_parameters_setup.set_terminal_res( terminal_res_ );

		// could use this later to minimize more residues...
		//stepwise_rna_job_parameters_setup.set_global_sample_res_list( option[ global_sample_res_list ]() ); //March 20, 2011

		// NOT SURE ABOUT THIS. false by default, but shows up as true in 'normal' erraser runs.
		stepwise_rna_job_parameters_setup.set_allow_chain_boundary_jump_partner_right_at_fixed_BP ( true );

		// NOT SURE ABOUT THIS...
		stepwise_rna_job_parameters_setup.set_add_virt_res_as_root( true );

		// ignore fold tree setup which is hopelessly complicated in stepwise_rna_job_parameters_setup.
		stepwise_rna_job_parameters_setup.force_fold_tree( f );
		stepwise_rna_job_parameters_setup.apply();
		job_parameters = stepwise_rna_job_parameters_setup.job_parameters();
		job_parameters->set_working_native_pose( get_native_pose() );
		job_parameters->set_force_syn_chi_res_list( syn_chi_res_list_ ); // will get automatically converted to working numbering.

		if ( is_internal ) runtime_assert( job_parameters->is_internal() );
		TR.Debug << "past job_parameters initialization " << std::endl;
	}

	// If setup_job_parameters_for_stepwise_) is called, then the user has not supplied their own StepWiseRNA_JobParameters object to the modeler. This means that StepWiseRNAMinimizer will not be generating a move map on its own, so we need to supply it with an AllowInsert object in order to handle the possibility of variable bond geometries.
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
	//figure_out_stepwise_rna_movemap( *minimize_move_map_, pose, minimize_res_ );
	figure_out_stepwise_rna_movemap( *minimize_move_map_, pose, allow_insert_ );

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

///////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Modeler::setup_root_based_on_full_model_info( pose::Pose & pose, StepWiseRNA_JobParametersCOP & job_parameters ){

	using namespace core::kinematics;

	FullModelInfo const & full_model_info = const_full_model_info( pose );
	utility::vector1< Size > const & res_list = get_res_list_from_full_model_info( pose );
	utility::vector1< Size > const & cutpoint_open_in_full_model = full_model_info.cutpoint_open_in_full_model();
	std::string const full_sequence = full_model_info.full_sequence();

	utility::vector1< Size > root_partition_res;
	if ( job_parameters_->working_moving_res()  > 0 ){
		ObjexxFCL::FArray1D < bool > const & partition_definition = job_parameters->partition_definition();
		utility::vector1< Size > partition_res0, partition_res1;
		for ( Size i = 1; i <= pose.total_residue(); i++ ){
			if ( partition_definition( i ) ) partition_res0.push_back( i );
			else partition_res1.push_back( i );
		}
		if ( partition_res0.size() == partition_res1.size() ) return;

		bool const root_partition = ( partition_res0.size() > partition_res1.size() );
		root_partition_res = root_partition ? partition_res0 : partition_res1;
	} else {
		for ( Size i = 1; i <= pose.total_residue(); i++ ) root_partition_res.push_back( i );
	}

	Size new_root( 0 ), possible_root( 0 );
	for ( Size n = 1; n <= root_partition_res.size(); n++ ){
		Size const i = root_partition_res[ n ];
		if ( !pose.fold_tree().possible_root( i ) ) continue;
		if ( res_list[ i ] == 1 ||
				 cutpoint_open_in_full_model.has_value( res_list[ i ] - 1 ) ){ // great, nothing will ever get prepended here.
			new_root = i; break;
		}
		if ( res_list[ i ] == full_sequence.size() ||
				 cutpoint_open_in_full_model.has_value( res_list[ i ] ) ){ // great, nothing will ever get appended here.
			new_root = i; break;
		}
		if ( possible_root == 0) possible_root = i;
	}
	if ( new_root == 0 ){
		runtime_assert( possible_root > 0 );
		new_root = possible_root;
	}

	FoldTree f = pose.fold_tree();
	if ( new_root == f.root() ) return;
	f.reorder( new_root );
	pose.fold_tree( f );
}

/////////////////////////////////////////////////////
StepWiseRNA_ModelerOptionsCOP
StepWiseRNA_Modeler::options(){
	return options_;
}

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_Modeler::set_options( StepWiseRNA_ModelerOptionsCOP options ){
	options_ = options;
}

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_Modeler::set_skip_sampling( bool const & setting ){
	if( options_->skip_sampling() != setting ){
		// needed to get aroud constant status of options -- which is desirable in most contexts.
		enumerate::rna::StepWiseRNA_ModelerOptionsOP new_options = options_->clone();
		new_options->set_skip_sampling( setting );
		set_options( new_options );
	}
}


} //rna
} //enumerate
} //stepwise
} //protocols
