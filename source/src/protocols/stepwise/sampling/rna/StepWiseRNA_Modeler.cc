// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/rna/SWA_Modeler.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu

#include <protocols/stepwise/sampling/rna/StepWiseRNA_Modeler.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_ModelerOptions.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_ResidueSampler.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_JobParametersUtil.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_Minimizer.hh>
#include <protocols/stepwise/sampling/rna/checker/RNA_BaseCentroidChecker.hh>
#include <protocols/stepwise/sampling/rna/checker/RNA_VDW_BinChecker.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_Util.hh>
#include <protocols/stepwise/StepWiseUtil.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/util.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/stream_util.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.rna.StepWiseRNA_Modeler" );

using utility::tools::make_vector1;

////////////////////////////////////////////////////////////////////////////////
//
// This is meant to be a simple 'master interface' to StepWiseAssembly and
//  StepWiseMonteCarlo functions, requiring only a pose and the residue to be
//  sampled.
//
// All of the complexities of setting up the JobParameters, etc. are hidden
//  in a single wrapper function in StepWiseRNA_JobParametersUtil.
//
//  -- Rhiju
////////////////////////////////////////////////////////////////////////////////


namespace protocols {
namespace stepwise {
namespace sampling {
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
		scorefxn_( scorefxn ),
		moving_res_( sample_res ),
		options_( new StepWiseRNA_ModelerOptions ),
		num_sampled_( 0 )
{}

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
	moving_res_ = src.moving_res_;
	working_fixed_res_ = src.working_fixed_res_;
	working_minimize_res_ = src.working_minimize_res_;
	terminal_res_ = src.terminal_res_;
	scorefxn_ = src.scorefxn_;
	options_ = src.options_;
	stepwise_rna_minimizer_ = src.stepwise_rna_minimizer_;
	minimize_move_map_ = src.minimize_move_map_;
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
	moving_res_ = moving_res;
	job_parameters_ = 0; // Important: will trigger reset of job parameters when we get pose.
}

//////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Modeler::apply( core::pose::Pose & pose ){

	using namespace core::pose;

	initialize_job_parameters_and_root( pose );

	utility::vector1< PoseOP > pose_list;
	do_residue_sampling( pose, pose_list );

	if ( sampling_successful( pose_list ) ) do_minimizing( pose, pose_list );
	job_parameters_ = 0; // Important: make sure that the next time this is used, job parameters is set explicitly -- or it will be reset.
}


//////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Modeler::initialize_job_parameters_and_root( pose::Pose & pose ){
	if ( job_parameters_ != 0 ) return;
	revise_root_and_moving_res( pose, moving_res_ );
	job_parameters_ = setup_job_parameters_for_stepwise_with_full_model_info( pose );
}

//////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Modeler::do_residue_sampling( pose::Pose & pose,
																					utility::vector1< PoseOP > & pose_list ){

	if ( options_->skip_sampling() || moving_res_ == 0 ) {
		add_to_pose_list( pose_list, pose, "input_pose" );
		return;
	}

	StepWiseRNA_ResidueSampler stepwise_rna_residue_sampler( job_parameters_ );
	stepwise_rna_residue_sampler.set_sampling_silent_file ( silent_file_ + "_sampling" );
	stepwise_rna_residue_sampler.set_scorefxn ( scorefxn_ );
	stepwise_rna_residue_sampler.set_options( options_ );

	base_centroid_checker_ = new checker::RNA_BaseCentroidChecker ( pose, job_parameters_,
																																	options_->tether_jump());
	base_centroid_checker_->set_floating_base( job_parameters_->floating_base() &&
																						 job_parameters_->working_moving_partition_pos().size() == 1  );
	base_centroid_checker_->set_allow_base_pair_only_screen( options_->allow_base_pair_only_centroid_screen() );
	stepwise_rna_residue_sampler.set_base_centroid_checker( base_centroid_checker_ );

	user_input_VDW_bin_checker_ = new checker::RNA_VDW_BinChecker();
	if ( VDW_rep_screen_info_.size() > 0 ) {
		options_->setup_options_for_VDW_bin_checker( user_input_VDW_bin_checker_ );
		user_input_VDW_bin_checker_->setup_using_user_input_VDW_pose( VDW_rep_screen_info_, pose, job_parameters_ );
	}
	stepwise_rna_residue_sampler.set_user_input_VDW_bin_checker ( user_input_VDW_bin_checker_ );

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
	stepwise_rna_minimizer_->set_base_centroid_checker( base_centroid_checker_ );
	stepwise_rna_minimizer_->set_user_input_VDW_bin_checker( user_input_VDW_bin_checker_ );

	StepWiseRNA_ModelerOptionsOP minimizer_options = options_->clone();
	minimizer_options->set_sampler_native_screen_rmsd_cutoff( options_->sampler_native_screen_rmsd_cutoff() + 1.0 );
	stepwise_rna_minimizer_->set_options( minimizer_options );

	if ( minimize_move_map_ ) {
		stepwise_rna_minimizer_->set_move_map( minimize_move_map_ );
		stepwise_rna_minimizer_->set_allow_insert( allow_insert_ );
	}
	stepwise_rna_minimizer_->set_extra_minimize_res( const_full_model_info( pose ).extra_minimize_res() );

	stepwise_rna_minimizer_->apply ( pose );

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Wrapper for routine below (which requires a const pose). If you call this,
// note that the pose's full_model_info object will be initialized based on its
// current fold_tree, cutpoint_variants, and any chain/residue-numbering information in
// PDBInfo.
StepWiseRNA_JobParametersOP
StepWiseRNA_Modeler::setup_job_parameters_for_stepwise_with_full_model_info( core::pose::Pose & pose )
{
	pose::full_model_info::make_sure_full_model_info_is_setup( pose );
	StepWiseRNA_JobParametersOP job_parameters;
	job_parameters = setup_job_parameters_for_swa( moving_res_, pose,
																								 get_native_pose(),
																								 calc_rms_res_, terminal_res_,
																								 syn_chi_res_list_, const_full_model_info( pose ).extra_minimize_res(),
																								 working_fixed_res_, working_minimize_res_,
																								 allow_insert_, minimize_move_map_ );
	return job_parameters;
}

///////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_Modeler::set_job_parameters( StepWiseRNA_JobParametersCOP job_parameters ){ job_parameters_ = job_parameters;	}

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
		sampling::rna::StepWiseRNA_ModelerOptionsOP new_options = options_->clone();
		new_options->set_skip_sampling( setting );
		set_options( new_options );
	}
}


} //rna
} //sampling
} //stepwise
} //protocols
