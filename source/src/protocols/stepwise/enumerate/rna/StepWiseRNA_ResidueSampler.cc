// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_ResidueSampler
/// @brief Master sampler for StepWise Assembly of RNA.
/// @detailed
/// @author Parin Sripakdeevong
/// @author Rhiju Das

//////////////////////////////////
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_ResidueSampler.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_FloatingBaseSampler.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_StandardResidueSampler.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_ModelerOptions.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_VirtualSugarSamplerWrapper.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_OutputData.hh>
#include <protocols/stepwise/enumerate/rna/screener/StepWiseRNA_BaseCentroidScreener.hh>
#include <protocols/stepwise/enumerate/rna/screener/StepWiseRNA_VDW_BinScreener.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_Util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <time.h>

using namespace core;

//////////////////////////////////////////////////////////////////////////
// Core routine for stepwise sampling of RNA
//
//  This carries out standard enumerative sampling of a nucleotide and its
//   connections to the next base OR floating base sampling, which skips
//   bulges (a.k.a., 'dinucleotide move')
//  Prior to these moves, any virtualized sugars at the residue to be sampled
//   or at residues involved in chain closure are instantiated.
//
//////////////////////////////////////////////////////////////////////////

static basic::Tracer TR( "protocols.stepwise.rna.StepWiseRNA_ResidueSampler" ) ;

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace rna {

//////////////////////////////////////////////////////////////////////////
//constructor!
StepWiseRNA_ResidueSampler::StepWiseRNA_ResidueSampler( StepWiseRNA_JobParametersCOP & job_parameters ):
	job_parameters_( job_parameters ),
	scorefxn_( core::scoring::ScoreFunctionFactory::create_score_function( "rna_hires.wts" ) ), // can be replaced from the outside
	options_( new StepWiseRNA_ModelerOptions ),
	sampling_silent_file_( "default.out" )
{
	set_native_pose( job_parameters_->working_native_pose() );
}

//////////////////////////////////////////////////////////////////////////
//destructor
StepWiseRNA_ResidueSampler::~StepWiseRNA_ResidueSampler()
{}

/////////////////////
std::string
StepWiseRNA_ResidueSampler::get_name() const {
	return "StepWiseRNA_ResidueSampler";
}

////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::set_options( StepWiseRNA_ModelerOptionsCOP options ){
	StepWiseRNA_ModelerOptionsOP sampler_options = options->clone();
	if ( sampler_options->choose_random() )	sampler_options->set_cluster_rmsd( 0.0 ); // don't cluster.
	if ( sampler_options->integration_test_mode() ){
		sampler_options->set_sampler_num_pose_kept( 2 );
		sampler_options->set_sampler_native_rmsd_screen( false ); // will be switched to true mid-way through sampling.
		sampler_options->set_sampler_native_screen_rmsd_cutoff( 1.0 );
	}
	options_ = sampler_options;
}

////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::apply( core::pose::Pose & pose ) {
	using namespace ObjexxFCL;

	output_title_text( "Enter StepWiseRNA_ResidueSampler::apply", TR.Debug );
	clock_t const time_start( clock() );
	output_options();

	Pose const pose_save = pose; pose = pose_save; //this recopy is actually useful for triggering graphics.
	if ( !instantiate_any_virtual_sugars( pose ) ) return;
	if ( job_parameters_->floating_base() ){
		floating_base_sampling( pose );
	} else{
		standard_sampling_WRAPPER( pose );
	}

	pose = pose_save;

	TR.Debug << "Total time in StepWiseRNA_ResidueSampler::apply " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::floating_base_sampling( pose::Pose & pose ){

	StepWiseRNA_FloatingBaseSampler floating_base_sampler( job_parameters_ );
	floating_base_sampler.set_base_centroid_screener( base_centroid_screener_ );
	floating_base_sampler.set_user_input_VDW_bin_screener ( user_input_VDW_bin_screener_ );
	floating_base_sampler.set_silent_file ( sampling_silent_file_ );
	floating_base_sampler.set_scorefxn ( scorefxn_ );
	floating_base_sampler.set_anchor_sugar_modeling( virtual_sugar_sampler_wrapper_->anchor_sugar_modeling() );
	floating_base_sampler.set_options( options_ );
	runtime_assert( !options_->use_green_packer() );
	runtime_assert( !options_->combine_long_loop_mode() );

	floating_base_sampler.apply( pose );

	pose_list_ = floating_base_sampler.get_pose_list();

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
StepWiseRNA_ResidueSampler::standard_sampling_WRAPPER( core::pose::Pose & pose ){

	using namespace ObjexxFCL;
	using namespace core::io::silent;
	using namespace core::id;
	using namespace core::scoring;

	// extraneous?
	pose::Pose const pose_save = pose;

	StepWiseRNA_StandardResidueSampler standard_residue_sampler( job_parameters_ );
	standard_residue_sampler.set_pose_list( pose_list_ ); // allows for accumulation of poses, if desired.
	standard_residue_sampler.set_base_centroid_screener( base_centroid_screener_ );
	standard_residue_sampler.set_user_input_VDW_bin_screener ( user_input_VDW_bin_screener_ );
	standard_residue_sampler.set_scorefxn ( scorefxn_ );
	standard_residue_sampler.set_options( options_ );
	standard_residue_sampler.set_silent_file( sampling_silent_file_ );

	if ( !virtual_sugar_sampler_wrapper_->sampling_sugar() ){
		standard_residue_sampler.apply( pose );
	} else{ //Case where we have to sample virtual sugar...
		utility::vector1< PoseOP > pose_list;
		virtual_sugar_sampler_wrapper_->prepare_from_prior_sampled_sugar_jobs( pose, pose_list );
		for ( Size n = 1; n <= pose_list.size(); n++ ){
			TR << TR.Blue << "Running on sugar job " << n << " out of " << pose_list.size() << TR.Reset << std::endl;
			pose = ( *pose_list[n] ); //set viewer_pose;
			standard_residue_sampler.set_extra_tag( tag_from_pose( *pose_list[n] )  );
			standard_residue_sampler.apply( pose );
		}
	}

	pose_list_ = standard_residue_sampler.pose_list();
	pose = pose_save;

}


//////////////////////Build previously virtualize sugar/////////////////////
bool
StepWiseRNA_ResidueSampler::instantiate_any_virtual_sugars( pose::Pose & pose ){
	Pose pose_save = pose;
	virtual_sugar_sampler_wrapper_ = new StepWiseRNA_VirtualSugarSamplerWrapper( job_parameters_ );
	virtual_sugar_sampler_wrapper_->set_scorefxn( scorefxn_ );
	virtual_sugar_sampler_wrapper_->set_options( options_ );
	virtual_sugar_sampler_wrapper_->apply( pose );
	pose = pose_save;
	return virtual_sugar_sampler_wrapper_->success();
}


////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::set_user_input_VDW_bin_screener( screener::StepWiseRNA_VDW_BinScreenerOP const & user_input_VDW_bin_screener ){
	user_input_VDW_bin_screener_ = user_input_VDW_bin_screener;
}

////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< PoseOP > &
StepWiseRNA_ResidueSampler::get_pose_list(){
	return pose_list_;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn ){
	scorefxn_ = scorefxn;
}

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::output_pose_list( std::string const final_sampler_output_silent_file ) const{
	using namespace core::io::silent;

	if ( options_->verbose() == false ){ //consistency check Apr 3, 2010
		utility_exit_with_message( "options_->verbose() == false, but StepWiseRNA_ResidueSampler::output_pose_list is still called?!" );
	}

	SilentFileData silent_file_data;

	for ( Size n = 1; n <= pose_list_.size(); n++ ) {
		output_data( silent_file_data, final_sampler_output_silent_file, tag_from_pose( *pose_list_[n] ), false, *( pose_list_[n] ), get_native_pose(), job_parameters_ );
	}

}

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::set_base_centroid_screener( screener::StepWiseRNA_BaseCentroidScreenerOP & screener ){
	base_centroid_screener_ = screener;
}

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::output_options(){
	//output screen options
	TR.Debug << "--------SCREEN OPTIONS---------- " << std::endl;
	output_boolean( "integration_test_mode_ = ", options_->integration_test_mode(), TR.Debug ); TR.Debug << std::endl;
	output_boolean( "native_rmsd_screen = ", options_->sampler_native_rmsd_screen(), TR.Debug ); TR.Debug << std::endl;
	TR.Debug << "native_screen_rmsd_cutoff = " << options_->sampler_native_screen_rmsd_cutoff() << std::endl;
	output_boolean( "perform_o2prime_pack = ", options_->sampler_perform_o2prime_pack(), TR.Debug ); TR.Debug << std::endl;
	output_seq_num_list( "working_moving_partition_pos = ", job_parameters_->working_moving_partition_pos(), TR.Debug );
	output_boolean( "centroid_screen = ", (base_centroid_screener_ != 0), TR.Debug ); TR.Debug << std::endl;
	output_boolean( "VDW_atr_rep_screen = ", options_->VDW_atr_rep_screen(), TR.Debug ); TR.Debug << std::endl;
	output_boolean( "sample_both_sugar_base_rotamer_ = ", options_->sample_both_sugar_base_rotamer(), TR.Debug ); TR.Debug << std::endl;
	output_boolean( "do_not_sample_multiple_virtual_sugar_ = ", options_->do_not_sample_multiple_virtual_sugar(), TR.Debug ); TR.Debug << std::endl;
	output_boolean( "sample_ONLY_multiple_virtual_sugar_ = ", options_->sample_ONLY_multiple_virtual_sugar(), TR.Debug ); TR.Debug << std::endl;
	output_boolean( "assert_no_virt_sugar_sampling_ = ", options_->sampler_assert_no_virt_sugar_sampling(), TR.Debug ); TR.Debug << std::endl;
	output_boolean( "distinguish_pucker_ ", options_->distinguish_pucker(), TR.Debug ); TR.Debug << std::endl;
	TR.Debug << "--------------------------------" << std::endl;
}


} //rna
} //enumerate
} //stepwise
} //protocols
