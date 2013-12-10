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
#include <protocols/swa/rna/StepWiseRNA_ResidueSampler.hh>
#include <protocols/swa/rna/StepWiseRNA_FloatingBaseSampler.hh>
#include <protocols/swa/rna/StepWiseRNA_StandardResidueSampler.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/swa/rna/StepWiseRNA_VirtualSugarSamplerWrapper.hh>
#include <protocols/swa/rna/StepWiseRNA_OutputData.hh>
#include <protocols/swa/rna/screener/StepWiseRNA_BaseCentroidScreener.hh>
#include <protocols/swa/rna/screener/StepWiseRNA_VDW_BinScreener.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
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

static basic::Tracer TR( "protocols.swa.rna.StepWiseRNA_ResidueSampler" ) ;

namespace protocols {
namespace swa {
namespace rna {

//////////////////////////////////////////////////////////////////////////
//constructor!
StepWiseRNA_ResidueSampler::StepWiseRNA_ResidueSampler( StepWiseRNA_JobParametersCOP & job_parameters ):
	job_parameters_( job_parameters ),
	scorefxn_( core::scoring::ScoreFunctionFactory::create_score_function( "rna_hires.wts" ) ), // can be replaced from the outside
	silent_file_( "silent_file.txt" ),
	output_filename_( "data.txt" ),
	num_pose_kept_( 108 ),
	cluster_rmsd_( 0.5001 ),
	verbose_( false ),
	native_rmsd_screen_( false ),
	native_screen_rmsd_cutoff_( 2.0 ),
	perform_o2prime_pack_( true ),
	use_green_packer_( false ),
	allow_bulge_at_chainbreak_( false ),
	integration_test_mode_( false ), //March 16, 2012
	centroid_screen_( true ),
	VDW_atr_rep_screen_( true ),
	allow_syn_pyrimidine_( false ), //New option Nov 15, 2010
	distinguish_pucker_( true ),
	finer_sampling_at_chain_closure_( false ), //New option Jun 10 2010
	PBP_clustering_at_chain_closure_( false ), //New option Aug 15 2010
	reinitialize_CCD_torsions_( false ), //New option Aug 15 2010 //Reinitialize_CCD_torsion to zero before every CCD chain closure
	extra_epsilon_rotamer_( false ), //New option Aug 30, 2010
	extra_beta_rotamer_( false ), //New option Aug 30, 2010
	extra_chi_( false ),
	sample_both_sugar_base_rotamer_( false ), //New option Nov 12, 2010 (mainly for square_RNA)
	include_torsion_value_in_tag_( false ), //For checking if the extra rotamer are important
	combine_long_loop_mode_( false ), //in this mode, the moving_residues must contact the last residue built from the other side.
	do_not_sample_multiple_virtual_sugar_( false ), //Nov 13, 2010, optimize the chain closure step speed
	sample_ONLY_multiple_virtual_sugar_( false ), //Nov 13, 2010, optimize the chain closure step speed
	assert_no_virt_sugar_sampling_( false ), //July 28 2011
	output_pdb_( false ), //Sept 24, 2011
	choose_random_( false ), // Rhiju, Jul 2013
	num_random_samples_( 1 ),
	force_centroid_interaction_( false ),  // Rhiju, Jul 2013
	use_phenix_geo_( false ),
	virtual_sugar_legacy_mode_( false ),
	kic_sampling_( false )
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
StepWiseRNA_ResidueSampler::apply( core::pose::Pose & pose ) {
	using namespace ObjexxFCL;

	output_title_text( "Enter StepWiseRNA_ResidueSampler::apply", TR.Debug );
	clock_t const time_start( clock() );
	output_options();

	Pose const pose_save = pose; pose = pose_save; //this recopy is actually useful for triggering graphics.
	instantiate_any_virtual_sugars( pose );
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
	floating_base_sampler.set_silent_file ( silent_file_ );
	floating_base_sampler.set_scorefxn ( scorefxn_ );
	floating_base_sampler.set_num_pose_kept ( num_pose_kept_ );
	floating_base_sampler.set_native_rmsd_screen ( native_rmsd_screen_ );
	floating_base_sampler.set_native_screen_rmsd_cutoff ( native_screen_rmsd_cutoff_ );
	floating_base_sampler.set_perform_o2prime_pack ( perform_o2prime_pack_ );
	floating_base_sampler.set_verbose ( verbose_ );
	floating_base_sampler.set_cluster_rmsd (	cluster_rmsd_	);
	floating_base_sampler.set_distinguish_pucker ( distinguish_pucker_ );
	floating_base_sampler.set_PBP_clustering_at_chain_closure ( PBP_clustering_at_chain_closure_ );
	floating_base_sampler.set_use_phenix_geo ( use_phenix_geo_  );
	floating_base_sampler.set_centroid_screen ( centroid_screen_ );
	floating_base_sampler.set_VDW_atr_rep_screen ( VDW_atr_rep_screen_ );
	floating_base_sampler.set_integration_test_mode( integration_test_mode_ ); //Should set after setting sampler_native_screen_rmsd_cutoff, fast, medium_fast options.
	floating_base_sampler.set_reinitialize_CCD_torsions( reinitialize_CCD_torsions_ );
	floating_base_sampler.set_anchor_sugar_modeling( virtual_sugar_sampler_wrapper_->anchor_sugar_modeling() );
	floating_base_sampler.set_choose_random( choose_random_ );
	floating_base_sampler.set_num_random_samples( num_random_samples_ );
	floating_base_sampler.set_try_sugar_instantiation( try_sugar_instantiation_ );
	runtime_assert( !use_green_packer_ );
	runtime_assert( !combine_long_loop_mode_ );

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
	standard_residue_sampler.set_silent_file ( silent_file_ );
	standard_residue_sampler.set_scorefxn ( scorefxn_ );
	standard_residue_sampler.set_num_pose_kept ( num_pose_kept_ );
	standard_residue_sampler.set_native_rmsd_screen ( native_rmsd_screen_ );
	standard_residue_sampler.set_native_screen_rmsd_cutoff ( native_screen_rmsd_cutoff_ );
	standard_residue_sampler.set_perform_o2prime_pack ( perform_o2prime_pack_ );
	standard_residue_sampler.set_use_green_packer ( use_green_packer_ );
	standard_residue_sampler.set_verbose ( verbose_ );
	standard_residue_sampler.set_cluster_rmsd (	cluster_rmsd_	);
	standard_residue_sampler.set_distinguish_pucker ( distinguish_pucker_ );
	standard_residue_sampler.set_finer_sampling_at_chain_closure ( finer_sampling_at_chain_closure_ );
	standard_residue_sampler.set_PBP_clustering_at_chain_closure ( PBP_clustering_at_chain_closure_ );
	standard_residue_sampler.set_allow_syn_pyrimidine( allow_syn_pyrimidine_ );
	standard_residue_sampler.set_extra_chi( extra_chi_ );
	standard_residue_sampler.set_use_phenix_geo ( use_phenix_geo_  );
	standard_residue_sampler.set_kic_sampling( kic_sampling_ );
	standard_residue_sampler.set_centroid_screen ( centroid_screen_ );
	standard_residue_sampler.set_VDW_atr_rep_screen ( VDW_atr_rep_screen_ );
	standard_residue_sampler.set_force_centroid_interaction ( force_centroid_interaction_ );
	standard_residue_sampler.set_integration_test_mode( integration_test_mode_ ); //Should set after setting sampler_native_screen_rmsd_cutoff, fast, medium_fast options.
	standard_residue_sampler.set_allow_bulge_at_chainbreak( allow_bulge_at_chainbreak_ );
	standard_residue_sampler.set_parin_favorite_output( parin_favorite_output_ );
	standard_residue_sampler.set_reinitialize_CCD_torsions( reinitialize_CCD_torsions_ );
	standard_residue_sampler.set_extra_epsilon_rotamer( extra_epsilon_rotamer_ );
	standard_residue_sampler.set_extra_beta_rotamer( extra_beta_rotamer_ );
	standard_residue_sampler.set_sample_both_sugar_base_rotamer( sample_both_sugar_base_rotamer_ ); //Nov 12, 2010
	standard_residue_sampler.set_include_torsion_value_in_tag( include_torsion_value_in_tag_ );
	standard_residue_sampler.set_combine_long_loop_mode( combine_long_loop_mode_ );
	standard_residue_sampler.set_choose_random( choose_random_ );
	standard_residue_sampler.set_num_random_samples( num_random_samples_ );

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
void
StepWiseRNA_ResidueSampler::instantiate_any_virtual_sugars( pose::Pose & pose ){
	Pose pose_save = pose;
	virtual_sugar_sampler_wrapper_ = new StepWiseRNA_VirtualSugarSamplerWrapper( job_parameters_ );
	virtual_sugar_sampler_wrapper_->set_scorefxn( scorefxn_ );
	virtual_sugar_sampler_wrapper_->set_do_not_sample_multiple_virtual_sugar( do_not_sample_multiple_virtual_sugar_ );
	virtual_sugar_sampler_wrapper_->set_sample_ONLY_multiple_virtual_sugar( sample_ONLY_multiple_virtual_sugar_ );
	virtual_sugar_sampler_wrapper_->set_assert_no_virt_sugar_sampling( assert_no_virt_sugar_sampling_ );
	virtual_sugar_sampler_wrapper_->set_use_phenix_geo ( use_phenix_geo_  );
	virtual_sugar_sampler_wrapper_->set_legacy_mode ( virtual_sugar_legacy_mode_  );
	virtual_sugar_sampler_wrapper_->set_integration_test_mode( integration_test_mode_ ); //Should set after setting sampler_native_screen_rmsd_cutoff, fast, medium_fast options.
	virtual_sugar_sampler_wrapper_->apply( pose );
	pose = pose_save;
}



////////////////////////////////////////////////////////////////////////////////////////

void
StepWiseRNA_ResidueSampler::set_user_input_VDW_bin_screener( screener::StepWiseRNA_VDW_BinScreenerOP const & user_input_VDW_bin_screener ){ user_input_VDW_bin_screener_ = user_input_VDW_bin_screener; }




////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< PoseOP > &
StepWiseRNA_ResidueSampler::get_pose_list(){
	return pose_list_;
}

///////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::set_silent_file( std::string const & silent_file ){
	silent_file_ = silent_file;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::set_num_pose_kept( core::Size const & num_pose_kept ){
	num_pose_kept_ = num_pose_kept ;
}
//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::set_native_rmsd_screen( bool const & setting ){
	native_rmsd_screen_ = setting;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::set_native_screen_rmsd_cutoff( core::Real const & setting ){
	native_screen_rmsd_cutoff_ = setting;
}
//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::set_integration_test_mode( bool const & setting ){
	integration_test_mode_ = setting;
	if ( integration_test_mode_ ){
		num_pose_kept_ = 2;
		native_rmsd_screen_ = false; // will be switched to true mid-way through sampling.
		native_screen_rmsd_cutoff_ = 1.0;
	}
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::set_verbose( bool const & setting ){
	verbose_ = setting;
}
//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::set_perform_o2prime_pack( bool const & setting ){
	perform_o2prime_pack_ = setting;
}


//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::set_allow_bulge_at_chainbreak( bool const & setting ){
	allow_bulge_at_chainbreak_ = setting;
}


//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::set_output_filename( std::string const & output_filename ){
	output_filename_ = output_filename;
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

	if ( verbose_ == false ){ //consistency check Apr 3, 2010
		utility_exit_with_message( "verbose_ == false, but StepWiseRNA_ResidueSampler::output_pose_list is still called?!" );
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
StepWiseRNA_ResidueSampler::set_cluster_rmsd( Real const & setting ){
	cluster_rmsd_ = setting;
	TR.Debug << "Set cluster_rmsd to " << cluster_rmsd_ << std::endl;
}


//////////////////////////////////////////////////////////////////
void
StepWiseRNA_ResidueSampler::output_options(){
	//output screen options
	TR.Debug << "--------SCREEN OPTIONS---------- " << std::endl;
	output_boolean( "integration_test_mode_ = ", integration_test_mode_, TR.Debug ); TR.Debug << std::endl;
	output_boolean( "native_rmsd_screen = ", native_rmsd_screen_, TR.Debug ); TR.Debug << std::endl;
	TR.Debug << "native_screen_rmsd_cutoff = " << native_screen_rmsd_cutoff_ << std::endl;
	output_boolean( "perform_o2prime_pack = ", perform_o2prime_pack_, TR.Debug ); TR.Debug << std::endl;
	output_seq_num_list( "working_moving_partition_pos = ", job_parameters_->working_moving_partition_pos(), TR.Debug );
	output_boolean( "centroid_screen = ", centroid_screen_, TR.Debug ); TR.Debug << std::endl;
	output_boolean( "VDW_atr_rep_screen = ", VDW_atr_rep_screen_, TR.Debug ); TR.Debug << std::endl;
	output_boolean( "sample_both_sugar_base_rotamer_ = ", sample_both_sugar_base_rotamer_, TR.Debug ); TR.Debug << std::endl;
	output_boolean( "do_not_sample_multiple_virtual_sugar_ = ", do_not_sample_multiple_virtual_sugar_, TR.Debug ); TR.Debug << std::endl;
	output_boolean( "sample_ONLY_multiple_virtual_sugar_ = ", sample_ONLY_multiple_virtual_sugar_, TR.Debug ); TR.Debug << std::endl;
	output_boolean( "assert_no_virt_sugar_sampling_ = ", assert_no_virt_sugar_sampling_, TR.Debug ); TR.Debug << std::endl;
	output_boolean( "distinguish_pucker_ ", distinguish_pucker_, TR.Debug ); TR.Debug << std::endl;
	TR.Debug << "--------------------------------" << std::endl;
}

}
}
}
