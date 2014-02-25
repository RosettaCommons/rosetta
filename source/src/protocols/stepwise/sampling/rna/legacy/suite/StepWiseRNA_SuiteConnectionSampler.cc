// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/rna/legacy/suite/StepWiseRNA_SuiteConnectionSampler.cc
/// @brief Sampler of a single residue, the heart of StepWise Assembly of RNA
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu
/// @author Parin Sripakdeevong, sripakpa@stanford.edu


#include <protocols/stepwise/sampling/rna/legacy/suite/StepWiseRNA_SuiteConnectionSampler.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_ModelerOptions.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_PoseSelection.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_OutputData.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_Util.hh>
#include <protocols/stepwise/sampling/rna/o2prime/O2PrimePacker.hh>
#include <protocols/stepwise/sampling/rna/phosphate/MultiPhosphateSampler.hh>
#include <protocols/stepwise/sampling/rna/phosphate/PhosphateUtil.hh>
#include <protocols/stepwise/sampling/rna/checker/RNA_VDW_BinChecker.hh>
#include <protocols/stepwise/sampling/rna/checker/AtrRepChecker.hh>
#include <protocols/stepwise/sampling/rna/checker/ChainClosableGeometryChecker.hh>
#include <protocols/stepwise/sampling/rna/checker/ChainClosureChecker.hh>
#include <protocols/stepwise/sampling/rna/checker/RNA_BaseCentroidChecker.hh>
#include <protocols/stepwise/StepWiseSampleAndScreen.hh>
#include <protocols/stepwise/screener/AtrRepScreener.hh>
#include <protocols/stepwise/screener/BaseCentroidScreener.hh>
#include <protocols/stepwise/screener/BulgeApplier.hh>
#include <protocols/stepwise/screener/ChainClosableGeometryScreener.hh>
#include <protocols/stepwise/screener/ChainClosureScreener.hh>
#include <protocols/stepwise/screener/IntegrationTestBreaker.hh>
#include <protocols/stepwise/screener/NativeRMSD_Screener.hh>
#include <protocols/stepwise/screener/O2PrimeScreener.hh>
#include <protocols/stepwise/screener/PhosphateScreener.hh>
#include <protocols/stepwise/screener/PoseSelectionScreener.hh>
#include <protocols/stepwise/screener/ResidueContactScreener.hh>
#include <protocols/stepwise/screener/SampleApplier.hh>
#include <protocols/stepwise/screener/StepWiseScreener.hh>
#include <protocols/stepwise/screener/TagDefinition.hh>
#include <protocols/stepwise/screener/VDW_BinScreener.hh>
#include <protocols/rotamer_sampler/RotamerBase.hh>
#include <protocols/rotamer_sampler/rna/RNA_RotamerSamplerUtil.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/rna/RNA_Util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_Util.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

using ObjexxFCL::string_of;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Enumerates over suite backbone, sugar pucker, and base chi angle; and packs 2' hydroxyls, for a single RNA residue
// at the terminus of an RNA pose.
//
//                 epsilon,zeta,alpha,
//                    beta,gamma      pucker
//
//     5'   --Sugar --  Phosphate -- Sugar  ... -- Phosphate -- Sugar -- Phosphate -- 3'
//              |                      | chi                       |
//          Reference               Moving                       Distal
//             Base                 Residue                      Base
//                                                        |
//                                       If no residues to distal, close_chain_to_distal
//
//                                     |<------ gap_size + 1 ----->|
//
//
// Will soon be unified with rigid body sampling!
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static basic::Tracer TR( "protocols.stepwise.rna.StepWiseRNA_SuiteConnectionSampler" ) ;

namespace protocols {
namespace stepwise {
namespace sampling {
namespace rna {
namespace legacy {
namespace suite {

	//Constructor
	StepWiseRNA_SuiteConnectionSampler::StepWiseRNA_SuiteConnectionSampler( StepWiseRNA_JobParametersCOP & job_parameters ):
		job_parameters_( job_parameters ),
		moving_res_(  job_parameters_->working_moving_res() ), // Might not corresponds to user input.
		moving_suite_(  job_parameters_->working_moving_suite() ), // dofs betweeen this value and value+1 actually move.
		is_prepend_(  job_parameters_->is_prepend() ),
		is_internal_(  job_parameters_->is_internal() ), // no cutpoints before or after moving_res_.
		actually_moving_res_( job_parameters_->actually_moving_res() ), //Now same as moving_res_
		gap_size_( job_parameters_->gap_size() ), /* If this is zero or one, need to screen or closable chain break */
		working_moving_partition_pos_( job_parameters_->working_moving_partition_pos() ),
		num_nucleotides_(  job_parameters_->working_moving_res_list().size() ),
		is_dinucleotide_( num_nucleotides_ == 2 ),
		close_chain_to_distal_( gap_size_ == 0 ),
		five_prime_chain_break_res_( job_parameters_->five_prime_chain_break_res() ),
		last_append_res_( ( is_prepend_ ) ? moving_res_ - 1: moving_res_ ),
		last_prepend_res_( ( is_prepend_ ) ? moving_res_: moving_res_ + 1 ),
		atom_atom_overlap_dist_cutoff_(-1.0 ),
		extra_tag_( "" ),
		build_pose_from_scratch_( job_parameters_->working_sequence().length() == ( num_nucleotides_ + 1 ) ), // somewhat hacky, used for rna puzzle
		kic_sampling_( false ), // will be updated below.
		rebuild_bulge_mode_( job_parameters_->rebuild_bulge_mode() ),
		scorefxn_( core::scoring::ScoreFunctionFactory::create_score_function( "rna_hires.wts" ) ), // can be replaced from the outside
		silent_file_( "silent_file.txt" ),
		native_rmsd_screen_( false ), // will be updated below
		bin_size_( 20 )
	{
		set_native_pose( job_parameters_->working_native_pose() );
		if ( num_nucleotides_ != 1 && num_nucleotides_ != 2 )		utility_exit_with_message( "num_nucleotides_ != 1 and num_nucleotides_ != 2" );
	}

	//Destructor
	StepWiseRNA_SuiteConnectionSampler::~StepWiseRNA_SuiteConnectionSampler()
	{}

		/////////////////////
	std::string
	StepWiseRNA_SuiteConnectionSampler::get_name() const {
		return "StepWiseRNA_SuiteConnectionSampler";
	}

	////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_SuiteConnectionSampler::apply( core::pose::Pose & pose ) {

		using namespace core::scoring;
		using namespace core::pose;
		using namespace core::io::silent;
		using namespace core::id;

		output_title_text( "Enter StepWiseRNA_SuiteConnectionSampler::standard_sampling", TR.Debug );
		clock_t const time_start( clock() );

		initialize_poses_and_checkers( pose );
		initialize_screeners( pose );
		StepWiseSampleAndScreen sample_and_screen( rotamer_sampler_, screeners_ );
		sample_and_screen.set_max_ntries( get_max_ntries() );

		// Do it!
		TR << "KICKING OFF SAMPLE AND SCREEN!!!!!!!!!!" << std::endl;
		sample_and_screen.run();
		sample_and_screen.output_counts();
		sample_and_screen.output_info_on_random_trials();

		pose_selection_->finalize( !build_pose_from_scratch_ /*do_clustering*/ );
		pose_list_ = pose_selection_->pose_list();

		TR.Debug << "Total time in StepWiseRNA_SuiteConnectionSampler: " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;
	}

	//////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////
	// Create screeners. This actually defines the 'main loop'.
	//////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_SuiteConnectionSampler::initialize_screeners( pose::Pose & pose ){

		using namespace screener;

		screeners_.clear();

		screeners_.push_back( new SampleApplier( *screening_pose_ ) );

		NativeRMSD_ScreenerOP native_rmsd_screener = new NativeRMSD_Screener( *get_native_pose(), *screening_pose_,
																																					 job_parameters_, options_->sampler_native_screen_rmsd_cutoff(),
																																					 native_rmsd_screen_ /*do_screen*/ );
		screeners_.push_back( native_rmsd_screener );

		// For KIC closure, immediate check that a closed loop solution was actually found.
		if ( kic_sampling_ ) {
			screeners_.push_back( new ChainClosureScreener( chain_closure_checker_, *screening_pose_, true /*just do closure check*/ ) );
		}

		if ( options_->combine_long_loop_mode() && !close_chain_to_distal_ ) { //residue-residue contact screen;
			screeners_.push_back( new ResidueContactScreener( *screening_pose_, last_append_res_,  last_prepend_res_, atom_atom_overlap_dist_cutoff_ ) );
		}

		if ( base_centroid_checker_ ){
			bool const force_centroid_interaction = (options_->force_centroid_interaction() || !close_chain_to_distal_);
			screeners_.push_back( new BaseCentroidScreener( base_centroid_checker_, screening_pose_, force_centroid_interaction ) );
		}

		if ( close_chain_to_distal_ && !kic_sampling_ ) {
			screeners_.push_back( new ChainClosableGeometryScreener( chain_closable_geometry_checker_, screening_pose_,
																											 options_->finer_sampling_at_chain_closure() ) );
		}

		AtrRepScreenerOP atr_rep_screener = ( options_->VDW_atr_rep_screen() ) ? new AtrRepScreener( atr_rep_checker_, *screening_pose_ ) : 0;
		screeners_.push_back( atr_rep_screener );

		// shouldn't this be before atr/rep screener?
		if ( ( user_input_VDW_bin_checker_->user_inputted_VDW_screen_pose() ) && ( gap_size_ != 0 ) && ( !is_internal_ ) ){
			screeners_.push_back( new VDW_BinScreener( user_input_VDW_bin_checker_, *screening_pose_, moving_res_ ) );
		}

		if ( close_chain_to_distal_ && !kic_sampling_ ){
			screeners_.push_back( new ChainClosureScreener( chain_closure_checker_ ) );
		}

		if ( options_->sampler_perform_phosphate_pack() ) {
			screeners_.push_back( new PhosphateScreener( phosphate_sampler_ ) );
		}

		if ( options_->sampler_perform_o2prime_pack() ) {
			screeners_.push_back( new O2PrimeScreener( o2prime_packer_ ) );
		}

		screeners_.push_back( new SampleApplier( pose ) );

		TagDefinitionOP tag_definition = new TagDefinition( pose, screeners_[1], options_->sampler_include_torsion_value_in_tag(),
																						moving_res_, is_prepend_, extra_tag_ );
		screeners_.push_back( tag_definition );

		if ( !rebuild_bulge_mode_ && options_->allow_bulge_at_chainbreak() && working_moving_partition_pos_.size() == 1
				 && close_chain_to_distal_ ) screeners_.push_back( new BulgeApplier( atr_rep_checker_, base_centroid_checker_, moving_res_ ) ); // apply bulge at the last minute.

		screeners_.push_back( new PoseSelectionScreener( pose_selection_, pose /*const reference*/, tag_definition,
																										 options_->verbose(), silent_file_, get_native_pose(), job_parameters_ ) );

		if ( options_->integration_test_mode() ) {
			screeners_.insert( screeners_.begin() /*right at beginning!*/,
												 new IntegrationTestBreaker( atr_rep_screener, screeners_[ screeners_.size() ], native_rmsd_screener ) );
		}

	}


	////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_SuiteConnectionSampler::initialize_poses_and_checkers( pose::Pose & pose ){

		if ( options_->combine_long_loop_mode() && !close_chain_to_distal_ ){ //residue-residue contact screen;
			TR.Debug << "options_->combine_long_loop_mode() && gap_size_ == 0" << std::endl;
			TR.Debug << "Enforcing contact between LAST_APPEND_RES: " << last_append_res_ << " and LAST_PREPEND_RES: " << last_prepend_res_  << std::endl;
			TR.Debug << "atom_atom_overlap_dist_cutoff " << atom_atom_overlap_dist_cutoff_ << std::endl;
		}

		if ( options_->sampler_perform_phosphate_pack() ){
			phosphate_sampler_ = new phosphate::MultiPhosphateSampler( pose );
			phosphate_sampler_->set_moving_partition_res( working_moving_partition_pos_ );
		}

		/////////////////////////////// O2prime sampling/virtualization //////////////////////////
		Pose pose_with_virtual_O2prime_hydrogen = pose;
		add_virtual_O2Prime_hydrogen( pose_with_virtual_O2prime_hydrogen );

		if ( options_->sampler_perform_o2prime_pack() ) {
			o2prime_packer_ = new o2prime::O2PrimePacker( pose, scorefxn_, working_moving_partition_pos_ /* moving_res_*/ );
			o2prime_packer_->set_use_green_packer( options_->use_green_packer() ); // green_packer not on by default
			o2prime_packer_->set_partition_definition( job_parameters_->partition_definition() ); // would be needed for green packer
		} else {
			// Otherwise, virtualize the 2-OH.
			pose = pose_with_virtual_O2prime_hydrogen;
		}

		kic_sampling_ = ( options_->kic_sampling_if_relevant() && close_chain_to_distal_ );

		screening_pose_ = pose_with_virtual_O2prime_hydrogen.clone(); //Hard copy
		phosphate::remove_terminal_phosphates( *screening_pose_ );

		if ( close_chain_to_distal_ )	pose::add_variant_type_to_pose_residue( *screening_pose_,      "VIRTUAL_PHOSPHATE",    five_prime_chain_break_res_ + 1 );

		atr_rep_checker_ = new checker::AtrRepChecker( *screening_pose_, job_parameters_ );
		atr_rep_checker_->set_kic_sampling( kic_sampling_ );

		chain_closable_geometry_checker_ = new checker::ChainClosableGeometryChecker( five_prime_chain_break_res_, gap_size_ );

		chain_closure_checker_ = new checker::ChainClosureChecker( pose_with_virtual_O2prime_hydrogen, five_prime_chain_break_res_ ); //Hard copy
		chain_closure_checker_->set_reinitialize_CCD_torsions( options_->reinitialize_CCD_torsions() );

		if ( close_chain_to_distal_ ) add_harmonic_chain_break_constraint( pose, five_prime_chain_break_res_ );

		pose_selection_ = new StepWiseRNA_PoseSelection( job_parameters_, scorefxn_ );
		pose_selection_->set_num_pose_kept( get_num_pose_kept() );
		pose_selection_->set_cluster_rmsd( options_->cluster_rmsd() );
		pose_selection_->set_PBP_clustering_at_chain_closure( options_->PBP_clustering_at_chain_closure() );
		pose_selection_->set_distinguish_pucker( options_->distinguish_pucker() );
		// Allows for accumulation and clustering of poses across multiple jobs with e.g. different conformations of sampled sugars.
		pose_selection_->set_pose_list( pose_list_ );

		native_rmsd_screen_ = options_->sampler_native_rmsd_screen();

		/////Get the Rotamer Sampler/////
		runtime_assert( base_centroid_checker_ );
		rotamer_sampler_ = rotamer_sampler::rna::setup_rotamer_sampler( *screening_pose_, options_,
																																		job_parameters_, build_pose_from_scratch_,
																																		kic_sampling_, close_chain_to_distal_ );
		set_extra_tag( tag_from_pose( pose ) );
	}



	/////////////////////////////////////////////////////////////////////////////////
	Size
	StepWiseRNA_SuiteConnectionSampler::get_num_pose_kept(){

		Size num_pose_kept = options_->sampler_num_pose_kept();

		if ( build_pose_from_scratch_ ){
			TR.Debug << "Since build_pose_from_scratch, choose to increase NUM_POSE_KEPT by 36 fold. ";
			TR.Debug << " Somewhat hacky..since sample both sugar..want to make sure that we keep are good energy score states " << std::endl;
			TR.Debug << "Old_num_pose_kept = " << num_pose_kept  << std::endl;
			num_pose_kept = 36* num_pose_kept;
			TR.Debug << "New_num_pose_kept = " << num_pose_kept  << std::endl;
		}

		if ( job_parameters_->sample_both_sugar_base_rotamer() ){
			TR.Debug << "Since build_pose_from_scratch, choose to increase NUM_POSE_KEPT by 12 fold. ";
			TR.Debug << " Somewhat hacky..since sample both sugar..want to make sure that we keep are good energy score states " << std::endl;
			TR.Debug << "Old_num_pose_kept = " << num_pose_kept  << std::endl;
			num_pose_kept = 12* num_pose_kept;
			TR.Debug << "New_num_pose_kept = " << num_pose_kept  << std::endl;
		}

		return num_pose_kept;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	Size
	StepWiseRNA_SuiteConnectionSampler::get_max_ntries(){
		Size max_ntries = std::max( 10000, 100 * int( options_->num_random_samples() ) );
		if ( kic_sampling_ ) max_ntries = 5 * options_->num_random_samples(); // some chains just aren't closable.
		return max_ntries;
	}

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_SuiteConnectionSampler::set_base_centroid_checker( checker::RNA_BaseCentroidCheckerOP & checker ){
	base_centroid_checker_ = checker;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_SuiteConnectionSampler::set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn ){
	scorefxn_ = scorefxn;
}

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_SuiteConnectionSampler::set_user_input_VDW_bin_checker( checker::RNA_VDW_BinCheckerOP const & user_input_VDW_bin_checker ){ user_input_VDW_bin_checker_ = user_input_VDW_bin_checker; }

//////////////////////////////////////////////////////////////////
utility::vector1< pose::PoseOP > &
StepWiseRNA_SuiteConnectionSampler::pose_list(){
	return pose_list_;
}

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_SuiteConnectionSampler::set_pose_list( utility::vector1< pose::PoseOP > &	pose_list ){
	pose_list_ = pose_list;
}

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_SuiteConnectionSampler::set_options( StepWiseRNA_ModelerOptionsCOP options ){
	options_ = options;
}

} //suite
} //legacy
} //rna
} //sampling
} //stepwise
} //protocols
