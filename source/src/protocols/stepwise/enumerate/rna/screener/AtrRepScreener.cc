// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/enumerate/rna/screener/AtrRepScreener.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/enumerate/rna/screener/AtrRepScreener.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_Util.hh>
#include <protocols/stepwise/StepWiseUtil.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>

#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

static basic::Tracer TR( "protocols.stepwise.rna.screener.AtrRepScreener" ) ;

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace rna {
namespace screener {

	//Constructor
	AtrRepScreener::AtrRepScreener( pose::Pose const & pose,
																	StepWiseRNA_JobParametersCOP & job_parameters ):
		working_moving_res_(     job_parameters->working_moving_res() ),
		working_reference_res_(  job_parameters->working_reference_res() ),
		gap_size_(    job_parameters->gap_size() ),
		is_prepend_(  job_parameters->is_prepend() ),
		is_internal_(  job_parameters->is_internal() ),
		sample_both_sugar_base_rotamer_( job_parameters->sample_both_sugar_base_rotamer() ),
		separate_moving_residue_to_estimate_baseline_( true )
	{
		initialize_parameters();
		initialize_scorefxn();
		get_base_atr_rep_score( pose );
	}

	//Constructor
	AtrRepScreener::AtrRepScreener( pose::Pose const & pose,
																	Size const moving_res,
																	Size const reference_res,
																	Size const gap_size,
																	bool const is_internal /* = false */,
																	bool const separate_moving_residue_to_estimate_baseline, /* = true */
																	bool const sample_both_sugar_base_rotamer /* = false */
																	):
		working_moving_res_( moving_res    ),
		working_reference_res_( reference_res ),
		gap_size_( gap_size ),
		is_prepend_(  working_reference_res_ > working_moving_res_ ),
		is_internal_( is_internal ),
		sample_both_sugar_base_rotamer_( sample_both_sugar_base_rotamer ),
		separate_moving_residue_to_estimate_baseline_( separate_moving_residue_to_estimate_baseline  )
	{
		initialize_parameters();
		initialize_scorefxn();
		get_base_atr_rep_score( pose );
	}

	//Destructor
	AtrRepScreener::~AtrRepScreener()
	{}

	///////////////////////////////////////////
	void
	AtrRepScreener::initialize_scorefxn(){
		// Bare minimum to check for contact (fa_atr) but not clash (fa_rep)
		atr_rep_screening_scorefxn_ =  new core::scoring::ScoreFunction;
		atr_rep_screening_scorefxn_->set_weight( core::scoring::fa_atr , 0.23 );
		atr_rep_screening_scorefxn_->set_weight( core::scoring::fa_rep , 0.12 );
	}

	///////////////////////////////////////////
	void
	AtrRepScreener::initialize_parameters(){
		rep_cutoff_ = 4.0; // const, but can be updated below.
		base_atr_score_ = 0.0;
		base_rep_score_ = 0.0;
		delta_atr_score_ = 0.0;
		delta_rep_score_ = 0.0;
		verbose_ = false;
		output_pdb_ = false;
		kic_sampling_ = false;
	}

	///////////////////////////////////////////
	void
	AtrRepScreener::get_base_atr_rep_score( core::pose::Pose const & pose ){

		using namespace core::conformation;
		using namespace core::scoring;
		using namespace core::pose;
		using namespace ObjexxFCL;

		Size const nres = pose.total_residue();

		///////////////////////////////Old_way////////////////////////////////////////////

		pose::Pose base_pose_screen = pose; //hard copy

		if ( output_pdb_ ) base_pose_screen.dump_pdb( "base_atr_rep_before.pdb" );

		pose::remove_variant_type_from_pose_residue( base_pose_screen, "FIVE_PRIME_PHOSPHATE", working_moving_res_ );
		pose::add_variant_type_to_pose_residue( base_pose_screen, "VIRTUAL_PHOSPHATE", working_moving_res_ ); //May 7...
		if ( ( working_moving_res_ + 1 ) <= nres ){
			pose::remove_variant_type_from_pose_residue( base_pose_screen, "FIVE_PRIME_PHOSPHATE", working_moving_res_ + 1 );
			pose::add_variant_type_to_pose_residue( base_pose_screen, "VIRTUAL_PHOSPHATE", working_moving_res_ + 1 ); //May 7...
		}

		if ( sample_both_sugar_base_rotamer_ ){ //Nov 15, 2010
			Size const extra_sample_sugar_base_res = ( is_prepend_ ) ? ( working_moving_res_ + 1 ) : ( working_moving_res_ - 1 );
			TR << "base_pose_screen extra_sample_sugar_base_res = " << extra_sample_sugar_base_res << std::endl;
			pose::add_variant_type_to_pose_residue( base_pose_screen, "VIRTUAL_RIBOSE", extra_sample_sugar_base_res );
			//pose::add_variant_type_to_pose_residue( base_pose_screen, "VIRTUAL_RNA_RESIDUE", extra_sample_sugar_base_res );
			TR << base_pose_screen.annotated_sequence() << std::endl;
			TR << base_pose_screen.fold_tree() << std::endl;
		}

		// I think this should work... push apart different parts of the structure so that whatever fa_atr, fa_rep is left is
		// due to "intra-domain" interactions.
		// Crap this doesn't work when building 2 or more nucleotides.
		if ( separate_moving_residue_to_estimate_baseline_ ){
			Size jump_at_moving_suite = base_pose_screen.fold_tree().jump_nr( working_moving_res_, working_reference_res_ );
			if ( jump_at_moving_suite == 0 ){
				//				runtime_assert( std::abs( int( working_moving_res_ ) - int( working_reference_res_ ) ) == 1 );
				Size const working_moving_suite = ( working_moving_res_ > working_reference_res_ ) ? (working_moving_res_ - 1) : working_moving_res_;
				jump_at_moving_suite = make_cut_at_moving_suite( base_pose_screen, working_moving_suite );
				TR.Debug << "Made new cutpoint at suite " << working_moving_suite << " found Jump " << jump_at_moving_suite <<  std::endl;
			}
			kinematics::Jump j = base_pose_screen.jump( jump_at_moving_suite );
			j.set_translation( Vector( 1.0e4, 0.0, 0.0 ) );
			base_pose_screen.set_jump( jump_at_moving_suite, j );
		}

		( *atr_rep_screening_scorefxn_ )( base_pose_screen );

		EnergyMap const & energy_map = base_pose_screen.energies().total_energies();
		base_atr_score_ = atr_rep_screening_scorefxn_->get_weight( fa_atr ) * energy_map[ scoring::fa_atr ]; //
		base_rep_score_ = atr_rep_screening_scorefxn_->get_weight( fa_rep ) * energy_map[ scoring::fa_rep ];
		TR.Debug << "base_rep = " << base_rep_score_ << " base_atr = " << base_atr_score_ << std::endl;

		if ( output_pdb_ ) base_pose_screen.dump_pdb( "base_atr_rep_after.pdb" );

	}

	///////////////////////////////////////////////////////////////////////////////
	bool
	AtrRepScreener::check_screen( pose::Pose & current_pose_screen ){

		using namespace core::scoring;
		using namespace ObjexxFCL;

		bool close_chain = ( gap_size_ == 0 );
		if ( close_chain && is_internal_ ) return true; //Don't screen at all Mar 1, 2010

		( *atr_rep_screening_scorefxn_ )( current_pose_screen );

		EnergyMap const & energy_map = current_pose_screen.energies().total_energies();

		Real rep_score = atr_rep_screening_scorefxn_->get_weight( fa_rep ) * energy_map[scoring::fa_rep];
		Real atr_score = atr_rep_screening_scorefxn_->get_weight( fa_atr ) * energy_map[scoring::fa_atr];

		delta_rep_score_ = rep_score - base_rep_score_;
		delta_atr_score_ = atr_score - base_atr_score_;

		if ( delta_rep_score_ < (  -0.1 ) ){
			// changed from -0.01 after triggering in SWM runs -- rhiju, sep. 2013.
			std::string const message = "delta_rep_score_ = " + string_of( delta_rep_score_ ) + " rep_score = " + string_of( rep_score ) + " base_rep_score = " + string_of( base_rep_score_ );
			std::cerr << "WORKING MOVING RES " << working_moving_res_ << "  WORKING REFERENCE RES " << working_reference_res_ << "  GAP_SIZE " << gap_size_ << "  IS_PREPEND" << is_prepend_ << "  IS_INTERNAL" << is_internal_ << "  SEPARATE " << separate_moving_residue_to_estimate_baseline_ << std::endl;
			std::cerr << current_pose_screen.fold_tree() << std::endl;
			std::cerr << current_pose_screen.annotated_sequence() << std::endl;
			current_pose_screen.dump_pdb( "PROBLEM.pdb" );
			utility_exit_with_message( "delta_rep_score_ < (  -0.1 ), " + message );
		}

		if ( delta_atr_score_ > (  +0.1 ) ){
			// changed from +0.01 after triggering in SWM runs -- rhiju, sep. 2013.
			std::string const message = "delta_atr_score_ = " + string_of( delta_atr_score_ ) + " atr_score = " + string_of( atr_score ) + " base_atr_score = " + string_of( base_atr_score_ );
			std::cerr << "WORKING MOVING RES " << working_moving_res_ << "  WORKING REFERENCE RES " << working_reference_res_ << "  GAP_SIZE " << gap_size_ << "  IS_PREPEND" << is_prepend_ << "  IS_INTERNAL" << is_internal_ << "  SEPARATE " << separate_moving_residue_to_estimate_baseline_ << std::endl;
			std::cerr << current_pose_screen.fold_tree() << std::endl;
			std::cerr << current_pose_screen.annotated_sequence() << std::endl;
			current_pose_screen.dump_pdb( "PROBLEM.pdb" );
			utility_exit_with_message( "delta_atr_score_ > (  +0.1 ), " + message );
		}

		Real actual_rep_cutoff = rep_cutoff_; //default
		if ( close_chain ) {
			if ( kic_sampling_ ) {
				actual_rep_cutoff = 200.0; // KIC needs a much higher cutoff -- atoms can get really close
			} else {
				actual_rep_cutoff = 10.0; //Parin's old parameter
			}
		}
		if ( is_internal_ ) actual_rep_cutoff = 200; //Bigger moving_element..easier to crash (before May 4 used to be (close_chain && is_internal) actual_rep_cutoff=200

		bool pass_rep_screen = false;

		if ( delta_rep_score_ < actual_rep_cutoff ){
			pass_rep_screen = true;
			count_data_.good_rep_rotamer_count++;
		}

		if ( delta_atr_score_ < (  - 1 ) || close_chain ) count_data_.good_atr_rotamer_count++;

		bool pass_atr_rep_screen = false;

		if ( close_chain ){
			pass_atr_rep_screen = pass_rep_screen;
		} else if ( is_internal_ ){
			if ( delta_atr_score_ < (  - 1 ) && ( delta_rep_score_ + delta_atr_score_ ) < ( actual_rep_cutoff - rep_cutoff_ ) ) pass_atr_rep_screen = true;
		} else{
			if ( delta_atr_score_ < (  - 1 ) && ( delta_rep_score_ + delta_atr_score_ ) < 0 ) pass_atr_rep_screen = true;
		}


		if ( pass_atr_rep_screen ) {
			//	if((delta_atr_score_<(-1)) && ((delta_rep_score_+delta_atr_score_) < 200) ) { //This causes about 5times more pose to pass the screen (50,000 poses vs 10,000 poses)
			count_data_.both_count++;
			if ( verbose_ ) {
				TR.Debug << " rep = " << delta_rep_score_ << " atr = " << delta_atr_score_;
				TR.Debug << "  stack_n = " << count_data_.base_stack_count << " pair_n = " << count_data_.base_pairing_count;
				TR.Debug << "  strict_pair_n = " << count_data_.strict_base_pairing_count;
				TR.Debug << "  centroid_n = " << count_data_.pass_base_centroid_screen;
				TR.Debug << "  bin_rep_n = " << count_data_.good_bin_rep_count;
				TR.Debug << "  atr_n = " << count_data_.good_atr_rotamer_count;
				TR.Debug << "  rep_n = " << count_data_.good_rep_rotamer_count;
				TR.Debug << "  both = " << count_data_.both_count << " tot = " << count_data_.tot_rotamer_count << std::endl;
			}
			return true;
		} else {
			return false;
		}

	}

} //screener
} //rna
} //enumerate
} //stepwise
} //protocols

