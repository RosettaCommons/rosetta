// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/rna/checker/RNA_AtrRepChecker.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/rna/checker/RNA_AtrRepChecker.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/sampler/rna/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>

#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/ScoringManager.hh>

#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

static thread_local basic::Tracer TR( "protocols.stepwise.modeler.rna.checker.RNA_AtrRepChecker" );

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace checker {

	//Constructor
	RNA_AtrRepChecker::RNA_AtrRepChecker( pose::Pose const & pose,
																working_parameters::StepWiseWorkingParametersCOP & working_parameters,
																bool loose_rep_cutoff ):
		moving_res_(     working_parameters->working_moving_res() ),
		reference_res_(  working_parameters->working_reference_res() ),
		gap_size_(    working_parameters->gap_size() ),
		is_prepend_(  working_parameters->is_prepend() ),
		is_internal_(  working_parameters->is_internal() ),
		sample_both_sugar_base_rotamer_( working_parameters->sample_both_sugar_base_rotamer() ),
		separate_moving_residue_to_estimate_baseline_( true ),
		loose_rep_cutoff_( loose_rep_cutoff ),
		extra_loose_rep_cutoff_( false )
	{
		initialize_parameters();
		initialize_scorefxn();
		get_base_atr_rep_score( pose );
	}

	//Constructor
	RNA_AtrRepChecker::RNA_AtrRepChecker( pose::Pose const & pose,
																	Size const moving_res,
																	Size const reference_res,
																	Size const gap_size,
																	bool const is_internal /* = false */,
																	bool const separate_moving_residue_to_estimate_baseline, /* = true */
																	bool const sample_both_sugar_base_rotamer /* = false */
																	):
		moving_res_( moving_res    ),
		reference_res_( reference_res ),
		gap_size_( gap_size ),
		is_prepend_(  reference_res_ > moving_res_ ),
		is_internal_( is_internal ),
		sample_both_sugar_base_rotamer_( sample_both_sugar_base_rotamer ),
		separate_moving_residue_to_estimate_baseline_( separate_moving_residue_to_estimate_baseline  ),
		loose_rep_cutoff_( false ),
		extra_loose_rep_cutoff_( false )
	{
		initialize_parameters();
		initialize_scorefxn();
		get_base_atr_rep_score( pose );
	}

	//Destructor
	RNA_AtrRepChecker::~RNA_AtrRepChecker()
	{}

	///////////////////////////////////////////
	void
	RNA_AtrRepChecker::initialize_scorefxn(){
		// Bare minimum to check for contact (fa_atr) but not clash (fa_rep)
		atr_rep_screening_scorefxn_ = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction );
		atr_rep_screening_scorefxn_->set_weight( core::scoring::fa_atr , 0.23 );
		atr_rep_screening_scorefxn_->set_weight( core::scoring::fa_rep , 0.12 );
	}

	///////////////////////////////////////////
	void
	RNA_AtrRepChecker::initialize_parameters(){
		rep_cutoff_ = 4.0; // const, but can be updated below.
		base_atr_score_ = 0.0;
		base_rep_score_ = 0.0;
		delta_atr_score_ = 0.0;
		delta_rep_score_ = 0.0;
		verbose_ = false;
		output_pdb_ = false;
	}

	///////////////////////////////////////////
	void
	RNA_AtrRepChecker::get_base_atr_rep_score( core::pose::Pose const & pose ){

		using namespace core::conformation;
		using namespace core::scoring;
		using namespace core::pose;
		using namespace ObjexxFCL;

		///////////////////////////////Old_way////////////////////////////////////////////
		pose::Pose base_pose_screen = pose; //hard copy
		base_pose_screen.remove_constraints(); // floating point errors if coordinate constraints are in there.
		if ( output_pdb_ )		base_pose_screen.dump_pdb( "base_atr_rep_before.pdb" );

		Size jump_at_moving_suite = base_pose_screen.fold_tree().jump_nr( moving_res_, reference_res_ );
		Size moving_suite( 0 );
		if ( jump_at_moving_suite == 0 ){
			runtime_assert( ( moving_res_ == reference_res_ + 1 ) || ( moving_res_ == reference_res_ - 1 ) );
			moving_suite = ( moving_res_ < reference_res_ ) ? moving_res_ : reference_res_;
		}

		if ( moving_suite > 0 ){
			// suite atoms will be sampled -- not in correct conformation yet.
			pose::remove_variant_type_from_pose_residue( base_pose_screen, core::chemical::FIVE_PRIME_PHOSPHATE, moving_suite + 1 );
			pose::add_variant_type_to_pose_residue( base_pose_screen, core::chemical::VIRTUAL_PHOSPHATE, moving_suite + 1 ); //May 7...
			if ( sampler::rna::modeler_sugar_at_five_prime( base_pose_screen, moving_suite ) ) {
				add_variant_type_to_pose_residue( base_pose_screen, core::chemical::VIRTUAL_RIBOSE, moving_suite );
			}
			if ( sampler::rna::modeler_sugar_at_three_prime( base_pose_screen, moving_suite ) ){
				add_variant_type_to_pose_residue( base_pose_screen, core::chemical::VIRTUAL_RIBOSE, moving_suite+1 );
			}
		}

		// I think this should work... push apart different parts of the structure so that whatever fa_atr, fa_rep is left is
		// due to "intra-domain" interactions.
		if ( separate_moving_residue_to_estimate_baseline_ ) jump_at_moving_suite = split_pose( base_pose_screen, moving_res_, reference_res_ );

		( *atr_rep_screening_scorefxn_ )( base_pose_screen );
		EnergyMap const & energy_map = base_pose_screen.energies().total_energies();
		base_atr_score_ = atr_rep_screening_scorefxn_->get_weight( fa_atr ) * energy_map[ scoring::fa_atr ]; //
		base_rep_score_ = atr_rep_screening_scorefxn_->get_weight( fa_rep ) * energy_map[ scoring::fa_rep ];

		if ( output_pdb_ )		base_pose_screen.dump_pdb( "base_atr_rep_after.pdb" );
	}

	///////////////////////////////////////////////////////////////////////////////
	bool
	RNA_AtrRepChecker::check_screen( pose::Pose & current_pose_screen ){

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
			std::cerr << "WORKING MOVING RES " << moving_res_ << "  WORKING REFERENCE RES " << reference_res_ << "  GAP_SIZE " << gap_size_ << "  IS_PREPEND" << is_prepend_ << "  IS_INTERNAL" << is_internal_ << "  SEPARATE " << separate_moving_residue_to_estimate_baseline_ << std::endl;
			std::cerr << current_pose_screen.fold_tree() << std::endl;
			std::cerr << current_pose_screen.annotated_sequence() << std::endl;
			current_pose_screen.dump_pdb( "PROBLEM.pdb" );
			//			output_rep( current_pose_screen, "CURRENT_POSE" );
			//			get_base_atr_rep_score( current_pose_screen );
			utility_exit_with_message( "delta_rep_score_ < (  -0.1 ), " + message );
		}

		if ( delta_atr_score_ > (  +0.1 ) ){
			// changed from +0.01 after triggering in SWM runs -- rhiju, sep. 2013.
			std::string const message = "delta_atr_score_ = " + string_of( delta_atr_score_ ) + " atr_score = " + string_of( atr_score ) + " base_atr_score = " + string_of( base_atr_score_ );
			std::cerr << "WORKING MOVING RES " << moving_res_ << "  WORKING REFERENCE RES " << reference_res_ << "  GAP_SIZE " << gap_size_ << "  IS_PREPEND" << is_prepend_ << "  IS_INTERNAL" << is_internal_ << "  SEPARATE " << separate_moving_residue_to_estimate_baseline_ << std::endl;
			std::cerr << current_pose_screen.fold_tree() << std::endl;
			std::cerr << current_pose_screen.annotated_sequence() << std::endl;
			current_pose_screen.dump_pdb( "PROBLEM.pdb" );
			utility_exit_with_message( "delta_atr_score_ > (  +0.1 ), " + message );
		}

		Real actual_rep_cutoff = rep_cutoff_; //default
		if ( is_internal_ || loose_rep_cutoff_ ) {
			actual_rep_cutoff = 200.0; // KIC needs a much higher cutoff -- atoms can get really close
		} else if ( close_chain ){
			actual_rep_cutoff = 10.0; //Parin's old parameter
		}
		if ( extra_loose_rep_cutoff_ ) actual_rep_cutoff = 2000;

		bool pass_rep_screen = false;
		if ( delta_rep_score_ < actual_rep_cutoff ){
			pass_rep_screen = true;
			count_data_.good_rep_rotamer_count++;
		}

		if ( delta_atr_score_ < (  - 1 ) || close_chain ) count_data_.good_atr_rotamer_count++;

		bool pass_atr_rep_screen = false;

		if ( close_chain ){
			pass_atr_rep_screen = pass_rep_screen;
		} else if ( is_internal_ || loose_rep_cutoff_ ){
			if ( delta_atr_score_ < (  - 1.0 ) &&
					 ( delta_rep_score_ + delta_atr_score_ ) < ( actual_rep_cutoff - rep_cutoff_ ) ) pass_atr_rep_screen = true;
		} else{
			if ( delta_atr_score_ < (  - 1.0 ) &&
					 ( delta_rep_score_ + delta_atr_score_ ) < 0 ) pass_atr_rep_screen = true;
		}

		//		TR << delta_rep_score_ << " " << actual_rep_cutoff << "   " << delta_atr_score_ << " " << pass_atr_rep_screen << std::endl;

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

	///////////////////////////////////////////
	void
	RNA_AtrRepChecker::output_rep( pose::Pose const & pose, std::string const tag ){

		core::scoring::methods::EnergyMethodOptions options(atr_rep_screening_scorefxn_->energy_method_options());
		core::scoring::etable::EtableCOP etable(core::scoring::ScoringManager::get_instance()->etable( options ));
		core::scoring::etable::AnalyticEtableEvaluator eval(*etable);
		for ( Size i = 1; i <= pose.total_residue(); i++ ){
			for ( Size ii = 1; ii <= pose.residue( i ).natoms(); ii++ ){
				for ( Size j = i+1; j <= pose.total_residue(); j++ ){
					for ( Size jj = 1; jj <= pose.residue( j ).natoms(); jj++ ){
						Real atrE, repE, solE, d2;
						eval.atom_pair_energy_v( pose.residue( i ).atom( ii ),
																		 pose.residue( j ).atom( jj ),
																		 1.0, atrE, repE, solE, d2 );
						if ( repE != 0.0 ) {
							TR << tag << " " << i << " " << pose.residue( i ).atom_name( ii ) << " -- "
								 << j << " " << pose.residue( j ).atom_name( jj ) << " " << repE << std::endl;
						}
					}
				}
			}
		}
	}

} //checker
} //rna
} //modeler
} //stepwise
} //protocols

