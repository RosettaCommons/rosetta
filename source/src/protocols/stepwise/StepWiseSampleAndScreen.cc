// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/StepWiseSampleAndScreen.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/StepWiseSampleAndScreen.hh>
#include <protocols/stepwise/screener/StepWiseScreener.hh>
#include <protocols/stepwise/screener/StepWiseScreenerType.hh>
#include <protocols/stepwise/screener/AnchorSugarScreener.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerBase.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerWithResidueAlternatives.hh>
#include <protocols/stepwise/sampler/copy_dofs/ResidueAlternativeStepWiseSamplerComb.hh>
#include <protocols/moves/CompositionMover.hh>
#include <utility/string_util.hh>
#include <ObjexxFCL/format.hh>
#include <basic/Tracer.hh>
#include <algorithm>

static thread_local basic::Tracer TR( "protocols.stepwise.StepWiseSampleAndScreen" );

using namespace core;

//////////////////////////////////////////////////////////////////////////////////////////////
//
// Unify all stepwise modeler
//
// Our previous work made it clear that information needs to be passed from early screens
//  to later ones, but this communication was not clear in our original modeler code, despite
//  extensive encapsulation.
//
// It also became clear that suite modeler (for RNA), rigid body modeler (ligands, 'floating base' for RNA),
//  protein main-chain modeler, and more could all go into one unified framework, with screens set up in a modular
//  fashion (see StepWiseScreenerType).
//
// So this is that grand framework.
//
// Note that there is still legacy code around with other classes called 'screener'. These will be deprecated or
//  renamed soon in favor of the StepWiseScreener class.
//
//  -- Rhiju, Feb 2014
//
//////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace stepwise {

	//Constructor
	StepWiseSampleAndScreen::StepWiseSampleAndScreen( sampler::StepWiseSamplerBaseOP sampler,
																										utility::vector1< screener::StepWiseScreenerOP > screeners ):
		sampler_( sampler ),
		screeners_( screeners ),
		max_ntries_( 0 ),
		num_random_samples_( 0 ),
		verbose_( false )
	{
		runtime_assert( num_screeners() > 0 );
		reset();
	}

	//Destructor
	StepWiseSampleAndScreen::~StepWiseSampleAndScreen()
	{}

	//////////////////////////////////////////////////////////////////////
	void
	StepWiseSampleAndScreen::run()
	{
		using namespace protocols::moves;
		using namespace protocols::stepwise::screener;

		if ( verbose_ ) {
			TR << "Running SampleAndScreen... " << std::endl;
			TR << std::endl; sampler_->show( TR, 0 ); TR << std::endl;
		}

		Size n( 0 );
		CompositionMoverOP update_movers( new CompositionMover ), restore_movers( new CompositionMover );
		reset();
		for ( sampler_->reset(); sampler_->not_end(); ++( *sampler_ ) ) {

			if ( sampler_->random() && ( num_tries() >= max_ntries_ || num_successes() >= num_random_samples_ ) ) break;
			update_movers->clear();
			restore_movers->clear();
			set_ok_to_increment();

			for ( n = 1; n <= num_screeners(); n++ ){

				StepWiseScreenerOP screener = screeners_[ n ];
				screener->get_update( sampler_ );
				if ( n > 1 ) screener->apply_mover( update_movers, 1, n - 1 ); //info from previous screeners.

				if ( screener->check_screen() ){
					screener->increment_count();
					screener->add_mover( update_movers, restore_movers );
					early_exit_check( n ); // for debugging.
				} else {
					screener->fast_forward( sampler_ );	break;
				}
			} // check screens

			Size const last_passed_screener = n - 1;
			for ( Size m = 2; m <= last_passed_screener; m++ ) screeners_[ m ]->apply_mover( restore_movers, m - 1, 1 );

		} // sampler

		if ( verbose_ ) output_counts();
		if ( sampler_->random() ) output_info_on_random_trials();

	}

	//////////////////////////////////////////////////////////////////////
	void
	StepWiseSampleAndScreen::reset(){
		for ( Size n = 1; n <= num_screeners(); n++ ) screeners_[ n ]->reset();
	}

	//////////////////////////////////////////////////////////////////////
	void
	StepWiseSampleAndScreen::output_counts() const{
		Size const num_digits = std::max( utility::get_num_digits( num_tries() ), Size(1) );
		for ( Size n = 1; n <= num_screeners(); n++ ) {
			TR << ObjexxFCL::format::I( num_digits, screeners_[ n ]->count() ) << " " << screeners_[ n ]->name() << std::endl;
		}
	}

	//////////////////////////////////////////////////////////////////////
	Size const &
	StepWiseSampleAndScreen::num_tries() const{
		return screeners_[1]->count();
	}

	//////////////////////////////////////////////////////////////////////
	Size const &
	StepWiseSampleAndScreen::num_successes() const{
		return screeners_[ num_screeners() ]->count();
	}

	//////////////////////////////////////////////////////////////////////
	Size
	StepWiseSampleAndScreen::num_screeners() const{
		return screeners_.size();
	}


	//////////////////////////////////////////////////////////////////////
	void
	StepWiseSampleAndScreen::output_info_on_random_trials() const{
		if ( !sampler_->random() ) return;
		TR << "Number of tries: " << num_tries()  << ". Number of successes: " << num_successes() <<  std::endl;
		TR.Debug << "Was shooting for max_tries: " << max_ntries_ << ". num_random_samples: " << num_random_samples_ << std::endl;
	}

	///////////////////////////////////////////////////////////////////////////
	// kind of a hack, but needed right now for comparison to prior work.
	// did not want to keep track of some 'inner loops' (sugar modeler),
	// to be closer in line to classic SWA sampler code.
	void
	StepWiseSampleAndScreen::set_ok_to_increment(){

		bool ok_to_increment_screeners( true );

		using namespace protocols::stepwise::sampler;
		using namespace protocols::stepwise::sampler::rigid_body;
		if ( sampler_->type() == RIGID_BODY_WITH_RESIDUE_ALTERNATIVES ){
			RigidBodyStepWiseSamplerWithResidueAlternatives & rigid_body_rotamer_with_residue_alternatives = *( static_cast< RigidBodyStepWiseSamplerWithResidueAlternatives * >( sampler_.get() ) );
			ok_to_increment_screeners = ( rigid_body_rotamer_with_residue_alternatives.residue_alternatives_rotamer()->id() == 1 );
		}
		if ( ok_to_increment_screeners ){
			for ( Size n = 1; n <= num_screeners(); n++ ) screeners_[ n ]->set_ok_to_increment( true );
		}
	}


	/////////////////////////////////////////////////////////////////
	void
	StepWiseSampleAndScreen::early_exit_check( Size const n ) {

		return;

		/////////////////////////////////////////
		// following just for debugging.
		/////////////////////////////////////////
		if ( n < num_screeners() ) return;

		bool early_exit( false );
		Size anchor_sugar_solution_number( 0 );
		using namespace screener;
		for ( Size j = 1; j <= num_screeners(); j++ ){
			screener::StepWiseScreenerOP screener = screeners_[ j ];
			if ( screeners_[j]->type() == ANCHOR_SUGAR ){
				AnchorSugarScreener & anchor_sugar_screener = *( static_cast< AnchorSugarScreener * >( screener.get() ) );
				anchor_sugar_solution_number = anchor_sugar_screener.anchor_sugar_solution_number();
				break;
			}
		}


		using namespace protocols::stepwise::sampler;
		using namespace protocols::stepwise::sampler::rigid_body;
		screener::StepWiseScreenerOP screener = screeners_[ n ];
		if ( sampler_->type() == RIGID_BODY_WITH_RESIDUE_ALTERNATIVES ){
			RigidBodyStepWiseSamplerWithResidueAlternatives & rigid_body_rotamer_with_residue_alternatives = *( static_cast< RigidBodyStepWiseSamplerWithResidueAlternatives * >( sampler_.get() ) );
			TR << "Rigid body ID " << rigid_body_rotamer_with_residue_alternatives.rigid_body_rotamer()->id() << ": ";
			TR << " ID overall " << rigid_body_rotamer_with_residue_alternatives.residue_alternatives_rotamer()->id() << "; ";
			TR << " ID at 3 " << rigid_body_rotamer_with_residue_alternatives.residue_alternatives_rotamer()->id_for_resnum( 3 ) << "; ";
			bool anchor_sugar_screener_legacy = ( anchor_sugar_solution_number  > 0 );
			if ( !anchor_sugar_screener_legacy ){
				anchor_sugar_solution_number = rigid_body_rotamer_with_residue_alternatives.residue_alternatives_rotamer()->id_for_resnum( 5 );
			}
			TR << " ID at 5 " << anchor_sugar_solution_number;
			if ( anchor_sugar_solution_number ==  2 )	early_exit = true;
			TR << " [legacy: " << anchor_sugar_screener_legacy << "]" << std::endl;
		}


		if ( early_exit ){
			std::cout << screener->name() <<  " NTRIES " << num_tries() << std::endl;
			output_counts();
			exit( 0 );
		}

	}


} //stepwise
} //protocols
