// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/canonical_sampling/BiasedMonteCarloMover.cc
/// @brief BiasedMonteCarlo methods implemented
/// @author


// Unit Headers
#include <protocols/canonical_sampling/BiasedMonteCarlo.hh>
//#include <protocols/canonical_sampling/BiasedMonteCarloCreator.hh>
#include <protocols/canonical_sampling/BiasEnergy.hh>
#include <protocols/moves/MonteCarloStatus.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>

#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

#ifdef BOINC_GRAPHICS
#include <protocols/boinc/boinc.hh>
#endif

static THREAD_LOCAL basic::Tracer tr( "protocols.canonical_sampling.BiasedMonteCarlo" );

namespace protocols {
namespace canonical_sampling {

using namespace core;

// std::string
// BiasedMonteCarloCreator::keyname() const {
//  return BiasedMonteCarloCreator::mover_name();
// }

// protocols::moves::MoverOP
// BiasedMonteCarloCreator::create_mover() const {
//  return new BiasedMonteCarlo;
// }

// std::string
// BiasedMonteCarloCreator::mover_name() {
//  return "BiasedMonteCarlo";
// }
BiasedMonteCarlo::BiasedMonteCarlo(
	Pose const & init_pose, // PoseCOP init_pose,
	ScoreFunction const & scorefxn, // ScoreFunctionCOP scorefxn,
	Real const temperature,
	BiasEnergyOP bias_energy
) : Parent( init_pose, scorefxn, temperature ),
	bias_energy_( bias_energy )
{
	runtime_assert( bias_energy_ != 0 );
	bias_energy_->set_temperature( temperature );
}

BiasedMonteCarlo::BiasedMonteCarlo(
	ScoreFunction const & scorefxn, // ScoreFunctionCOP scorefxn,
	Real const temperature,
	BiasEnergyOP bias_energy
) : Parent( scorefxn, temperature ),
	bias_energy_( bias_energy )
{
	runtime_assert( bias_energy_ != 0 );
	bias_energy_->set_temperature( temperature );
}

BiasedMonteCarlo::BiasedMonteCarlo( BiasedMonteCarlo const & src ) :
	Parent( src ),
	bias_energy_( src.bias_energy_ )
{
	std::cerr << "BiasedMonteCarlo copy-cstor" << std::endl;
}


void
BiasedMonteCarlo::set_temperature( Real const temp )
{
	Parent::set_temperature( temp );
	bias_energy_->set_temperature( temp );
}

bool
BiasedMonteCarlo::boltzmann(
	Pose & pose,
	std::string const & move_type, // = unk
	core::Real const proposal_density_ratio, // = 1
	core::Real const inner_score_delta_over_temperature // = 0
)
{

	// Work around a current bug in the pose observer classes..
#ifdef BOINC_GRAPHICS
	if ( get_update_boinc() ) {
		boinc::Boinc::update_graphics_current( pose );
	}
#endif

	// score the pose:
	Real const score( Parent::score_function().score( pose ) );

	// obtain the bias score
	Real const bias_score( bias_energy_->evaluate( pose ) );

	// the biased_score is used for boltzmann acceptance
	Real const biased_score( score + bias_score );

	// keep the lowest scoring pose ( based on the naked score )
	Real const last_lowest_score( lowest_score() );

	tr.Trace << "BiasedMC: (score/biased/last/lowest) " << score << " " << biased_score
		<< " " << last_accepted_score() << " " << last_lowest_score << std::endl;
	// now delegate deciscion making...
	bool const accept( Parent::boltzmann(
		biased_score,
		move_type,
		proposal_density_ratio,
		inner_score_delta_over_temperature,
		false /*don't check lowest_score */ )
	);

	// keep the lowest scoring pose ( based on the naked score )
	if ( score < lowest_score() ) { //want to keep the best scoring pose by the unbiased score
		set_lowest_score_pose( pose, score );
		set_mc_accepted( protocols::moves::MCA_accepted_score_beat_low ); //3;
#ifdef BOINC_GRAPHICS
		if ( get_update_boinc() ) {
			boinc::Boinc::update_graphics_low_energy( pose, lowest_score() );
		}
#endif
	} //MCA_accepted_score_beat_low

	// rejected ? then recover the last_accepted_pose
	if ( !accept ) {
		evaluate_convergence_checks( pose, true /*reject*/, false /* not final*/ );
		pose = last_accepted_pose();
		tr.Trace << "reject!" << std::endl;
	};

	// update the bias energy at current conformation
	// do this now, as the changed bias-energy also might influence the score of the last-accepted pose
	bias_energy_->update( pose );
	Real const updated_biased_score( pose.energies().total_energy() + bias_energy_->evaluate( pose ) );
	tr.Trace << "updated biased score: "  << score << " " << updated_biased_score << std::endl;

	if ( accept ) {
		//accepted !
		set_last_accepted_pose( pose, updated_biased_score );
	} else {
		set_last_accepted_score( updated_biased_score );
	}

#ifdef BOINC_GRAPHICS
	if ( get_update_boinc() ) {
		boinc::Boinc::update_graphics_last_accepted( pose, score );
	}
#endif

	return accept; // accept!
}

void
BiasedMonteCarlo::reset( Pose const & pose )
{
	Parent::reset( pose );
	// score the pose:
	Real const score( last_accepted_pose().energies().total_energy() );
	// obtain the bias score
	Real const bias_score( bias_energy_->evaluate( last_accepted_pose() ) );
	set_last_accepted_score( score+bias_score );
}

/// set the scorefxn,  re-scores last-accepted and lowest-score pose
void
BiasedMonteCarlo::score_function( ScoreFunction const & scorefxn )
{
	using namespace scoring;
	Parent::score_function( scorefxn );

	Real const score( last_accepted_pose().energies().total_energy() );
	// obtain the bias score
	Real const bias_score( bias_energy_->evaluate( last_accepted_pose() ) );
	set_last_accepted_score( score+bias_score );

}
} //moves
} //protocols

