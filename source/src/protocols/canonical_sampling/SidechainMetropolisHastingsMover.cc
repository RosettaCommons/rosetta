// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/canonical_sampling/SidechainMetropolisHastingsMover.cc
/// @brief SidechainMetropolisHastingsMover methods implemented
/// @author


// Unit Headers
#include <protocols/canonical_sampling/SidechainMetropolisHastingsMover.hh>
#include <protocols/canonical_sampling/SidechainMetropolisHastingsMoverCreator.hh>

#include <protocols/simple_moves/sidechain_moves/SidechainMoverBase.hh>

// protocols headers
#include <protocols/backrub/BackrubMover.hh>
#include <basic/datacache/DataMap.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/Job.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMCMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/canonical_sampling/ThermodynamicMover.hh>
#include <protocols/canonical_sampling/ThermodynamicObserver.hh>
#include <protocols/rosetta_scripts/util.hh>

// core headers
#include <core/kinematics/MoveMap.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <core/pose/Pose.hh>
#include <core/pack/interaction_graph/SimpleInteractionGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/conformation/Residue.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>


// numeric headers
#include <numeric/random/random.hh>

// utility headers
#include <utility/file/file_sys_util.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tag/Tag.hh>

// C++ Headers

using basic::T;
using basic::Error;
using basic::Warning;

static thread_local basic::Tracer tr( "protocols.canonical_sampling.SidechainSidechainMetropolisHastingsMover" );

namespace protocols {
namespace canonical_sampling {

using namespace core;
using namespace scoring;

std::string
SidechainMetropolisHastingsMoverCreator::keyname() const {
	return SidechainMetropolisHastingsMoverCreator::mover_name();
}

protocols::moves::MoverOP
SidechainMetropolisHastingsMoverCreator::create_mover() const {
	return new SidechainMetropolisHastingsMover;
}

std::string
SidechainMetropolisHastingsMoverCreator::mover_name() {
	return "SidechainMetropolisHastings";
}

SidechainMetropolisHastingsMover::SidechainMetropolisHastingsMover() :
	stride_( 1000 )
{}

SidechainMetropolisHastingsMover::SidechainMetropolisHastingsMover( core::Size stride ) :
	stride_( stride )
{}

SidechainMetropolisHastingsMover::SidechainMetropolisHastingsMover(
	SidechainMetropolisHastingsMover const & other
) :
	protocols::canonical_sampling::MetropolisHastingsMover(other),
	stride_( other.stride_ )
{}

SidechainMetropolisHastingsMover::~SidechainMetropolisHastingsMover(){}

bool
SidechainMetropolisHastingsMover::pass_metropolis(core::Real delta_energy, core::Real last_proposal_density_ratio ) const {
	core::Real boltz_factor = delta_energy / monte_carlo()->temperature();
	if ( tr.Trace.visible() ) {
		tr.Trace << " temperature: " << monte_carlo()->temperature()
						 << " deltaE= " << delta_energy
						 << " boltzman=" << boltz_factor
						 << " lpd= " << last_proposal_density_ratio << std::endl;
	}
	core::Real probability = std::exp( std::min( 40.0, std::max( -40.0, boltz_factor ))) *  last_proposal_density_ratio ;
	if ( probability < 1 && numeric::random::rg().uniform() >= probability ) {
		return false;
	} else {
		return true;
	}
}

void
SidechainMetropolisHastingsMover::apply( core::pose::Pose & pose )
{
	prepare_simulation( pose );

	scoring::ScoreFunction const& sfxn = monte_carlo()->score_function();
	pack::interaction_graph::SimpleInteractionGraphOP ig;
	ig = new pack::interaction_graph::SimpleInteractionGraph(); //commented out debug
	ig->set_scorefunction( sfxn );

	utility::vector1< Real > new_chi;
	Real current_energy = sfxn(pose);

	//	std::string const traj_file_tag( jd2::current_output_name() );
	//  counters_.reset();

	jd2::JobOP job;
	if ( jd2::jd2_used() ) {
		job = jd2::get_current_job();
	}

	//ek for fast sidechain sampling and internal mc trials
	utility::vector1< conformation::ResidueOP > current;
	//	utility::vector1< conformation::ResidueOP > previous;
	//	utility::vector1< pack::dunbrack::ChiVector > chi_vectors;
	//	utility::vector1< pack::dunbrack::RotVector > rot_vectors;


	current.resize(pose.total_residue());
	//	previous.resize(pose.total_residue());

	//	rot_vectors.resize( pose.total_residue() );
	//	chi_vectors.resize( pose.total_residue() );

	for ( core::Size i = 1; i <= pose.total_residue(); i++ ){
		current[ i ] = new core::conformation::Residue( pose.residue( i ) );
	}


	runtime_assert( ig != 0 );

	ig->initialize( pose );
	Real last_accepted_prop_density( 1.0 );
	Real last_accepted_dE( 0.0 );
	for ( Size ct = 1; ct <= ntrials(); ct++) {
		protocols::simple_moves::sidechain_moves::SidechainMoverBaseOP move = dynamic_cast< protocols::simple_moves::sidechain_moves::SidechainMoverBase* >( random_mover().get() );
		runtime_assert( move != 0 ); //fow now only Sidechain Movers...

		Size resid = move->suggest_residue_number( pose );
		conformation::ResidueOP new_state( new conformation::Residue( pose.residue( resid ) ) );
		new_state = move->make_move( new_state );
		set_last_move( move );

		Real delta_energy = ig->consider_substitution( resid, new_state );
		if ( pass_metropolis( delta_energy, move->last_proposal_density_ratio() ) ) { //ek
			ig->commit_change( resid );
			current_energy -= delta_energy;
			current[ resid ] = new_state;
			set_last_accepted( true );
			last_accepted_prop_density = move->last_proposal_density_ratio();
			last_accepted_dE = delta_energy;
		} else { //rejected metropolis criterion
			ig->reject_change( resid );
			set_last_accepted( false );
		}

		tempering()->temperature_move( current_energy );
		move->observe_after_metropolis( *this );

		Size model_count( output_count( ct ) );
		if ( model_count ) {
			for( Size res_i = 1; res_i <= current.size(); res_i++ ){
				pose.replace_residue( res_i, (*current[ res_i ]), true );
			}
			core::Real const score( sfxn( pose ) );
			if ( std::abs( score-current_energy ) > 1 ) { //threshold 0.1 gives a couple warnings -- but it never drifts apart
				tr.Warning << "Energy mismatch!!! score=" << score << " ig->energy " << current_energy << std::endl;
			}
			nonconst_monte_carlo().set_last_accepted_pose( pose );
			if ( score < nonconst_monte_carlo().lowest_score() ) {
				nonconst_monte_carlo().set_lowest_score_pose( pose );
			}
			if ( job ) {
				job->add_string_real_pair( "prop_density", move->last_proposal_density_ratio() );
				job->add_string_real_pair( "prop_density_accept", last_accepted_prop_density );
				job->add_string_real_pair( "move_dE", delta_energy );
				job->add_string_real_pair( "move_dE_accept", last_accepted_dE );
				job->add_string_string_pair( "move_type", move->type() );
			}
		}

		for (Size i = 1; i <= observers().size(); ++i) {
			if ( observers()[ i ]->requires_pose() && !model_count ) continue;
			observers()[i]->observe_after_metropolis(*this);
		}

	} //for ntrials

	wind_down_simulation( pose );
}

core::Size
SidechainMetropolisHastingsMover::output_count( Size ct ) const {
	if ( ct % stride_ == 0 ) {
		return ct / stride_;
	}	else return 0;
}

std::string
SidechainMetropolisHastingsMover::get_name() const
{
	return "SidechainMetropolisHastingsMover";
}

protocols::moves::MoverOP
SidechainMetropolisHastingsMover::clone() const
{
	return new protocols::canonical_sampling::SidechainMetropolisHastingsMover(*this);
}

protocols::moves::MoverOP
SidechainMetropolisHastingsMover::fresh_instance() const
{
	return new SidechainMetropolisHastingsMover;
}

void
SidechainMetropolisHastingsMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	pose::Pose const & pose
) {
	stride_ = tag->getOption< Size >( "stride", stride_ );
	Parent::parse_my_tag( tag, data, filters, movers, pose );
}


} //moves
} //protocols

