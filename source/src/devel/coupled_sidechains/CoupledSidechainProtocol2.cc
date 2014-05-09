// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/sidechain_moves/SidechainMCMover.cc
/// @brief implementation of SidechainMCMover class and functions
/// @author Oliver Lange


#include <devel/coupled_sidechains/CoupledSidechainProtocol2.hh>
//#include <devel/coupled_sidechains/CoupledSidechainProtocolCreator.hh>

#include <basic/prof.hh>

// Core Headers
#include <core/types.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>

#include <core/pack/interaction_graph/SimpleInteractionGraph.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>
#include <basic/basic.hh>
//Auto Headers
#include <basic/options/option_macros.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <basic/datacache/DataMap.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>
#include <protocols/canonical_sampling/SidechainMetropolisHastingsMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/canonical_sampling/SilentTrajectoryRecorder.hh>
#include <protocols/canonical_sampling/SimulatedTempering.hh>
#include <protocols/simple_moves/sidechain_moves/PerturbRotamerSidechainMover.hh>
#include <protocols/simple_moves/sidechain_moves/JumpRotamerSidechainMover.hh>
#include <protocols/simple_moves/sidechain_moves/PerturbChiSidechainMover.hh>

#include <protocols/jd2/Job.hh>
#include <protocols/jd2/util.hh>

#include <utility/tag/Tag.hh>

// Numeric Headers
#include <numeric/angle.functions.hh>
#include <numeric/constants.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>

#include <ObjexxFCL/format.hh>
// Auto-header: duplicate removed #include <basic/Tracer.hh>

// C++ Headers
#include <ostream>
#include <sstream>
#include <fstream>
#include <utility/fixedsizearray1.hh>

#ifdef WIN_PYROSETTA
	#include <protocols/canonical_sampling/ThermodynamicObserver.hh>
#endif


using namespace core;
using namespace core::pose;
using namespace basic;
using namespace protocols::moves;
using namespace protocols;

static numeric::random::RandomGenerator Rg(38621127);
static basic::Tracer tr("devel.coupled_sidechains.CoupledSidechainProtocol");

OPT_1GRP_KEY(Integer,rotamers,observer_stride)
OPT_1GRP_KEY(Real,rotamers,unif)
OPT_1GRP_KEY(Real,rotamers,pert)
OPT_1GRP_KEY(Real,rotamers,within)
OPT_1GRP_KEY(Real,rotamers,pert_magnitude)


bool devel::coupled_sidechains::CoupledSidechainProtocol::options_registered_( false );

void devel::coupled_sidechains::CoupledSidechainProtocol::register_options() {
  using namespace basic::options;
  using namespace OptionKeys;
  if ( options_registered_ ) return;
  options_registered_ = true;
	protocols::canonical_sampling::SimulatedTempering::register_options();
	protocols::canonical_sampling::SilentTrajectoryRecorder::register_options();
	//	protocols::canonical_sampling::SilentTrajectoryRecorder::register_options();
	//protocols::canonical_sampling::SilentTrajectoryRecorder::register_options();
	OPT( score::weights );
	OPT( score::patch );
	OPT( run::n_cycles );
	OPT( packing::resfile ); // is handled directly by SidechainMoverBase
	NEW_OPT( rotamers::observer_stride," how often should the entire score-line be written to file", 100);

	NEW_OPT( rotamers::unif,"",0.1);
	NEW_OPT( rotamers::pert,"",0);
	NEW_OPT( rotamers::pert_magnitude, "deviation for perturb_chi moves in degree", 10 );
	NEW_OPT( rotamers::within,"",0);
}

namespace devel {
namespace coupled_sidechains {
/*
std::string
CoupledSidechainProtocolCreator::keyname() const {
	return CoupledSidechainProtocolCreator::mover_name();
}

protocols::moves::MoverOP
CoupledSidechainProtocolCreator::create_mover() const {
	return new CoupledSidechainProtocol;
}

std::string
CoupledSidechainProtocolCreator::mover_name() {
	return "SidechainMC";
}
*/
CoupledSidechainProtocol::CoupledSidechainProtocol()  : protocols::moves::Mover()
{
	//	set_defaults();
	init_from_options();
	setup_objects();
}

protocols::moves::MoverOP
CoupledSidechainProtocol::clone() const {
  return( protocols::moves::MoverOP( new CoupledSidechainProtocol( *this ) ) );
}

protocols::moves::MoverOP
CoupledSidechainProtocol::fresh_instance() const {
	return (protocols::moves::MoverOP( new CoupledSidechainProtocol ));
}

void
CoupledSidechainProtocol::init_from_options() {
	if ( !options_registered_ ) return;

	using namespace options;
	using namespace options::OptionKeys;
	// core::Real unif = option[ rotamers::unif ]; // Unused variable causes warning.
	core::Real pert = option[ rotamers::pert ];
	core::Real within = option[ rotamers::within ];
	ntrials_=option[ run::n_cycles ]();


	prob_withinrot_ = within;
	prob_pert_chi_ = pert;
	prob_jump_rot_ = 1.0-pert-within;
	stride_ = option[ rotamers::observer_stride ];
	pert_magnitude_ = option[ rotamers::pert_magnitude ];

}

void
CoupledSidechainProtocol::setup_objects() {

	scoring::ScoreFunctionOP scorefunction( core::scoring::getScoreFunction() );

	core::pack::task::TaskFactoryOP task_factory( new core::pack::task::TaskFactory );
	task_factory->push_back( new core::pack::task::operation::RestrictToRepacking );

	protocols::canonical_sampling::SimulatedTemperingOP tempering  = new protocols::canonical_sampling::SimulatedTempering();
	moves::MonteCarloOP mc_object = new moves::MonteCarlo( *scorefunction, 0.6 );
	tempering->set_monte_carlo( mc_object );

	sampler_ = new protocols::canonical_sampling::SidechainMetropolisHastingsMover( stride_ );
	sampler_->set_monte_carlo( mc_object );
	sampler_->set_tempering( tempering );
	sampler_->add_observer( new protocols::canonical_sampling::SilentTrajectoryRecorder );
	sampler_->set_ntrials( ntrials_ );
	if ( prob_withinrot_ > 0.0 ) {
		protocols::simple_moves::sidechain_moves::SidechainMoverBaseOP mover( new protocols::simple_moves::sidechain_moves::PerturbRotamerSidechainMover );
		//		mover->set_task_factory( task_factory );
		sampler_->add_mover( mover, prob_withinrot_ );
	}

	if ( prob_pert_chi_ > 0.0 ) {
		protocols::simple_moves::sidechain_moves::PerturbChiSidechainMoverOP mover = new protocols::simple_moves::sidechain_moves::PerturbChiSidechainMover();
		mover->set_magnitude( pert_magnitude_ );
		//		mover->set_task_factory( task_factory );
		sampler_->add_mover( mover, prob_pert_chi_ );
	}

	if ( prob_jump_rot_ > 0.0 ) {
		protocols::simple_moves::sidechain_moves::SidechainMoverBaseOP mover = new protocols::simple_moves::sidechain_moves::JumpRotamerSidechainMover();
		//		mover->set_task_factory( task_factory );
		sampler_->add_mover( mover, prob_jump_rot_ );
	}

}

/// @detailed
void
CoupledSidechainProtocol::apply(
	Pose & pose
)
{
	sampler_->apply( pose );
}

std::string
CoupledSidechainProtocol::get_name() const {
	return "CoupledSidechainProtocol";
}

void
CoupledSidechainProtocol::parse_my_tag( utility::tag::TagCOP const /*tag*/, basic::datacache::DataMap & /*data*/, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & /*pose*/) {
	/*	ntrials_ = tag->getOption<core::Size>( "ntrials", 10000 );
	set_prob_uniform( tag->getOption<core::Real>( "prob_uniform", 0.0 ) );
	set_prob_withinrot( tag->getOption<core::Real>( "prob_withinrot", 0.0 ) );
	set_prob_random_pert_current( tag->getOption<core::Real>( "prob_random_pert_current", 0.0 ) );
	core::Real between_rot = 1.0 - prob_uniform() - prob_withinrot () - prob_random_pert_current();
	set_scorefunction(
		protocols::rosetta_scripts::parse_score_function(tag, data)->clone() );
	*/
}

} // moves
} // protocols
