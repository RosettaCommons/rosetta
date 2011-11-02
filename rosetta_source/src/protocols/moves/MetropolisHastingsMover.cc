// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/moves/MetropolisHastingsMover.cc
/// @brief MetropolisHastingsMover methods implemented
/// @author


// Unit Headers
#include <protocols/moves/MetropolisHastingsMover.hh>
#include <protocols/moves/MetropolisHastingsMoverCreator.hh>


// protocols headers
#include <protocols/moves/BackrubMover.hh>
#include <protocols/moves/DataMap.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/moves/SidechainMover.hh>
#include <protocols/moves/SidechainMCMover.hh>
#include <protocols/moves/BackboneMover.hh>
#include <protocols/moves/ThermodynamicMover.hh>
#include <protocols/moves/ThermodynamicObserver.hh>
#include <protocols/rosetta_scripts/util.hh>

// core headers
#include <core/kinematics/MoveMap.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh> // REQUIRED FOR WINDOWS

// numeric headers
#include <numeric/random/random.hh>

// utility headers
#include <utility/file/file_sys_util.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tag/Tag.hh>

#include <core/chemical/ResidueType.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


// C++ Headers

using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR( "protocols.moves.MetropolisHastingsMover" );
static numeric::random::RandomGenerator RG(638767547);

namespace protocols {
namespace moves {

std::string
MetropolisHastingsMoverCreator::keyname() const {
	return MetropolisHastingsMoverCreator::mover_name();
}

protocols::moves::MoverOP
MetropolisHastingsMoverCreator::create_mover() const {
	return new MetropolisHastingsMover;
}

std::string
MetropolisHastingsMoverCreator::mover_name() {
	return "MetropolisHastings";
}

MetropolisHastingsMover::MetropolisHastingsMover() :
	monte_carlo_(0),
	ntrials_(1000),
	trial_(0)
{}

MetropolisHastingsMover::MetropolisHastingsMover(
	MetropolisHastingsMover const & metropolis_hastings_mover
) :
	//utility::pointer::ReferenceCount(),
	Mover(metropolis_hastings_mover),
	monte_carlo_(new protocols::moves::MonteCarlo(*metropolis_hastings_mover.monte_carlo_)),
	ntrials_(metropolis_hastings_mover.ntrials_),
	trial_(metropolis_hastings_mover.trial_),
	weighted_sampler_(metropolis_hastings_mover.weighted_sampler_),
	output_name_(metropolis_hastings_mover.output_name_)
{
	for (core::Size i = 1; i <= metropolis_hastings_mover.movers_.size(); ++i) {
		movers_.push_back(reinterpret_cast<protocols::moves::ThermodynamicMover *>(metropolis_hastings_mover.movers_[i]->clone()()));
	}

	for (core::Size i = 1; i <= metropolis_hastings_mover.observers_.size(); ++i) {
		observers_.push_back(reinterpret_cast<protocols::moves::ThermodynamicObserver *>(metropolis_hastings_mover.observers_[i]->clone()()));
	}
}

MetropolisHastingsMover::~MetropolisHastingsMover(){}

void
MetropolisHastingsMover::apply( core::pose::Pose & pose )
{
	bool output_name_from_job_distributor(false);

	if (output_name() == "") {
		set_output_name(protocols::jd2::JobDistributor::get_instance()->current_output_name());
		output_name_from_job_distributor = true;
	}

	for (core::Size i = 1; i <= movers_.size(); ++i) {
		TR << "Initializing " << movers_[i]->get_name() << std::endl;
		movers_[i]->set_metropolis_hastings_mover(this);
		movers_[i]->initialize_simulation(pose, *this);
	}

	monte_carlo_->reset(pose);
	monte_carlo_->reset_counters();

	for (core::Size i = 1; i <= observers_.size(); ++i) {
		TR << "Initializing " << observers_[i]->get_name() << std::endl;
		observers_[i]->initialize_simulation(pose, *this);
	}

	TR << "Initial Score:" << std::endl;
	monte_carlo_->score_function().show(TR, pose);

	TR << "Running " << ntrials_ << " trials..." << std::endl;

	for (trial_ = 1; trial_ <= ntrials_; ++trial_) {

		protocols::moves::ThermodynamicMoverOP mover(random_mover());
		mover->apply(pose);
		monte_carlo_->boltzmann(pose, mover->type(), mover->last_proposal_density_ratio(), mover->last_inner_score_temperature_delta());
		mover->observe_after_metropolis(*this);
		
		for (core::Size i = 1; i <= observers_.size(); ++i) {
			observers_[i]->observe_after_metropolis(*this);
		}
	}

	for (core::Size i = 1; i <= movers_.size(); ++i) {
		TR << "Finalizing " << movers_[i]->get_name() << std::endl;
		movers_[i]->finalize_simulation(pose, *this);
	}

	for (core::Size i = 1; i <= observers_.size(); ++i) {
		TR << "Finalizing " << observers_[i]->get_name() << std::endl;
		observers_[i]->finalize_simulation(pose, *this);
	}

	TR << "Final Score:" << std::endl;
	monte_carlo_->score_function().show(TR, pose);
	TR.flush();

	monte_carlo_->show_counters();

	if (output_name_from_job_distributor) set_output_name("");
}

std::string
MetropolisHastingsMover::get_name() const
{
	return "MetropolisHastingsMover";
}

protocols::moves::MoverOP
MetropolisHastingsMover::clone() const
{
	return new MetropolisHastingsMover(*this);
}

protocols::moves::MoverOP
MetropolisHastingsMover::fresh_instance() const
{
	return new MetropolisHastingsMover;
}

bool
MetropolisHastingsMover::reinitialize_for_each_job() const { return false; }

bool
MetropolisHastingsMover::reinitialize_for_new_input() const { return false; }

void
MetropolisHastingsMover::parse_my_tag(
	utility::tag::TagPtr const tag,
	protocols::moves::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
)
{
	core::Real const temperature(tag->getOption< core::Real >( "temperature", 0.6 ) );
	ntrials_ = tag->getOption< core::Size >( "trials", ntrials_ );

	core::scoring::ScoreFunctionOP score_fxn(protocols::rosetta_scripts::parse_score_function(tag, data));

	monte_carlo_ = new protocols::moves::MonteCarlo(*score_fxn, temperature);

	protocols::moves::MoverFactory *mover_factory(protocols::moves::MoverFactory::get_instance());

	utility::vector0< utility::tag::TagPtr > const subtags( tag->getTags() );

	for( utility::vector0< utility::tag::TagPtr >::const_iterator subtag_it = subtags.begin(); subtag_it != subtags.end(); ++subtag_it ) {

		TagPtr const subtag = *subtag_it;

		protocols::moves::MoverOP mover;

		if (subtag->getName() == "Add") {

			std::string mover_name = subtag->getOption<std::string>( "mover_name", "null" );
			Movers_map::const_iterator mover_iter( movers.find( mover_name ) );
			if ( mover_iter == movers.end() ) {
				TR << "Mover not found for XML tag:\n" << subtag << std::endl;
				utility_exit();
			}
			mover = mover_iter->second;

		} else {

			mover = mover_factory->newMover(subtag, data, filters, movers, pose);
		}

		protocols::moves::ThermodynamicMoverCOP thermomover( dynamic_cast<protocols::moves::ThermodynamicMover *>(mover()) );
		protocols::moves::ThermodynamicObserverCOP thermoobserver( dynamic_cast<protocols::moves::ThermodynamicObserver *>(mover()) );

		if (thermomover) {
			core::Real const weight( subtag->getOption< core::Real >( "sampling_weight", 1 ) );
			add_mover( reinterpret_cast<protocols::moves::ThermodynamicMover *>(thermomover->clone()()), weight);
		} else if (thermoobserver) {
			add_observer( reinterpret_cast<protocols::moves::ThermodynamicObserver *>(thermoobserver->clone()()));
		} else {
			TR << "Mover is not a ThermodynamicMover or ThermodynamicObserver for XML tag:\n" << subtag << std::endl;
			utility_exit();
		}
	}
}

protocols::moves::MonteCarloCOP
MetropolisHastingsMover::monte_carlo() const
{
	return monte_carlo_;
}

void
MetropolisHastingsMover::set_monte_carlo(
	protocols::moves::MonteCarloOP monte_carlo
)
{
	monte_carlo_ = monte_carlo;
}

core::Size
MetropolisHastingsMover::ntrials() const
{
	return ntrials_;
}

void
MetropolisHastingsMover::set_ntrials(
	core::Size ntrials
)
{
	ntrials_ = ntrials;
}

bool
MetropolisHastingsMover::finished() const
{
	return trial_ > ntrials_;
}

std::string const &
MetropolisHastingsMover::output_name() const
{
	return output_name_;
}

void
MetropolisHastingsMover::set_output_name(
	std::string const & output_name
)
{
	output_name_ = output_name;
}

protocols::moves::ThermodynamicMoverOP
MetropolisHastingsMover::random_mover()
{
	return movers_[weighted_sampler_.random_sample(RG)];
}

void
MetropolisHastingsMover::add_mover(
	protocols::moves::ThermodynamicMoverOP mover,
	core::Real weight
)
{
	mover->set_preserve_detailed_balance(true);
	movers_.push_back(mover);
	weighted_sampler_.add_weight(weight);
}

void
MetropolisHastingsMover::add_backrub_mover(
	core::Real weight
)
{
	if (!weight) return;

	protocols::moves::BackrubMoverOP backrub_mover(new protocols::moves::BackrubMover);
	backrub_mover->branchopt().read_database();

	add_mover(backrub_mover, weight);
}

void
MetropolisHastingsMover::add_small_mover(
	core::Real weight
)
{
	if (!weight) return;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	protocols::moves::SmallMoverOP small_mover(new protocols::moves::SmallMover);
	small_mover->nmoves(1);
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	if (utility::file::file_exists(option[ in::file::movemap ])) {
		movemap->init_from_file(option[ in::file::movemap ]);
	} else {
		movemap->set_bb(true);
	}
	small_mover->movemap(movemap);

	add_mover(small_mover, weight);
}

void
MetropolisHastingsMover::add_shear_mover(
	core::Real weight
)
{
	if (!weight) return;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	protocols::moves::ShearMoverOP shear_mover(new protocols::moves::ShearMover);
	shear_mover->nmoves(1);
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	if (utility::file::file_exists(option[ in::file::movemap ])) {
		movemap->init_from_file(option[ in::file::movemap ]);
	} else {
		movemap->set_bb(true);
	}
	shear_mover->movemap(movemap);

	add_mover(shear_mover, weight);
}

void
MetropolisHastingsMover::add_sidechain_mover(
	core::Real weight,
	core::Real prob_uniform,
	core::Real prob_withinrot,
	bool preserve_cbeta
)
{
	if (!weight) return;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::pack::task::TaskFactoryOP main_task_factory = new core::pack::task::TaskFactory;
	main_task_factory->push_back( new core::pack::task::operation::InitializeFromCommandline );
	if ( option[ packing::resfile ].user() ) {
		main_task_factory->push_back( new core::pack::task::operation::ReadResfile );
	} else {
		main_task_factory->push_back( new core::pack::task::operation::RestrictToRepacking );
	}
	if (preserve_cbeta) main_task_factory->push_back( new core::pack::task::operation::PreserveCBeta );

	protocols::moves::SidechainMoverOP sidechain_mover(new protocols::moves::SidechainMover);
	sidechain_mover->set_task_factory(main_task_factory);
	sidechain_mover->set_prob_uniform(prob_uniform);
	sidechain_mover->set_prob_withinrot(prob_withinrot);

	add_mover(sidechain_mover, weight);
}

void
MetropolisHastingsMover::add_sidechain_mc_mover(
	core::Real weight,
	core::Real prob_uniform,
	core::Real prob_withinrot,
	bool preserve_cbeta,
	core::Size ntrials
)
{
	if (!weight) return;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::pack::task::TaskFactoryOP main_task_factory = new core::pack::task::TaskFactory;
	main_task_factory->push_back( new core::pack::task::operation::InitializeFromCommandline );
	if ( option[ packing::resfile ].user() ) {
		main_task_factory->push_back( new core::pack::task::operation::ReadResfile );
	} else {
		main_task_factory->push_back( new core::pack::task::operation::RestrictToRepacking );
	}
	if (preserve_cbeta) main_task_factory->push_back( new core::pack::task::operation::PreserveCBeta );

	protocols::moves::SidechainMCMoverOP sidechain_mc_mover(new protocols::moves::SidechainMCMover);
	sidechain_mc_mover->set_task_factory(main_task_factory);
	sidechain_mc_mover->set_prob_uniform(prob_uniform);
	sidechain_mc_mover->set_prob_withinrot(prob_withinrot);
	sidechain_mc_mover->setup(monte_carlo_->score_function().clone());
	sidechain_mc_mover->set_ntrials(ntrials);

	add_mover(sidechain_mc_mover, weight);
}

void
MetropolisHastingsMover::add_observer(
	protocols::moves::ThermodynamicObserverOP observer
)
{
	observers_.push_back(observer);
}

} //moves
} //protocols

