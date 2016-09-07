// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/canonical_sampling/TrajectoryRecorder.cc
///
/// @brief
/// @author


// Unit header or inline function header
#include <protocols/canonical_sampling/TrajectoryRecorder.hh>

// Other project headers or inline function headers
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <utility/tag/Tag.hh>

// External library headers

// C++ headers
#include <iomanip>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer tr( "protocols.canonical_sampling.TrajectoryRecorder" );

#include <utility/vector0.hh>
#include <utility/vector1.hh>
// Operating system headers

OPT_1GRP_KEY( Integer, trajectory, stride )
OPT_1GRP_KEY( Integer, trajectory, cache_limit )
OPT_1GRP_KEY( Boolean, trajectory, cumulate_jobs )
OPT_1GRP_KEY( Boolean, trajectory, cumulate_replicas )

bool protocols::canonical_sampling::TrajectoryRecorder::options_registered_( false );

void protocols::canonical_sampling::TrajectoryRecorder::register_options() {
	using namespace basic::options;
	using namespace OptionKeys;
	if ( options_registered_ ) return;
	options_registered_ = true;
	NEW_OPT( trajectory::stride, "how often should a snapshot be written to the trajectory", 1 );
	NEW_OPT( trajectory::cache_limit, "the maximum number of poses to cache before performing IO", 500 );
	NEW_OPT( trajectory::cumulate_jobs, "write structures from different jobs into the same trajectory file", false );
	NEW_OPT( trajectory::cumulate_replicas, "write structures from different replicas the same trajectory file", false );
}


namespace protocols {
namespace canonical_sampling {

TrajectoryRecorder::TrajectoryRecorder() :
	stride_(1),
	model_count_(0),
	step_count_(0),
	cache_limit_(500),
	cumulate_jobs_( false ),
	cumulate_replicas_( false )
{
	using namespace basic::options;
	using namespace OptionKeys;
	if ( options_registered_ ) {
		stride_ = option[ OptionKeys::trajectory::stride ]();
		cache_limit_ = option[ OptionKeys::trajectory::cache_limit ]();
		cumulate_jobs_ = option[ OptionKeys::trajectory::cumulate_jobs ]();
		cumulate_replicas_ = option[ OptionKeys::trajectory::cumulate_replicas ]();
	}
	file_name_ = "traj";
}

TrajectoryRecorder::~TrajectoryRecorder() = default;

TrajectoryRecorder::TrajectoryRecorder( TrajectoryRecorder const & ) = default;

TrajectoryRecorder&
TrajectoryRecorder::operator=( TrajectoryRecorder const & /* other */ )
{
	// assignment not allowed
	runtime_assert(false);
	return *this;
}

std::string
TrajectoryRecorder::get_name() const
{
	return "TrajectoryRecorder";
}

void
TrajectoryRecorder::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /* data */,
	protocols::filters::Filters_map const & /* filters */,
	protocols::moves::Movers_map const & /* movers */,
	core::pose::Pose const & /* pose */
)
{
	stride_ = tag->getOption< core::Size >( "stride", 100 );
	file_name_ = tag->getOption< std::string >( "filename", file_name_ );
	cache_limit_ = tag->getOption< core::Size >( "cache_limit", 500 );
	cumulate_jobs_= tag->getOption< bool > ("cumulate_jobs", false );
	cumulate_replicas_= tag->getOption< bool > ("cumulate_replicas", false );
}

void
TrajectoryRecorder::reset(
	protocols::moves::MonteCarlo const &,
	protocols::canonical_sampling::MetropolisHastingsMover const * //=0
) {
	model_count_ = 0;
	step_count_ = 0;
}

void
TrajectoryRecorder::update_after_boltzmann(
	core::pose::Pose const & pose,
	protocols::canonical_sampling::MetropolisHastingsMover const * metropolis_hastings_mover //= 0
) {
	++step_count_;

	if ( step_count_ % stride_ == 0 ) {
		++model_count_;
		write_model(pose, metropolis_hastings_mover);
	}
}

void
TrajectoryRecorder::update_after_boltzmann( protocols::moves::MonteCarlo const & mc ) {
	update_after_boltzmann(mc.last_accepted_pose());
}

void
TrajectoryRecorder::apply( core::pose::Pose & pose ) {
	update_after_boltzmann(pose);
}

void
TrajectoryRecorder::initialize_simulation(
	core::pose::Pose &,
	protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover,
	core::Size cycle //default=0; non-zero if trajectory is restarted
) {
	reset(
		*(metropolis_hastings_mover.monte_carlo()),
		&metropolis_hastings_mover
	);
	if ( cycle != 0 ) {
		step_count_ = cycle;
		model_count_ = step_count_ / stride();
	}
}

void
TrajectoryRecorder::observe_after_metropolis(
	protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover
)
{
	update_after_boltzmann(
		metropolis_hastings_mover.monte_carlo()->last_accepted_pose(),
		&metropolis_hastings_mover
	);
}

} // namespace canonical_sampling
} // namespace protocols
