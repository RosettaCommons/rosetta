// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/canonical_sampling/SilentTrajectoryRecorder.cc
///
/// @brief
/// @author Oliver Lange


// Unit header or inline function header
#include <protocols/canonical_sampling/SilentTrajectoryRecorder.hh>
#include <protocols/canonical_sampling/SilentTrajectoryRecorderCreator.hh>

// Other project headers or inline function headers
#include <core/io/raw_data/ScoreStruct.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/canonical_sampling/ThermodynamicMover.hh>  // required for Windows build
#include <protocols/jd2/ScoreMap.hh>
#include <utility/tag/Tag.hh>
#include <protocols/jd2/util.hh>

// just for tracer output
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/util.hh>

// External library headers
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <iomanip>

// Operating system headers

// Forward declarations
OPT_1GRP_KEY( Integer, trajectory, score_stride )

static basic::Tracer tr( "protocols.canonical_sampling.SilentTrajectoryRecorder" );

bool protocols::canonical_sampling::SilentTrajectoryRecorder::options_registered_( false );

void protocols::canonical_sampling::SilentTrajectoryRecorder::register_options() {
  using namespace basic::options;
  using namespace OptionKeys;
  if ( options_registered_ ) return;
  options_registered_ = true;
	Parent::register_options();
	NEW_OPT( trajectory::score_stride, "write x-times score-only before a decoy is written", 1 );
}



namespace protocols {
namespace canonical_sampling {

std::string
SilentTrajectoryRecorderCreator::keyname() const {
	return SilentTrajectoryRecorderCreator::mover_name();
}

protocols::moves::MoverOP
SilentTrajectoryRecorderCreator::create_mover() const {
	return new SilentTrajectoryRecorder;
}

std::string
SilentTrajectoryRecorderCreator::mover_name() {
	return "SilentTrajectoryRecorder";
}

SilentTrajectoryRecorder::SilentTrajectoryRecorder() {
	using namespace basic::options;
  using namespace OptionKeys;
	if ( options_registered_ ) {
		score_stride_ = option[ OptionKeys::trajectory::score_stride ]();
	}
}

SilentTrajectoryRecorder::SilentTrajectoryRecorder(
	SilentTrajectoryRecorder const & other
) :	
	TrajectoryRecorder(other),
	score_stride_(other.score_stride_)
{}

protocols::moves::MoverOP
SilentTrajectoryRecorder::clone() const
{
	return new protocols::canonical_sampling::SilentTrajectoryRecorder( *this );
}

protocols::moves::MoverOP
SilentTrajectoryRecorder::fresh_instance() const
{
	return new SilentTrajectoryRecorder;
}

std::string
SilentTrajectoryRecorder::get_name() const
{
	return "SilentTrajectoryRecorder";
}

void
SilentTrajectoryRecorder::parse_my_tag(
	utility::tag::TagPtr const tag,
	protocols::moves::DataMap& data /* data */,
	protocols::filters::Filters_map const& filters /* filters */,
	protocols::moves::Movers_map const& movers /* movers */,
	core::pose::Pose const& pose /* pose */
) {
	Parent::parse_my_tag( tag, data, filters, movers, pose );
	score_stride_ = tag->getOption< core::Size >( "score_stride", 100 );
}

void
SilentTrajectoryRecorder::write_model(
	core::pose::Pose const & pose,
	protocols::canonical_sampling::MetropolisHastingsMoverCAP metropolis_hastings_mover //= 0
) {
	runtime_assert( jd2::jd2_used() );
	core::Size mc = model_count();
	jd2::output_intermediate_pose( pose, current_output_name(), mc,  ( mc % score_stride_ ) != 0 && mc > 1 ); //write always first a structure
}

void
SilentTrajectoryRecorder::initialize_simulation(
  core::pose::Pose & pose,
	protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover )
{
	if ( !cumulate_replicas() && metropolis_hastings_mover.output_name() != "" ) {
		std::ostringstream filename;
		current_output_name_ = metropolis_hastings_mover.output_name();
		tr.Info << "obtained output name " << current_output_name_ << " from MetropolisHastings Object" << std::endl;
	} else {
		current_output_name_ = file_name();
		tr.Info << "no output name obtained because " << ( cumulate_replicas() ? " cumulate-mode " : " not available " ) << std::endl;
	}
	Parent::initialize_simulation(pose, metropolis_hastings_mover);
	tr.Info << std::setprecision( 3 );
	tr.Debug << std::setprecision( 3 );
}

void
SilentTrajectoryRecorder::observe_after_metropolis(
	protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover
)
{
	protocols::moves::MonteCarlo const& mc( *(metropolis_hastings_mover.monte_carlo()) );
	Pose const& pose( mc.last_accepted_pose() );
	if (step_count() % std::max(stride(),(core::Size)500) == 0) {
		if ( tr.Info.visible() ) {
			jd2::JobOP job( jd2::get_current_job() ) ;
			tr.Info << step_count() << " E=" << pose.energies().total_energy();
			//output what is in job-object (e.g. temperature )
			for ( jd2::Job::StringRealPairs::const_iterator it( job->output_string_real_pairs_begin()), end(job->output_string_real_pairs_end()); it != end; ++it ) {
				tr.Info << " " << it->first << "=" << it->second;
			}
			tr.Info << std::endl;
		}
		mc.show_counters();
	}
	
	Parent::observe_after_metropolis(metropolis_hastings_mover);
}


} // namespace canonical_sampling
} // namespace protocols
