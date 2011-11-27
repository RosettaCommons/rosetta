// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/TrajectoryRecorder.cc
///
/// @brief
/// @author


// Unit header or inline function header
#include <protocols/moves/TrajectoryRecorder.hh>

// Other project headers or inline function headers
#include <core/io/raw_data/ScoreStruct.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <protocols/moves/MetropolisHastingsMover.hh>
#include <protocols/moves/MonteCarlo.hh>
// AUTO-REMOVED #include <protocols/moves/ThermodynamicMover.hh>  // required for Windows build
#include <protocols/ScoreMap.hh>
#include <utility/tag/Tag.hh>

// just for tracer output
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/util.hh>
// External library headers

// C++ headers
#include <iomanip>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>

static basic::Tracer tr( "protocols.moves.TrajectoryRecorder" );

#include <utility/vector0.hh>
#include <utility/vector1.hh>
// Operating system headers

OPT_1GRP_KEY( Integer, trajectory, stride )
OPT_1GRP_KEY( Boolean, trajectory, cumulate )

bool protocols::moves::TrajectoryRecorder::options_registered_( false );

void protocols::moves::TrajectoryRecorder::register_options() {
  using namespace basic::options;
  using namespace OptionKeys;
  if ( options_registered_ ) return;
  options_registered_ = true;
	NEW_OPT( trajectory::stride, "how often should a snapshot be written to the trajectory", 1 );
	NEW_OPT( trajectory::cumulate, "write all decoys into the same trajectory file", false );
}


namespace protocols {
namespace moves {

TrajectoryRecorder::TrajectoryRecorder() :
	stride_(1),
	model_count_(0),
	step_count_(0),
	cumulate_( false )
{
  using namespace basic::options;
  using namespace OptionKeys;
	if ( options_registered_ ) {
		stride_ = option[ OptionKeys::trajectory::stride ]();
		cumulate_ = option[ OptionKeys::trajectory::cumulate ]();
	}
	file_name_ = "traj";
}

TrajectoryRecorder::~TrajectoryRecorder()
{}

TrajectoryRecorder::TrajectoryRecorder(
	TrajectoryRecorder const & other
) :
	ThermodynamicObserver(other),
	stride_(other.stride_),
	model_count_(other.model_count_),
	step_count_(other.step_count_),
	file_name_(other.file_name_),
	cumulate_( other.cumulate_ )
{}

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
	utility::tag::TagPtr const tag,
	protocols::moves::DataMap & /* data */,
	protocols::filters::Filters_map const & /* filters */,
	protocols::moves::Movers_map const & /* movers */,
	core::pose::Pose const & /* pose */
)
{
	stride_ = tag->getOption< core::Size >( "stride", 100 );
	file_name_ = tag->getOption< std::string >( "filename", "traj" );
	cumulate_= tag->getOption< bool > ("cumulate", false );
}

void
TrajectoryRecorder::reset( protocols::moves::MonteCarlo const & mc ) {
	model_count_ = 0;
	step_count_ = 0;
	//	write_model(mc.last_accepted_pose()); //why writing this model ? can fuck up score file
}

void
TrajectoryRecorder::update_after_boltzmann(	core::pose::Pose const & pose ) {
	++step_count_;

	if (step_count_ % stride_ == 0) {
		++model_count_;
		write_model(pose);
	}
}

void
TrajectoryRecorder::update_after_boltzmann(	protocols::moves::MonteCarlo const & mc ) {
	update_after_boltzmann(mc.last_accepted_pose());
}

void
TrajectoryRecorder::apply( core::pose::Pose & pose ) {
	update_after_boltzmann(pose);
}

void
TrajectoryRecorder::initialize_simulation(
  core::pose::Pose & pose,
	protocols::moves::MetropolisHastingsMover const & metropolis_hastings_mover )
{
	if ( !cumulate_ && metropolis_hastings_mover.output_name() != "" ) {
		std::ostringstream filename;
		current_output_name_ = metropolis_hastings_mover.output_name();
		tr.Info << "obtained output name " << current_output_name_ << " from MetropolisHastings Object" << std::endl;
	} else {
		current_output_name_ = file_name();
		tr.Info << "no output name obtained because " << ( cumulate_ ? " cumulate-mode " : " not available " ) << std::endl;
	}
	reset(*(metropolis_hastings_mover.monte_carlo()));
	tr.Info << std::setprecision( 3 );
	tr.Debug << std::setprecision( 3 );
}

void
TrajectoryRecorder::observe_after_metropolis(
	protocols::moves::MetropolisHastingsMover const & metropolis_hastings_mover
)
{
	MonteCarlo const& mc( *(metropolis_hastings_mover.monte_carlo()) );
	Pose const& pose( mc.last_accepted_pose() );
	if (step_count_ % std::max(stride_,(core::Size)500) == 0) {
		if ( tr.Info.visible() ) {
			jd2::JobOP job( jd2::get_current_job() ) ;
			tr.Info << step_count_ << " E=" << pose.energies().total_energy();
			//output what is in job-object (e.g. temperature )
			for ( jd2::Job::StringRealPairs::const_iterator it( job->output_string_real_pairs_begin()), end(job->output_string_real_pairs_end()); it != end; ++it ) {
				tr.Info << " " << it->first << "=" << it->second;
			}
			tr.Info << std::endl;
		}
		mc.show_counters();
	}
	update_after_boltzmann(*(metropolis_hastings_mover.monte_carlo()));

}

} // namespace moves
} // namespace protocols
