// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/canonical_sampling/SilentTrajectoryRecorder.cc
///
/// @brief
/// @author Oliver Lange


// Unit header or inline function header
#include <protocols/canonical_sampling/SilentTrajectoryRecorder.hh>
#include <protocols/canonical_sampling/SilentTrajectoryRecorderCreator.hh>

// Other project headers or inline function headers
#include <core/io/silent/SilentFileData.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/canonical_sampling/ThermodynamicMover.hh>  // required for Windows build
#include <protocols/jd2/ScoreMap.hh>
#include <utility/tag/Tag.hh>
#include <protocols/jd2/util.hh>
#include <utility/string_util.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <protocols/jd2/SilentFileJobOutputter.hh>
#include <protocols/jd2/JobDistributor.hh>

// just for tracer output
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/util.hh>

// External library headers
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>

// C++ headers
#include <iomanip>

// Operating system headers

// Forward declarations
OPT_1GRP_KEY( Integer, trajectory, score_stride )

static THREAD_LOCAL basic::Tracer tr( "protocols.canonical_sampling.SilentTrajectoryRecorder" );

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
using namespace core;

std::string
SilentTrajectoryRecorderCreator::keyname() const {
	return SilentTrajectoryRecorderCreator::mover_name();
}

protocols::moves::MoverOP
SilentTrajectoryRecorderCreator::create_mover() const {
	return protocols::moves::MoverOP( new SilentTrajectoryRecorder );
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
	return protocols::moves::MoverOP( new protocols::canonical_sampling::SilentTrajectoryRecorder( *this ) );
}

protocols::moves::MoverOP
SilentTrajectoryRecorder::fresh_instance() const
{
	return protocols::moves::MoverOP( new SilentTrajectoryRecorder );
}

std::string
SilentTrajectoryRecorder::get_name() const
{
	return "SilentTrajectoryRecorder";
}

void
SilentTrajectoryRecorder::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& data /* data */,
	protocols::filters::Filters_map const& filters /* filters */,
	protocols::moves::Movers_map const& movers /* movers */,
	core::pose::Pose const& pose /* pose */
) {
	Parent::parse_my_tag( tag, data, filters, movers, pose );
	score_stride_ = tag->getOption< core::Size >( "score_stride", 100 );
	runtime_assert( jd2::jd2_used() );
	job_outputter_ = jd2::JobDistributor::get_instance()->job_outputter();

	std::string silent_struct_type = tag->getOption< std::string >( "silent_struct_type", "any" );
	if ( silent_struct_type != "any" ) {
		jd2::SilentFileJobOutputterOP silent_job_outputter( new jd2::SilentFileJobOutputter() );
		silent_job_outputter->set_forced_silent_struct_type( silent_struct_type );
		silent_job_outputter->set_write_separate_scorefile( tag->getOption< bool >( "write_extra_scores", false ) );
		job_outputter_ = silent_job_outputter;
	}
}

void
SilentTrajectoryRecorder::write_model(
	core::pose::Pose const & pose,
	protocols::canonical_sampling::MetropolisHastingsMover const * metropolis_hastings_mover //= 0
) {
	runtime_assert( jd2::jd2_used() );
	std::string filename( metropolis_hastings_mover ? metropolis_hastings_mover->output_file_name(file_name(), cumulate_jobs(), cumulate_replicas()) : file_name() );
	core::Size mc = model_count();
	tr.Debug << "write model " << filename << " count: " << mc << std::endl;
	job_outputter_->other_pose( jd2::get_current_job(),  pose, filename, mc,  ( mc % score_stride_ ) != 0 && mc > 1 );
}

bool
SilentTrajectoryRecorder::restart_simulation(
	core::pose::Pose & pose,
	protocols::canonical_sampling::MetropolisHastingsMover& metropolis_hastings_mover,
	core::Size& cycle,
	core::Size& temp_level,
	core::Real& temperature
) {
	core::scoring::constraints::ConstraintCOPs csts = pose.constraint_set()->get_all_constraints(); // copy cst info before restart

	utility::file::FileName physical_filename( metropolis_hastings_mover.get_last_checkpoint()+".out" );
	tr.Debug << "restarting from checkpoing file " << physical_filename << std::endl;

	//check existence of file
	temp_level = 1;
	temperature = -1.0;
	cycle = 0;
	//co
	//check for correct tags in file
	tr.Info << "restarting from trajector file " << physical_filename << ". Reading tags now ..." << std::endl;
	io::silent::SilentFileData sfd( physical_filename );
	std::ostringstream replica_id_str;
	replica_id_str << std::setw(3) << std::setfill('0') << jd2::current_replica();
	std::string tag = jd2::current_output_name()+"_"+replica_id_str.str();
	tr.Info << "tag to match: " << tag << std::endl;
	utility::vector1< std::string > matched_tags_in_file;
	bool found = sfd.matched_tags( tag, "last", matched_tags_in_file );
	tr.Info << "if matched_tags?: " << found << std::endl;
	if ( found ) {
		sfd.read_file( sfd.filename(), matched_tags_in_file );
		runtime_assert( sfd.size() == 1 );
		sfd.begin()->fill_pose( pose );

		pose.data().clear(); // otherwise the order of score options will change after restart
		pose.add_constraints( csts ); // recover cst from the copy

		std::string decoy_tag=matched_tags_in_file.front();
		tr.Info << "decoy_tag matched: " << decoy_tag << std::endl;
		Size ind=decoy_tag.find_last_of( '_' );
		//  cycle = utility::string2int(decoy_tag.substr( ind+1 ) )*stride();
		cycle = utility::string2int( decoy_tag.substr( ind+1 )); // the actual trial_number is stored in the decoy tag of the checkpoint file
		tr.Info << "cycle number got from the decoy_tag: " << cycle << std::endl;
		if ( sfd.begin()->has_energy( "temp_level" ) ) {
			temp_level = (Size) sfd.begin()->get_energy( "temp_level" );
		}
		if ( sfd.begin()->has_energy( "temperature" ) ) {
			temperature = sfd.begin()->get_energy( "temperature" );
		}
		tr.Info << "temp_level: " << temp_level << std::endl;
		sfd.begin()->print_scores( tr.Info );
	}
	return found;
}

void
SilentTrajectoryRecorder::initialize_simulation(
	core::pose::Pose & pose,
	protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover,
	core::Size cycle //default=0; non-zero if trajectory is restarted
) {
	Parent::initialize_simulation(pose, metropolis_hastings_mover,cycle);
	tr.Info << std::setprecision( 3 );
	tr.Debug << std::setprecision( 3 );
}

void
SilentTrajectoryRecorder::observe_after_metropolis(
	protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover
)
{ ///this is only for tracer output !
	protocols::moves::MonteCarlo const& mc( *(metropolis_hastings_mover.monte_carlo()) );
	Pose const& pose( mc.last_accepted_pose() );
	if ( step_count() % std::max(stride(),(core::Size)500) == 0 ) {
		jd2::JobOP job( jd2::get_current_job() ) ;
		tr.Info << step_count() << " E=" << pose.energies().total_energy();
		//output what is in job-object (e.g. temperature )
		for ( auto it( job->output_string_real_pairs_begin()), end(job->output_string_real_pairs_end()); it != end; ++it ) {
			tr.Info << " " << it->first << "=" << it->second;
		}
		mc.show_counters();
		tr.Info << std::endl;
	}

	Parent::observe_after_metropolis(metropolis_hastings_mover);
}


} // namespace canonical_sampling
} // namespace protocols
