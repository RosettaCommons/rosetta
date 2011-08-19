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
#include <protocols/moves/TrajectoryRecorderCreator.hh>

// Other project headers or inline function headers
#include <core/io/raw_data/ScoreStruct.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <protocols/moves/MetropolisHastingsMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/ThermodynamicMover.hh>  // required for Windows build
#include <protocols/ScoreMap.hh>
#include <utility/tag/Tag.hh>

// External library headers

// C++ headers
#include <iomanip>

// Operating system headers

// Forward declarations


namespace protocols {
namespace moves {

std::string
TrajectoryRecorderCreator::keyname() const {
	return TrajectoryRecorderCreator::mover_name();
}

protocols::moves::MoverOP
TrajectoryRecorderCreator::create_mover() const {
	return new TrajectoryRecorder;
}

std::string
TrajectoryRecorderCreator::mover_name() {
	return "TrajectoryRecorder";
}

TrajectoryRecorder::TrajectoryRecorder() :
	stride_(1),
	model_count_(0),
	step_count_(0)
{}

TrajectoryRecorder::~TrajectoryRecorder()
{}

TrajectoryRecorder::TrajectoryRecorder(
	TrajectoryRecorder const & other
) :
	ThermodynamicObserver(other),
	stride_(other.stride_),
	model_count_(other.model_count_),
	step_count_(other.step_count_),
	file_name_(other.file_name_)
{}

TrajectoryRecorder&
TrajectoryRecorder::operator=( TrajectoryRecorder const & /* other */ )
{
	// assignment not allowed
	runtime_assert(false);
	return * this;
}

MoverOP
TrajectoryRecorder::clone() const
{
	return new TrajectoryRecorder( *this );
}

MoverOP
TrajectoryRecorder::fresh_instance() const
{
	return new TrajectoryRecorder;
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
	file_name_ = tag->getOption< std::string >( "filename", "traj.pdb" );
}

void
TrajectoryRecorder::write_model(
	core::pose::Pose const & pose
)
{
	if (trajectory_stream_.filename() == "") trajectory_stream_.open(file_name_);

	std::map < std::string, core::Real > score_map;
	std::map < std::string, std::string > string_map;
	protocols::ScoreMap::score_map_from_scored_pose(score_map, pose);
	core::io::raw_data::ScoreStruct score_struct;

	trajectory_stream_ << "MODEL     " << std::setw(4) << model_count_ << std::endl;
	trajectory_stream_ << "REMARK  99 " << step_count_ << std::endl;
	trajectory_stream_ << "REMARK  98 " << pose.energies().total_energy() << std::endl;
	trajectory_stream_ << "REMARK  97 ";
	score_struct.print_header(trajectory_stream_, score_map, string_map, false);
	trajectory_stream_ << "REMARK  97 ";
	score_struct.print_scores(trajectory_stream_, score_map, string_map);

	pose.dump_pdb(trajectory_stream_);
	trajectory_stream_ << "ENDMDL" << std::endl;
}

void
TrajectoryRecorder::reset(
	protocols::moves::MonteCarlo const & mc
)
{
	model_count_ = 0;
	step_count_ = 0;
	write_model(mc.last_accepted_pose());
}

void
TrajectoryRecorder::update_after_boltzmann(
	core::pose::Pose const & pose
)
{
	++step_count_;

	if (step_count_ % stride_ == 0) {
		++model_count_;
		write_model(pose);
	}
}

void
TrajectoryRecorder::update_after_boltzmann(
	protocols::moves::MonteCarlo const & mc
)
{
	update_after_boltzmann(mc.last_accepted_pose());
}

void
TrajectoryRecorder::apply(
	core::pose::Pose & pose
)
{
	update_after_boltzmann(pose);
}

void
TrajectoryRecorder::initialize_simulation(
	core::pose::Pose & pose,
	protocols::moves::MetropolisHastingsMover const & metropolis_hastings_mover
)
{
	std::string original_file_name(file_name());

	if (metropolis_hastings_mover.output_name() != "") {
		std::ostringstream filename;
		filename << metropolis_hastings_mover.output_name() << "_" << file_name();
		file_name(filename.str());
	}

	reset(*(metropolis_hastings_mover.monte_carlo()));

	file_name(original_file_name);
}

void
TrajectoryRecorder::observe_after_metropolis(
	protocols::moves::MetropolisHastingsMover const & metropolis_hastings_mover
)
{
	update_after_boltzmann(*(metropolis_hastings_mover.monte_carlo()));
}

void
TrajectoryRecorder::finalize_simulation(
	core::pose::Pose & pose,
	protocols::moves::MetropolisHastingsMover const & metropolis_hastings_mover
)
{
	trajectory_stream_.close();
}


} // namespace moves
} // namespace protocols
