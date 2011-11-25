// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/PDBTrajectoryRecorder.cc
///
/// @brief
/// @author


// Unit header or inline function header
#include <protocols/moves/PDBTrajectoryRecorder.hh>
#include <protocols/moves/PDBTrajectoryRecorderCreator.hh>

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
PDBTrajectoryRecorderCreator::keyname() const {
	return PDBTrajectoryRecorderCreator::mover_name();
}

protocols::moves::MoverOP
PDBTrajectoryRecorderCreator::create_mover() const {
	return new PDBTrajectoryRecorder;
}

std::string
PDBTrajectoryRecorderCreator::mover_name() {
	return "PDBTrajectoryRecorder";
}

PDBTrajectoryRecorder::PDBTrajectoryRecorder() {}

PDBTrajectoryRecorder::~PDBTrajectoryRecorder() {}

PDBTrajectoryRecorder::PDBTrajectoryRecorder(
	PDBTrajectoryRecorder const & other
) : TrajectoryRecorder( other ) {}

PDBTrajectoryRecorder&
PDBTrajectoryRecorder::operator=( PDBTrajectoryRecorder const & /* other */ )
{
	// assignment not allowed
	runtime_assert(false);
	return *this;
}

MoverOP
PDBTrajectoryRecorder::clone() const
{
	return new PDBTrajectoryRecorder( *this );
}

MoverOP
PDBTrajectoryRecorder::fresh_instance() const
{
	return new PDBTrajectoryRecorder;
}

std::string
PDBTrajectoryRecorder::get_name() const
{
	return "PDBTrajectoryRecorder";
}

void
PDBTrajectoryRecorder::parse_my_tag(
	utility::tag::TagPtr const tag,
	protocols::moves::DataMap& data /* data */,
	protocols::filters::Filters_map const& filter /* filters */,
	protocols::moves::Movers_map const& mover /* movers */,
	core::pose::Pose const& pose/* pose */
)
{
	Parent::parse_my_tag( tag, data, filter, mover, pose );
}

void
PDBTrajectoryRecorder::write_model(
	core::pose::Pose const & pose
)
{
	if (trajectory_stream_.filename() == "") {
		std::string filename( cumulate() ? file_name() : current_output_name()+"_"+file_name() );
		if ( filename.find( ".pdb" ) == std::string::npos ) {
			filename = filename+".pdb";
		}
		trajectory_stream_.open( filename );
	}

	std::map < std::string, core::Real > score_map;
	std::map < std::string, std::string > string_map;
	protocols::ScoreMap::score_map_from_scored_pose(score_map, pose);
	core::io::raw_data::ScoreStruct score_struct;

	trajectory_stream_ << "MODEL     " << std::setw(4) << model_count() << std::endl;
	trajectory_stream_ << "REMARK  99 " << step_count() << std::endl;
	trajectory_stream_ << "REMARK  98 " << pose.energies().total_energy() << std::endl;
	trajectory_stream_ << "REMARK  97 ";
	score_struct.print_header(trajectory_stream_, score_map, string_map, false);
	trajectory_stream_ << "REMARK  97 ";
	score_struct.print_scores(trajectory_stream_, score_map, string_map);

	pose.dump_pdb(trajectory_stream_);
	trajectory_stream_ << "ENDMDL" << std::endl;
}

void
PDBTrajectoryRecorder::finalize_simulation(
	core::pose::Pose & pose,
	protocols::moves::MetropolisHastingsMover const & metropolis_hastings_mover
)
{
	Parent::finalize_simulation( pose, metropolis_hastings_mover );
	trajectory_stream_.close();
}


} // namespace moves
} // namespace protocols
