// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/MetricRecorder.cc
///
/// @brief
/// @author


// Unit header or inline function header
#include <protocols/moves/MetricRecorder.hh>
#include <protocols/moves/MetricRecorderCreator.hh>

// Other project headers or inline function headers
#include <core/io/raw_data/ScoreStruct.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <protocols/moves/MetropolisHastingsMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/ThermodynamicMover.hh>  // required for Windows build
#include <protocols/rosetta_scripts/util.hh>
#include <utility/tag/Tag.hh>

// External library headers

// C++ headers
#include <algorithm>
#include <iomanip>
#include <sstream>

// Operating system headers

// Forward declarations


namespace protocols {
namespace moves {

std::string
MetricRecorderCreator::keyname() const {
	return MetricRecorderCreator::mover_name();
}

protocols::moves::MoverOP
MetricRecorderCreator::create_mover() const {
	return new MetricRecorder;
}

std::string
MetricRecorderCreator::mover_name() {
	return "MetricRecorder";
}

MetricRecorder::MetricRecorder() :
	stride_(1),
	step_count_(0)
{}

MetricRecorder::~MetricRecorder()
{}

MetricRecorder::MetricRecorder(
	MetricRecorder const & other
) :
	ThermodynamicObserver(other),
	stride_(other.stride_),
	step_count_(other.step_count_),
	file_name_(other.file_name_),
	torsion_ids_(other.torsion_ids_)
{}

MetricRecorder&
MetricRecorder::operator=( MetricRecorder const & /* other */ )
{
	// assignment not allowed
	runtime_assert(false);
	return * this;
}

MoverOP
MetricRecorder::clone() const
{
	return new MetricRecorder( *this );
}

MoverOP
MetricRecorder::fresh_instance() const
{
	return new MetricRecorder;
}

std::string
MetricRecorder::get_name() const
{
	return "MetricRecorder";
}

void
MetricRecorder::parse_my_tag(
	utility::tag::TagPtr const tag,
	protocols::moves::DataMap & /* data */,
	protocols::filters::Filters_map const & /* filters */,
	protocols::moves::Movers_map const & /* movers */,
	core::pose::Pose const & pose
)
{
	stride_ = tag->getOption< core::Size >( "stride", 100 );
	file_name_ = tag->getOption< std::string >( "filename", "metrics.txt" );

	utility::vector0< utility::tag::TagPtr > const subtags( tag->getTags() );

	for( utility::vector0< utility::tag::TagPtr >::const_iterator subtag_it = subtags.begin(); subtag_it != subtags.end(); ++subtag_it ) {

		TagPtr const subtag = *subtag_it;

		protocols::moves::MoverOP mover;

		if (subtag->getName() == "Torsion") {

			std::string rsd = subtag->getOption< std::string >( "rsd" );
			std::string type = subtag->getOption< std::string >( "type" );
			core::Size torsion = subtag->getOption< core::Size >( "torsion" );
			std::string name = subtag->getOption< std::string >( "name", "" );
			if (subtag->hasOption("name")) name = subtag->getOptions().find("name")->second;
				
			add_torsion(pose, rsd, type, torsion, name);

		} else {
			utility_exit_with_message("Parsed unknown metric type in MetricRecorder");
		}
	}
}

std::string const &
MetricRecorder::file_name() const
{
	return file_name_;
}

void
MetricRecorder::file_name(
	std::string const & file_name
)
{
	file_name_ = file_name;
}

core::Size
MetricRecorder::stride() const
{
	return stride_;
}

void
MetricRecorder::stride(
	core::Size stride
)
{
	stride_ = stride;
}

void
MetricRecorder::add_torsion(
	core::id::TorsionID const & torsion_id,
	std::string name // = ""
)
{
	runtime_assert(torsion_id.valid());

	if (name == "") {
	
		std::ostringstream name_stream;
		name_stream << torsion_id;
		name = name_stream.str();
	}

	torsion_ids_.push_back(std::pair<core::id::TorsionID, std::string>(torsion_id, name));
}

void
MetricRecorder::add_torsion(
	core::pose::Pose const & pose,
	std::string const & rsd,
	std::string type,
	core::Size torsion,
	std::string name // = ""
)
{
	core::Size parsed_rsd = protocols::rosetta_scripts::parse_resnum(rsd, pose);

	std::transform(type.begin(), type.end(), type.begin(), toupper);
	core::id::TorsionType parsed_type;
	if (type == "BB") {
		parsed_type = core::id::BB;
	} else if (type == "CHI") {
		parsed_type = core::id::CHI;
	} else if (type == "JUMP") {
		parsed_type = core::id::JUMP;
	} else {
		utility_exit_with_message("Unknown torsion type");
	}
	
	core::id::TorsionID torsion_id(parsed_rsd, parsed_type, torsion);
	
	add_torsion(torsion_id, name);
}

void
MetricRecorder::reset(
	core::pose::Pose const & pose
)
{
	step_count_ = 0;
	recorder_stream_.close();
	update_after_boltzmann(pose);
}

void
MetricRecorder::update_after_boltzmann(
	core::pose::Pose const & pose
)
{
	if (recorder_stream_.filename() == "") recorder_stream_.open(file_name_);

	if (step_count_ == 0) {

		recorder_stream_ << "Trial" << '\t' << "Score";

		for (core::Size i = 1; i <= torsion_ids_.size(); ++i) {
			recorder_stream_ << '\t' << torsion_ids_[i].second;
		}

		recorder_stream_ << std::endl;
	}

	if (step_count_ % stride_ == 0) {

		recorder_stream_ << step_count_ << '\t' << pose.energies().total_energy();
		
		for (core::Size i = 1; i <= torsion_ids_.size(); ++i) {
			recorder_stream_ << '\t' << pose.torsion(torsion_ids_[i].first);
		}
		
		recorder_stream_ << std::endl;
	}

	++step_count_;
}

void
MetricRecorder::apply(
	core::pose::Pose & pose
)
{
	update_after_boltzmann(pose);
}

void
MetricRecorder::initialize_simulation(
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

	reset(metropolis_hastings_mover.monte_carlo()->last_accepted_pose());

	file_name(original_file_name);
}

void
MetricRecorder::observe_after_metropolis(
	protocols::moves::MetropolisHastingsMover const & metropolis_hastings_mover
)
{
	update_after_boltzmann(metropolis_hastings_mover.monte_carlo()->last_accepted_pose());
}

void
MetricRecorder::finalize_simulation(
	core::pose::Pose & pose,
	protocols::moves::MetropolisHastingsMover const & metropolis_hastings_mover
)
{
	recorder_stream_.close();
}


} // namespace moves
} // namespace protocols
