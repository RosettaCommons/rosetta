// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/canonical_sampling/MetricRecorder.cc
///
/// @brief
/// @author


// Unit header or inline function header
#include <protocols/canonical_sampling/MetricRecorder.hh>
#include <protocols/canonical_sampling/MetricRecorderCreator.hh>

// Other project headers or inline function headers
#include <core/pose/Pose.hh>
//Gabe testing
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
//#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/types.hh>
#include <core/id/AtomID.hh>
//End gabe testing
#include <core/scoring/Energies.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/util.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>
#include <protocols/canonical_sampling/TemperingBase.hh>
#include <protocols/moves/MonteCarlo.hh>
// AUTO-REMOVED #include <protocols/canonical_sampling/ThermodynamicMover.hh>  // required for Windows build
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <utility/tag/Tag.hh>

// External library headers

// C++ headers
#include <algorithm>
// AUTO-REMOVED #include <iomanip>
#include <sstream>

#include <core/id/TorsionID.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


// Operating system headers

// Forward declarations


namespace protocols {
namespace canonical_sampling {

std::string
MetricRecorderCreator::keyname() const {
	return MetricRecorderCreator::mover_name();
}

protocols::moves::MoverOP
MetricRecorderCreator::create_mover() const {
	return protocols::moves::MoverOP( new MetricRecorder );
}

std::string
MetricRecorderCreator::mover_name() {
	return "MetricRecorder";
}

MetricRecorder::MetricRecorder() :
	stride_(1),
	cumulate_jobs_(false),
	cumulate_replicas_(false),
	prepend_output_name_(false),
	step_count_(0),
	last_flush_(0)
{}

MetricRecorder::~MetricRecorder()
{}

MetricRecorder::MetricRecorder(
	MetricRecorder const & other
) :
	protocols::canonical_sampling::ThermodynamicObserver(other),
	file_name_(other.file_name_),
	stride_(other.stride_),
	cumulate_jobs_(other.cumulate_jobs_),
	cumulate_replicas_(other.cumulate_replicas_),
	prepend_output_name_(other.prepend_output_name_),
	step_count_(other.step_count_),
	torsion_ids_(other.torsion_ids_),
	last_flush_(other.last_flush_)
{}

MetricRecorder&
MetricRecorder::operator=( MetricRecorder const & /* other */ )
{
	// assignment not allowed
	runtime_assert(false);
	return * this;
}

protocols::moves::MoverOP
MetricRecorder::clone() const
{
	return protocols::moves::MoverOP( new protocols::canonical_sampling::MetricRecorder( *this ) );
}

protocols::moves::MoverOP
MetricRecorder::fresh_instance() const
{
	return protocols::moves::MoverOP( new MetricRecorder );
}

std::string
MetricRecorder::get_name() const
{
	return "MetricRecorder";
}

void
MetricRecorder::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /* data */,
	protocols::filters::Filters_map const & /* filters */,
	protocols::moves::Movers_map const & /* movers */,
	core::pose::Pose const & pose
)
{
	stride_ = tag->getOption< core::Size >( "stride", 100 );
	cumulate_jobs_ = tag->getOption< bool >( "cumulate_jobs", false );
	cumulate_replicas_ = tag->getOption< bool >( "cumulate_replicas", false );
	prepend_output_name_ = tag->getOption< bool >( "prepend_output_name", false );
	file_name_ = tag->getOption< std::string >( "filename", "metrics.txt" );

	utility::vector0< utility::tag::TagCOP > const subtags( tag->getTags() );

	for( utility::vector0< utility::tag::TagCOP >::const_iterator subtag_it = subtags.begin(); subtag_it != subtags.end(); ++subtag_it ) {

		TagCOP const subtag = *subtag_it;

		protocols::moves::MoverOP mover;

		if (subtag->getName() == "Torsion") {

			std::string rsd = subtag->getOption< std::string >( "rsd" );
			std::string type = subtag->getOption< std::string >( "type" );
			core::Size torsion = subtag->getOption< core::Size >( "torsion" );
			std::string name = subtag->getOption< std::string >( "name", "" );
			if (subtag->hasOption("name")) name = subtag->getOptions().find("name")->second;

			add_torsion(pose, rsd, type, torsion, name);

		} else if (subtag->getName() == "AllChi") {
			for ( Size i = 1; i <= pose.total_residue(); ++i) {
				for (Size j = 1; j <= pose.residue_type(i).nchi(); ++j) {
					std::ostringstream name_stream;
					name_stream << pose.residue_type(i).name3() << "_" << i << "_Chi" << j ;
					std::ostringstream res_id_str;
					res_id_str << i;
					add_torsion(pose, res_id_str.str(), "CHI", j, name_stream.str());
				}
			}

		} else if (subtag->getName() == "AllBB") {
			for ( Size i = 1; i <= pose.total_residue(); ++i) {
				std::ostringstream res_id_str;
				res_id_str << i;
				add_torsion(pose, res_id_str.str(), "BB", 1, pose.residue_type(i).name3() + "_" + res_id_str.str() + "_phi");
				add_torsion(pose, res_id_str.str(), "BB", 2, pose.residue_type(i).name3() + "_" + res_id_str.str() + "_psi");
				add_torsion(pose, res_id_str.str(), "BB", 3, pose.residue_type(i).name3() + "_" + res_id_str.str() + "_omega");
			}

		} else {
			utility_exit_with_message("Parsed unknown metric type in MetricRecorder" + subtag->getName());
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

bool
MetricRecorder::cumulate_jobs() const
{
	return cumulate_jobs_;
}

void
MetricRecorder::cumulate_jobs(
	bool cumulate_jobs
)
{
	cumulate_jobs_ = cumulate_jobs;
}

bool
MetricRecorder::cumulate_replicas() const
{
	return cumulate_replicas_;
}

void
MetricRecorder::cumulate_replicas(
	bool cumulate_replicas
)
{
	cumulate_replicas_ = cumulate_replicas;
}

bool
MetricRecorder::prepend_output_name() const
{
	return prepend_output_name_;
}

void
MetricRecorder::prepend_output_name(
	bool prepend_output_name
)
{
	prepend_output_name_ = prepend_output_name;
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
	core::Size parsed_rsd = core::pose::parse_resnum(rsd, pose);

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
	core::pose::Pose const & pose,
	protocols::canonical_sampling::MetropolisHastingsMover const * metropolis_hastings_mover //= 0
)
{
	step_count_ = 0;
	recorder_stream_.close();
	update_after_boltzmann(pose, metropolis_hastings_mover);
}

void
MetricRecorder::update_after_boltzmann(
	core::pose::Pose const & pose,
	protocols::canonical_sampling::MetropolisHastingsMover const * metropolis_hastings_mover //= 0
)
{
	if (recorder_stream_.filename() == "") {
		std::ostringstream file_name_stream;
		if (prepend_output_name_ && !metropolis_hastings_mover) {
			file_name_stream << protocols::jd2::JobDistributor::get_instance()->current_output_name() << '_';
		}
		file_name_stream << file_name_;
		recorder_stream_.open(file_name_stream.str());
	}

	protocols::jd2::JobCOP job( protocols::jd2::get_current_job() );
	core::Size nstruct_index( job ? job->nstruct_index() : 1 );
	std::string output_name( metropolis_hastings_mover ? metropolis_hastings_mover->output_name() : "" );

	core::Size replica = protocols::jd2::current_replica();

	TemperingBaseCOP tempering = 0;
	if (metropolis_hastings_mover) {
		tempering = utility::pointer::dynamic_pointer_cast< TemperingBase const >( metropolis_hastings_mover->tempering() );
	}

	if (step_count_ == 0) {

		// output header if not cumulating, replica exchange inactive, or this is the first replica
		if ((!cumulate_jobs_ || nstruct_index == 1) && (!cumulate_replicas_ || replica <= 1)) {

			recorder_stream_ << "Trial";
			if (cumulate_jobs_ && output_name.length()) recorder_stream_ << '\t' << "Job";
			if (cumulate_replicas_ && replica) recorder_stream_ << '\t' << "Replica";
			if (tempering) recorder_stream_ << '\t' << "Temperature";
			recorder_stream_ << '\t' << "Score";

			for (core::Size i = 1; i <= torsion_ids_.size(); ++i) {
				recorder_stream_ << '\t' << torsion_ids_[i].second;
			}

			recorder_stream_ << std::endl;
			recorder_stream_.flush();
		}
		last_flush_ = time(NULL);
	}

	if (step_count_ % stride_ == 0) {

		recorder_stream_ << step_count_;
		if (cumulate_jobs_ && output_name.length()) recorder_stream_ << '\t' << output_name;
		if (cumulate_replicas_ && replica) recorder_stream_ << '\t' << replica;
		if (tempering) recorder_stream_ << '\t' << metropolis_hastings_mover->monte_carlo()->temperature();
		recorder_stream_ << '\t' << pose.energies().total_energy();

		for (core::Size i = 1; i <= torsion_ids_.size(); ++i) {
			recorder_stream_ << '\t' << pose.torsion(torsion_ids_[i].first);
		}

		recorder_stream_ << std::endl;

		time_t now = time(NULL);
		if (now-last_flush_ > 10) {
			recorder_stream_.flush();
			last_flush_ = now;
		}
	}

	++step_count_;
}

void
MetricRecorder::apply(
	core::pose::Pose & pose
)
{
	update_after_boltzmann(pose, 0);
}

void
MetricRecorder::initialize_simulation(
	core::pose::Pose &,
	protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover,
	core::Size //default=0; non-zero if trajectory is restarted
)
{
	std::string original_file_name(file_name());

	file_name(metropolis_hastings_mover.output_file_name(file_name(), cumulate_jobs_, cumulate_replicas_));

	reset(
		metropolis_hastings_mover.monte_carlo()->last_accepted_pose(),
		&metropolis_hastings_mover
	);

	file_name(original_file_name);
}

void
MetricRecorder::observe_after_metropolis(
	protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover
)
{
	update_after_boltzmann(
		metropolis_hastings_mover.monte_carlo()->last_accepted_pose(),
		&metropolis_hastings_mover
	);
}

void
MetricRecorder::finalize_simulation(
	core::pose::Pose &,
	protocols::canonical_sampling::MetropolisHastingsMover const &
)
{
	recorder_stream_.close();
}


} // namespace canonical_sampling
} // namespace protocols
