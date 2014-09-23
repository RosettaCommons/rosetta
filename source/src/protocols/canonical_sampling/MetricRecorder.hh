// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/canonical_sampling/MetricRecorder.hh
/// @brief
/// @author

#ifndef INCLUDED_protocols_canonical_sampling_MetricRecorder_hh
#define INCLUDED_protocols_canonical_sampling_MetricRecorder_hh

// Project forward headers
#include <protocols/canonical_sampling/MetricRecorder.fwd.hh>

// Project headers
#include <protocols/canonical_sampling/ThermodynamicObserver.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/io/ozstream.hh>

// C++ headers
#include <ctime>
#include <string>

#include <core/id/TorsionID.fwd.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace canonical_sampling {

/// @brief Periodically output miscellaneous information.
/// @details This class is capable of writing out a variety of data related to 
/// the trajectory.  This includes the job name, the replica, the temperature, 
/// and the score.  Any number of torsion angles can also be added to the 
/// report using add_torsion().  Methods are also provided for specifying an 
/// output filename.  Most of the IO work is done by update_after_boltzmann().  

class MetricRecorder : public ThermodynamicObserver
{
public: // Creation

	/// @brief Default constructor.
	MetricRecorder();

	/// @brief Default destructor.
	~MetricRecorder();

	/// @brief Copy constructor.
	MetricRecorder( MetricRecorder const & );

	/// @brief Assignment operator.
	MetricRecorder&
	operator=( MetricRecorder const & );

public: // Methods

	/// @brief Return a copy of this mover.
	virtual
	protocols::moves::MoverOP
	clone() const;

	/// @brief Return a newly instantiated mover.
	virtual
	protocols::moves::MoverOP
	fresh_instance() const;

	/// @brief Return the name of this mover.
	virtual
	std::string
	get_name() const;

	/// @brief Use a RosettaScripts tag to configure this mover.
	virtual
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);

	/// @brief Return true.  This mover needs to be reinitialized for each job.
	virtual
	bool
	reinitialize_for_each_job() const { return true; };

	/// @brief Return the name of the file being written to.
	std::string const &
	file_name() const;

	/// @brief Set the name of the file being written to.
	void
	file_name(
		std::string const & file_name
	);

	/// @brief Return the frequency with which data is written.
	core::Size
	stride() const;

	/// @brief Set the frequency with which data is written.
	void
	stride(
		core::Size stride
	);

	/// @brief Return true if every job is being reported to the same file.
	bool
	cumulate_jobs() const;

	/// @brief Indicate whether or not every job should be reported to the same 
	/// file.
	void
	cumulate_jobs(
		bool cumulate_jobs
	);

	/// @brief Return true if every replica is being reported to the same file.
	bool
	cumulate_replicas() const;

	/// @brief Indicate whether or not every replica should be reported to the 
	/// same file.
	void
	cumulate_replicas(
		bool cumulate_replicas
	);

	/// @brief Return true if the job name should be prepended onto the output 
	/// filename.
	bool
	prepend_output_name() const;

	/// @brief Indicate whether or not the job name should be prepended onto the 
	/// output filename.
	void
	prepend_output_name(
		bool prepend_output_name
	);

	/// @brief Include the given torsion in the output.
	void
	add_torsion(
		core::id::TorsionID const & torsion_id,
		std::string name = ""
	);

	/// @brief Include the given torsion in the output.
	void
	add_torsion(
		core::pose::Pose const & pose,
		std::string const & rsd,
		std::string type,
		core::Size torsion,
		std::string name = ""
	);

	/// @brief Truncate the output file and rewrite the output header.
	/// @details This method may not actually truncate the output file.  It 
	/// really just closes and reopens the file, and I'm not sure whether or not 
	/// it picks a new name when it does the reopening.
	void
	reset(
		core::pose::Pose const & pose,
		protocols::canonical_sampling::MetropolisHastingsMover const * metropolis_hastings_mover = 0
	);

	/// @brief Write information like temperature, score, and torsion angles to a 
	/// file.
	void
	update_after_boltzmann(
		core::pose::Pose const & pose,
		protocols::canonical_sampling::MetropolisHastingsMover const * metropolis_hastings_mover = 0
	);

	/// @brief Just invoke update_after_boltzmann() with a const pose.
	virtual
	void
	apply(
		core::pose::Pose & pose
	);

	virtual
	void
	initialize_simulation(
		core::pose::Pose & pose,
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover,
		core::Size cycle   //non-zero if trajectory is restarted
	);

	virtual
	void
	observe_after_metropolis(
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover
	);

	virtual
	void
	finalize_simulation(
		core::pose::Pose & pose,
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover
	);

private: // Fields

	std::string file_name_;
	core::Size stride_;
	bool cumulate_jobs_;
	bool cumulate_replicas_;
	bool prepend_output_name_;
	core::Size step_count_;
	utility::io::ozstream recorder_stream_;
	utility::vector1<std::pair<core::id::TorsionID, std::string> > torsion_ids_;
	time_t last_flush_;

}; // MetricRecorder


} // namespace canonical_sampling
} // namespace protocols


#endif // INCLUDED_protocols_canonical_sampling_MetricRecorder_HH
