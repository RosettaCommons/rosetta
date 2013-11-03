// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/canonical_sampling/TrajectoryRecorder.hh
///
/// @brief
/// @author


#ifndef INCLUDED_protocols_canonical_sampling_TrajectoryRecorder_hh
#define INCLUDED_protocols_canonical_sampling_TrajectoryRecorder_hh


// Project forward headers
#include <protocols/canonical_sampling/TrajectoryRecorder.fwd.hh>

// Project headers
#include <protocols/canonical_sampling/ThermodynamicObserver.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>

#include <utility/io/ozstream.hh>
#include <utility/vector1.hh>

// External library headers

// C++ headers
#include <string>

// Operating system headers

// Forward declarations

namespace protocols {
namespace canonical_sampling {

/// @brief
class TrajectoryRecorder : public protocols::canonical_sampling::ThermodynamicObserver {
public:
	/// @brief Associates relevant options with the TemperedDocking class
	static void register_options();

	/// @brief Constructor
	TrajectoryRecorder();

	/// @brief Destructor
	~TrajectoryRecorder();

	/// @brief Copy constructor
	TrajectoryRecorder( TrajectoryRecorder const & );

private:
	//assignment not allowed
	TrajectoryRecorder&
	operator=( TrajectoryRecorder const & );

public:
	virtual	std::string	get_name() const;

	virtual	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);

	std::string const &	file_name() const {
		return file_name_;
	}

	void file_name( std::string const & file_name )	{
		file_name_ = file_name;
	}

	core::Size model_count() const {
		return model_count_;
	}

	core::Size step_count() const {
		return step_count_;
	}

	core::Size stride() const {
		return stride_;
	}

	bool cumulate_jobs() const {
		return cumulate_jobs_;
	}

	bool cumulate_replicas() const {
		return cumulate_replicas_;
	}

	void stride( core::Size stride ) {
		stride_ = stride;
	}

	virtual void reset(
		protocols::moves::MonteCarlo const& mc,
		protocols::canonical_sampling::MetropolisHastingsMoverCAP metropolis_hastings_mover = 0
	);

	void update_after_boltzmann(
		core::pose::Pose const & pose,
		protocols::canonical_sampling::MetropolisHastingsMoverCAP metropolis_hastings_mover = 0
	);

	void update_after_boltzmann(
		protocols::moves::MonteCarlo const& mc
	);

	virtual void apply( core::pose::Pose& pose );

	virtual void initialize_simulation(
		core::pose::Pose & pose,
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover,
		core::Size cycle   //non-zero if trajectory is restarted
	);

	virtual	void observe_after_metropolis(
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover
	);

protected:

	virtual void 	write_model(
		core::pose::Pose const & pose,
		protocols::canonical_sampling::MetropolisHastingsMoverCAP metropolis_hastings_mover = 0
	) = 0;

private:
	core::Size stride_;
	core::Size model_count_;
	core::Size step_count_;
	std::string file_name_;
	bool cumulate_jobs_;
	bool cumulate_replicas_;

	static bool options_registered_;
}; // TrajectoryRecorder


} // namespace canonical_sampling
} // namespace protocols


#endif // INCLUDED_protocols_canonical_sampling_TrajectoryRecorder_HH
