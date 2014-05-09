// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /protocols/canonical_sampling/TrialCounterObserver.hh
/// @brief
/// @author

#ifndef INCLUDED_protocols_canonical_sampling_TrialCounterObserver_hh
#define INCLUDED_protocols_canonical_sampling_TrialCounterObserver_hh

// Unit Headers
//#include <protocols/canonical_sampling/TrialCounterObserver.fwd.hh>

// Project Headers
#include <protocols/canonical_sampling/MetropolisHastingsMover.fwd.hh>
#include <protocols/canonical_sampling/ThermodynamicObserver.hh>
#include <protocols/canonical_sampling/MultiTemperatureTrialCounter.hh>
#include <core/id/DOF_ID_Range.hh>
#include <core/id/TorsionID_Range.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

// Utility Headers
#include <core/types.hh>

namespace protocols {
namespace canonical_sampling {

/// @brief Record the acceptance rate for every type of move attempted.
/// @details Most of the work is done by observe_after_metropolis().  Separate 
/// statistics are recorded for moves made at different temperature levels.  
class TrialCounterObserver : public protocols::canonical_sampling::ThermodynamicObserver {

public:

	/// @brief Default constructor.
	TrialCounterObserver();

	/// @brief Default destructor.
	virtual
	~TrialCounterObserver();

	virtual std::string get_name() const;

	virtual protocols::moves::MoverOP	clone() const;

  void
	parse_my_tag(
		 utility::tag::TagCOP tag,
     basic::datacache::DataMap & data,
     protocols::filters::Filters_map const & filters,
     protocols::moves::Movers_map const & movers,
     core::pose::Pose const & pose
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

	/// @brief Return false, as a valid pose is not required for counting trials.
	virtual
	bool
	requires_pose() { return false; }

private:
	protocols::canonical_sampling::MultiTemperatureTrialCounter counters_;
	std::string file_;
	core::Size io_stride_;
}; //end protocols::canonical_sampling::TrialCounterObserver

} //namespace canonical_sampling
} //namespace protocols

#endif // INCLUDED_protocols_canonical_sampling_TrialCounterObserver_HH
