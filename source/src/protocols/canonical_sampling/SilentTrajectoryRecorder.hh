// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/canonical_sampling/SilentTrajectoryRecorder.hh
///
/// @brief
/// @author


#ifndef INCLUDED_protocols_canonical_sampling_SilentTrajectoryRecorder_hh
#define INCLUDED_protocols_canonical_sampling_SilentTrajectoryRecorder_hh


// Project forward headers
#include <protocols/canonical_sampling/SilentTrajectoryRecorder.fwd.hh>
#include <protocols/canonical_sampling/TrajectoryRecorder.hh>


// Project headers
#include <protocols/canonical_sampling/ThermodynamicObserver.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <utility/io/ozstream.hh>


// External library headers


// C++ headers
#include <string>


// Operating system headers


// Forward declarations


namespace protocols {
namespace canonical_sampling {


	/// @brief
class SilentTrajectoryRecorder : public protocols::canonical_sampling::TrajectoryRecorder {
public:
	typedef TrajectoryRecorder Parent;
	/// @brief Constructor
	SilentTrajectoryRecorder();

	/// @brief Copy constructor
	SilentTrajectoryRecorder( SilentTrajectoryRecorder const & );

public:
	/// @brief Associates relevant options with the TemperedDocking class
	static void register_options();

	virtual	protocols::moves::MoverOP	clone() const;

	virtual	protocols::moves::MoverOP	fresh_instance() const;

	virtual	std::string	get_name() const;

	virtual	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);

	virtual void initialize_simulation(
		core::pose::Pose & pose,
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover,
		core::Size cycle   //non-zero if trajectory is restarted
	);

	virtual	void observe_after_metropolis(
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover
	);

	virtual
	bool
	restart_simulation(
			 core::pose::Pose & pose,
			 protocols::canonical_sampling::MetropolisHastingsMover& metropolis_hastings_mover,
			 core::Size& cycle,
			 core::Size& temp_level,
			 core::Real& temperature
	);

protected:

	virtual void 	write_model(
		core::pose::Pose const & pose,
		protocols::canonical_sampling::MetropolisHastingsMoverCAP metropolis_hastings_mover = 0
	);

	core::Size score_stride_;

	static bool options_registered_;
}; // SilentTrajectoryRecorder


} // namespace canonical_sampling
} // namespace protocols


#endif // INCLUDED_protocols_canonical_sampling_SilentTrajectoryRecorder_HH
