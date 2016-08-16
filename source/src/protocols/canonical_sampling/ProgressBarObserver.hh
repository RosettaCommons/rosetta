// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @author Kale Kundert (kale.kundert@ucsf.edu)

#ifndef INCLUDED_protocols_canonical_sampling_ProgressBarObserver_hh
#define INCLUDED_protocols_canonical_sampling_ProgressBarObserver_hh

// Unit Headers
#include <protocols/canonical_sampling/ProgressBarObserver.fwd.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMover.fwd.hh>
#include <protocols/canonical_sampling/ThermodynamicObserver.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Protocol headers
#include <protocols/moves/Mover.hh>

namespace protocols {
namespace canonical_sampling {

/// @brief Print an progress bar that updated on each iteration.
/// @details The progress bar is currently written to stdout, which of course
/// is not a tracer.  I'm not really sure if this is a good idea or not.  I
/// think it'd be cool to detect whether or not rosetta is attached to a TTY
/// and change how the progress bar is drawn depending.  For example, when
/// writing to files it's nicer to not write carriage return '\r' characters.
class ProgressBarObserver : public ThermodynamicObserver {

public:

	/// @brief Default constructor.
	ProgressBarObserver();

	/// @brief Copy constructor.
	ProgressBarObserver( ProgressBarObserver const & );

	/// @brief Default destructor.
	~ProgressBarObserver();

	std::string get_name() const {
		return "ProgressBarObserver";
	}

	protocols::moves::MoverOP clone() const;

	void observe_after_metropolis(
		MetropolisHastingsMover const & metropolis_hastings_mover);

	/// @brief Return false, as a valid pose is not required for printing a
	/// progress bar.
	bool requires_pose() { return false; }

private:
	core::Size progress_;
};

} // namespace canonical_sampling
} // namespace protocols

#endif
