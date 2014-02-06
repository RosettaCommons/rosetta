// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @author Kale Kundert (kale.kundert@ucsf.edu)

#ifndef INCLUDED_protocols_canonical_sampling_FixedTemperatureController_hh
#define INCLUDED_protocols_canonical_sampling_FixedTemperatureController_hh

// Unit Headers
#include <protocols/canonical_sampling/FixedTemperatureController.fwd.hh>
#include <protocols/canonical_sampling/TemperatureController.hh>

// Core Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Protocols headers
#include <protocols/moves/Mover.fwd.hh>

namespace protocols {
namespace canonical_sampling {

/// @brief Maintain a constant temperature.
/// @details This is the default temperature controller used by 
/// MetropolisHastingsMover.
class FixedTemperatureController : public TemperatureController {
public:

	/// @brief Constructor with temperature parameter.
	FixedTemperatureController( core::Real temp );

	/// @brief Default destructor.
	~FixedTemperatureController();

	/// @brief Return a copy of this mover.
	protocols::moves::MoverOP clone() const;

	/// @brief Ignore any command-line options.
	void init_from_options() {}

	virtual
	std::string
	get_name() const { return "FixedTemperatureContoller"; }

	/// @brief Return the same constant temperature every time.
	virtual core::Real temperature_move(
			core::pose::Pose & pose,
			MetropolisHastingsMover & mover,
			core::Real score);
};

} //namespace canonical_sampling
} //namespace protocols

#endif
