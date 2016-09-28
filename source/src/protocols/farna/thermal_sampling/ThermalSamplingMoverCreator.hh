// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/farna/thermal_sampling/ThermalSamplingMoverCreator.hh
/// @brief Use a simulated tempering simulation to refine a pose
/// @author Andy Watkins (amw579@nyu.edu)

#ifndef INCLUDED_protocols_farna_thermal_sampling_ThermalSamplingMoverCreator_hh
#define INCLUDED_protocols_farna_thermal_sampling_ThermalSamplingMoverCreator_hh

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace farna {
namespace thermal_sampling {

class ThermalSamplingMoverCreator : public protocols::moves::MoverCreator {

public:

	virtual protocols::moves::MoverOP
	create_mover() const;

	virtual std::string
	keyname() const;

};

} //protocols
} //farna
} //thermal_sampling

#endif //INCLUDED_protocols/farna/thermal_sampling_ThermalSamplingMover_fwd_hh
