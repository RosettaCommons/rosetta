// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/thermal_sampling/RECCES_Mover.fwd.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_thermal_sampling_RECCES_Mover_FWD_HH
#define INCLUDED_protocols_thermal_sampling_RECCES_Mover_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace thermal_sampling {

class RECCES_Mover;
typedef utility::pointer::shared_ptr< RECCES_Mover > RECCES_MoverOP;
typedef utility::pointer::shared_ptr< RECCES_Mover const > RECCES_MoverCOP;

} //thermal_sampling
} //protocols

#endif
