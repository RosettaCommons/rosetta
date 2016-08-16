// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/flexpack/annealer/FlexbbSimAnnealer.fwd.hh
/// @brief  Forward declarations and typedefs for the FlexbbSimAnnealer class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_flexpack_annealer_FlexbbSimAnnealer_fwd_hh
#define INCLUDED_protocols_flexpack_annealer_FlexbbSimAnnealer_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace flexpack {
namespace annealer {

class FlexbbSimAnnealer;

typedef utility::pointer::shared_ptr< FlexbbSimAnnealer > FlexbbSimAnnealerOP;

}
}
}

#endif
