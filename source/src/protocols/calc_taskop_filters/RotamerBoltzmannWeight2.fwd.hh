// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/calc_taskop_filters/RotamerBoltzmannWeight2.fwd.hh
/// @brief Next-generation RotamerBoltzmannWeight filter
/// @author Tom Linsky (tlinsky@uw.edu)


#ifndef INCLUDED_protocols_simple_filters_RotamerBoltzmannWeight2_fwd_hh
#define INCLUDED_protocols_simple_filters_RotamerBoltzmannWeight2_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>



// Forward
namespace protocols {
namespace calc_taskop_filters {

class RotamerBoltzmannWeight2;

typedef utility::pointer::shared_ptr< RotamerBoltzmannWeight2 > RotamerBoltzmannWeight2OP;
typedef utility::pointer::shared_ptr< RotamerBoltzmannWeight2 const > RotamerBoltzmannWeight2COP;



} //protocols
} //simple_filters


#endif //INCLUDED_protocols_simple_filters_RotamerBoltzmannWeight2_fwd_hh





