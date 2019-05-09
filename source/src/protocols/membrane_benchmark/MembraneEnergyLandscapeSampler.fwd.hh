// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/membrane_benchmark/MembraneEnergyLandscapeSampler.fwd.hh
/// @brief Sample all points on a 2D membrane energy landscape given implicit model and protein dimensions
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_benchmark_MembraneEnergyLandscapeSampler_fwd_hh
#define INCLUDED_protocols_membrane_benchmark_MembraneEnergyLandscapeSampler_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace membrane_benchmark {

class MembraneEnergyLandscapeSampler;

typedef utility::pointer::shared_ptr< MembraneEnergyLandscapeSampler > MembraneEnergyLandscapeSamplerOP;
typedef utility::pointer::shared_ptr< MembraneEnergyLandscapeSampler const > MembraneEnergyLandscapeSamplerCOP;

} //protocols
} //membrane_benchmark

#endif //INCLUDED_protocols_membrane_benchmark_MembraneEnergyLandscapeSampler_fwd_hh
