// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/electron_density/DensityZscores.fwd.hh
/// @brief protocol to score local density-fit
/// @author Gabriella Reggiano (reggiano@uw.edu)

#ifndef INCLUDED_protocols_electron_density_DensityZscores_fwd_hh
#define INCLUDED_protocols_electron_density_DensityZscores_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace electron_density {

class DensityZscores;

using DensityZscoresOP = utility::pointer::shared_ptr< DensityZscores >;
using DensityZscoresCOP = utility::pointer::shared_ptr< DensityZscores const >;

} //electron_density
} //protocols

#endif //INCLUDED_protocols_electron_density_DensityZscores_fwd_hh
