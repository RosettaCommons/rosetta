// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/scoring/membrane/MembraneData.fwd.hh
///
/// @brief  Membrane Base Potential
/// @details Holds access to centroid rotamer pair potential and membrane
///    environemnt pair potential statistics. Loads data in from etables
///    and the database in construction. All immutable data - no setters.
///    Last Modified: 3/28/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_MembraneData_fwd_hh
#define INCLUDED_core_scoring_membrane_MembraneData_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace membrane {

/// @brief Mmebrane Base Potnetial Class
class MembraneData;
typedef utility::pointer::shared_ptr< MembraneData > MembraneDataOP;
typedef utility::pointer::shared_ptr< MembraneData const > MembraneDataCOP;

} // membrane
} // scoring
} // core


#endif // INCLUDED_core_scoring_membrane_MembraneData_fwd_hh

