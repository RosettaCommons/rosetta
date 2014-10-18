// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       LipoFileOptions.fwd.hh
///
/// @brief      Options class - load object that stores Lipid Exposure Data from Membrane Spanning Topology
/// @details    Lipid exposure data calculated using run_lips.pl from octopus spanning
///             topology data. Requires blast.
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_io_LipoFileOptions_fwd_hh
#define INCLUDED_protocols_membrane_io_LipoFileOptions_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace membrane {
namespace io {

/// @brief Options class for reading lipid acc data
/// @details Set options for lipid acc data obj
class LipoFileOptions;
typedef utility::pointer::shared_ptr< LipoFileOptions > LipoFileOptionsOP;
typedef utility::pointer::shared_ptr< LipoFileOptions const > LipoFileOptionsCOP;

} // io
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_io_LipoFileOptions_fwd_hh
