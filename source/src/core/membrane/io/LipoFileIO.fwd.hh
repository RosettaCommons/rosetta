// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       LipoFileIO.fwd.hh
///
/// @brief      Lipid Exposure Data from Membrane Spanning Topology
/// @details    Lipid exposure data calculated using run_lips.pl from octopus spanning
///             topology data. Requires blast.
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_LipidFileIO_fwd_hh
#define INCLUDED_core_membrane_io_LipidFileIO_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace membrane {
namespace io {

/// @brief Lipid Accessibility Data
/// @details Read lipid accessibility data into LipidAcc Obj
class LipoFileIO;
typedef utility::pointer::owning_ptr< LipoFileIO > LipoFileIOOP;
typedef utility::pointer::owning_ptr< LipoFileIO const > LipoFileIOCOP;

} // io
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_LipoFileIO_fwd_hh

