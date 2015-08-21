// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_core_conformation_membrane_MembraneParams_hh
#define INCLUDED_core_conformation_membrane_MembraneParams_hh

/// @file  core/conformation/membrane/MembraneParams.hh
///
/// @brief  Membrane Params Enum - map Atom indices to atom types
/// @details Map membrane atom indices to human readable parameters as
///    they are listed in MEM.params
///    Last Modified: 6/28/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

namespace core {
namespace conformation {
namespace membrane {

enum MEM {
	thickness = 1,
	center = 2,
	normal = 3
};

} // membrane
} // conformation
} // core

#endif // INCLDUED_core_conformation_membrane_MembraneParams_hh
