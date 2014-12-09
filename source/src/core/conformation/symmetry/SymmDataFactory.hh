// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief  conformation container
/// @file   core/conformation/symmetry/SymmDescription.hh
/// @author Ingemar Andre


#ifndef INCLUDED_core_conformation_symmetry_SymmDataFactory_hh
#define INCLUDED_core_conformation_symmetry_SymmDataFactory_hh

#include <core/conformation/symmetry/SymmData.fwd.hh>

// Unit headers

namespace core {
namespace conformation {
namespace symmetry {

class SymmDataFactory
{
	public:

	SymmDataOP
	create_symm_data();
};

} // core
} // conformation
} // core
#endif
