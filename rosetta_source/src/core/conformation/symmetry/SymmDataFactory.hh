// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

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
