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
/// @file   core/conformation/Conformation.hh
/// @author Ingemar Andre

// Unit headers
#include <core/conformation/symmetry/SymmDataFactory.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <basic/tracer.hh>


namespace core {
namespace conformation {
namespace symmetry {

static thread_local basic::Tracer TR( "core.conformation.symmetry.symmdatafactory" );

SymmDataOP
SymmDataFactory::create_symm_data()
{
	return new SymmData();
}

} // core
} // conformation
} // core
