// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
