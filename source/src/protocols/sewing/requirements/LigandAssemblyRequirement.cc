// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/sewing/requirements/LigandAssemblyRequirement.cc
/// @brief an interface for making Ligand specific Requirements for Assemblies
/// @author Minnie Langlois (minnie@email.unc.edu)


// Unit header or inline function header
#include <protocols/sewing/requirements/LigandAssemblyRequirement.hh>
#include <protocols/sewing/requirements/AssemblyRequirementCreator.hh>
// NOTE: This file should have NO dependencies other than its header.


namespace protocols {
namespace sewing {
namespace requirements {
/*
std::string
LigandAssemblyRequirementCreator::assembly_requirement_ct_namer( std::string tag_name ){
return "assembly_requirement_" + tag_name + "_complex_type";
}
*/
// Defined to prevent pure virtual destructor error at run time.
//LigandAssemblyRequirement::~LigandAssemblyRequirement(){}



} //protocols
} //sewing
} //requirements


