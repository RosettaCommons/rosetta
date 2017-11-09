// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/legacy_sewing/conformation/AssemblyFactory.cc
///
/// @brief A super-simple factory for generating Assembly classes
///
/// @author Tim Jacobs

//Unit headers
#include <protocols/legacy_sewing/conformation/AssemblyFactory.hh>

//Package headers
#include <protocols/legacy_sewing/conformation/Assembly.hh>
#include <protocols/legacy_sewing/conformation/ContinuousAssembly.hh>
#include <protocols/legacy_sewing/conformation/DisembodiedAssembly.hh>

#include <protocols/legacy_sewing/conformation/Model.hh>

namespace protocols {
namespace legacy_sewing  {

AssemblyOP
AssemblyFactory::create_assembly(
	std::string assembly_name
) {

	if ( assembly_name == "continuous" ) {
		return AssemblyOP( new ContinuousAssembly() );
	} else if ( assembly_name == "discontinuous" ) {
		return AssemblyOP( new DisembodiedAssembly() );
	} else {
		utility_exit_with_message("Invalid Assembly name given to AssemblyFactory!");
	}
	return 0;
}

} //legacy_sewing namespace
} //protocols namespace
