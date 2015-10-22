// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/conformation/ResidueFactory.cc
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/conformation/ResidueFactory.hh>

// The point of this class is to limit the number of files
// that #include Residue_
//
// now in flux
//
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.fwd.hh>


namespace core {
namespace conformation {

ResidueOP
ResidueFactory::create_residue( chemical::ResidueType const & rsd_type )
{
	return ResidueOP( new Residue( rsd_type, true /* this superfluous arg to prevent type conversions */) );
}

ResidueOP
ResidueFactory::create_residue(
	chemical::ResidueType const & rsd_type,
	Residue const & current_rsd,
	Conformation const & conf,
	bool preserve_c_beta
)
{
	return ResidueOP( new Residue( rsd_type, current_rsd, conf, preserve_c_beta ) );
}

} // namespace conformation
} // namespace core
