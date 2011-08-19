// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @author Ian W. Davis


#ifndef INCLUDED_core_chemical_CSDAtomType_hh
#define INCLUDED_core_chemical_CSDAtomType_hh


// C++ headers
#include <string>

namespace core {
namespace chemical {

/// @brief Basic "conformational" atom type derived from analysis of Cambridge Structure Database.
class CSDAtomType {

public:

	///  @brief Construct a new CSDAtomType with its name
	CSDAtomType( std::string const & name_in):
		name_( name_in )
	{}

	/// @brief Return the name of the CSDAtomType
	std::string const& name() const { return name_; };

	// data
private:

	/// name
	std::string const name_;
};


} // chemical
} // core



#endif // INCLUDED_core_chemical_CSDAtomType_HH
