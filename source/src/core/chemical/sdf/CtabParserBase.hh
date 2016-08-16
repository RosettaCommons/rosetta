// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/chemical/sdf/CtabParserBase.hh
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_chemical_CtabParserBase_hh
#define INCLUDED_core_chemical_CtabParserBase_hh

#include <core/chemical/sdf/CtabParserBase.fwd.hh>

#include <core/chemical/sdf/MolFileIOData.fwd.hh>

#include <string>

namespace core {
namespace chemical {
namespace sdf {

class CtabParserBase {
public:

	CtabParserBase() {}
	virtual ~CtabParserBase() {}

	virtual bool parse(std::istream & tablein, std::string const & headerline, MolFileIOMolecule & molecule) = 0;

};

}
}
}


#endif
