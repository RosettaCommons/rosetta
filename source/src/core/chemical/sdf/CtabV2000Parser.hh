// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/sdf/CtabV2000Parser.hh
///
/// @brief parse the CTAB table into a residue object
/// @author Sam DeLuca
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_CORE_chemical_sdf_CtabV2000Parser_HH
#define INCLUDED_CORE_chemical_sdf_CtabV2000Parser_HH

#include <core/chemical/sdf/CtabV2000Parser.fwd.hh>
#include <core/chemical/sdf/CtabParserBase.hh>

#include <core/chemical/sdf/MolFileIOData.hh>

#include <istream>

namespace core {
namespace chemical {
namespace sdf {

class CtabV2000Parser: public CtabParserBase {

public:

	virtual ~CtabV2000Parser() {}

	virtual bool parse(std::istream & tablein, std::string const & headerline, MolFileIOMolecule & molecule);

private:

	virtual bool parse_atom_line( std::string line, MolFileIOAtom & atom);
	virtual bool parse_bond_line( std::string line, MolFileIOBond & bond);
	virtual bool parse_property_line( std::string line, MolFileIOMolecule & mol);

private:

};

} // sdf
} // chemical
} // core


#endif
