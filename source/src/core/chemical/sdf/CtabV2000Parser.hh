// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   core/chemical/sdf/ctab_parsher.hh
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
#include <string>

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
