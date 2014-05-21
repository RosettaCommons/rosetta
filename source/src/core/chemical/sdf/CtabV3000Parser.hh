// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   core/chemical/sdf/CtabV3000Parser.hh
///
/// @brief parse the CTAB table into a residue object
/// @author Sam DeLuca
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_CORE_chemical_sdf_CtabV3000Parser_HH
#define INCLUDED_CORE_chemical_sdf_CtabV3000Parser_HH

#include <core/chemical/sdf/CtabV3000Parser.fwd.hh>
#include <core/chemical/sdf/CtabParserBase.hh>

#include <core/chemical/sdf/MolFileIOData.hh>

#include <string>
#include <istream>

namespace core {
namespace chemical {
namespace sdf {

class CtabV3000Parser: public CtabParserBase {

public:

	virtual ~CtabV3000Parser() {}

	virtual bool parse(std::istream & tablein, std::string const & headerline, MolFileIOMolecule & molecule);

private:

	virtual bool parse_atom_line( std::string line, MolFileIOAtom & atom);
	virtual bool parse_bond_line( std::string line, MolFileIOBond & bond);

	/// @brief Having a '-' at the end of the line is a continuation character for V3000 sdf files.
	void getline(std::istream & istream, std::string line) const;
	/// @brief Split a key value pair on the first equals sign.
	void splitkv(std::string const & kvpair, std::string & key, std::string & value) const;

private:

};

} // sdf
} // chemical
} // core


#endif
