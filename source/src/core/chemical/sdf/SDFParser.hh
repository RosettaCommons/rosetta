// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/chemical/sdf/SDFParser.hh
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_chemical_SDFParser_hh
#define INCLUDED_core_chemical_SDFParser_hh

#include <core/chemical/sdf/SDFParser.fwd.hh>

#include <core/chemical/sdf/MolFileIOData.fwd.hh>

#include <core/types.hh>

#include <utility/vector1.hh>

namespace core {
namespace chemical {
namespace sdf {

class SDFParser {
public:

	SDFParser() {}

	virtual ~SDFParser() {}

	/// @brief parse the given input stream.
	/// n_entries are the maximum number of entries to parse - of zero parse all the remaining ones.
	utility::vector1< MolFileIOMoleculeOP > parse( std::istream & filein, core::Size n_entries = 0 );

private:

	/// @brief Parse the optional data lines in sdf file.
	///  Will consume lines up to and including the EOF or the '$$$$' delimeter.
	void parse_optional_data( std::istream & filein, MolFileIOMolecule & molecule );
	/// @brief Ignore everything up to and including the next '$$$$' entry delimeter.
	void eat_until_delimiter( std::istream & filein ) const;


};

}
}
}


#endif
