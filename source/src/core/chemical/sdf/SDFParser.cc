// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file  core/chemical/sdf/SDFParser.cc
///
/// @brief
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <core/chemical/sdf/SDFParser.hh>

#include <core/chemical/sdf/CtabV2000Parser.hh>
#include <core/chemical/sdf/CtabV3000Parser.hh>
#include <core/chemical/sdf/MolFileIOData.hh>

#include <basic/Tracer.hh>
#include <utility/string_util.hh>

#include <string>
#include <sstream>
#include <istream>

namespace core {
namespace chemical {
namespace sdf {

static basic::Tracer TR("core.chemical.sdf.SDFParser");

utility::vector1< MolFileIOMoleculeOP >
SDFParser::parse(std::istream & filein ) {
	CtabV2000Parser V2000parser;
	CtabV3000Parser V3000parser;

	utility::vector1< MolFileIOMoleculeOP > molecules;
	std::string name, line2, comments, versionline;

	while( filein ) {
		std::getline(filein,name);
		utility::trim( name, " \n" ); //modify in place
		std::getline(filein,line2);
		std::getline(filein,comments);
		std::getline(filein,versionline);

		if( ! filein ) {
			break;
		}
		if( versionline.size() < 39 ) {
			TR.Warning << "Warning: SDF header line too short for: '" << name << "'" << std::endl;
			eat_until_delimiter( filein );
			continue;
		}
		std::string version(versionline,33,6);
		MolFileIOMoleculeOP molecule( new MolFileIOMolecule );
		molecule->name( name );
		if ( version == "V3000" ) {
			if( ! V3000parser.parse(filein, versionline, *molecule) ) {
				TR.Warning << "Skipping V3000 sdf file entry for " << name << "'" << std::endl;
				eat_until_delimiter( filein );
				continue;
			}
		} else if( version == "V2000" ) {
			if( ! V2000parser.parse(filein, versionline, *molecule) ) {
				TR.Warning << "Skipping V2000 sdf file entry for " << name << "'" << std::endl;
				eat_until_delimiter( filein );
				continue;
			}
		} else {
			// Try parsing as V2000 - sometimes these omit the version line.
			if( ! V2000parser.parse(filein, versionline, *molecule) ) {
				TR.Warning << "Attempted to parse '" << name << "' as V2000 but failed. Skipping." << std::endl;
				eat_until_delimiter( filein );
				continue;
			}
		}

		// Will go up to and including the '$$$$' delimeter.
		parse_optional_data( filein, *molecule );

		molecules.push_back( molecule );
	}

	return molecules;
}

void
SDFParser::parse_optional_data( std::istream & filein, MolFileIOMolecule & molecule ) {
	std::string header;
	std::string line;
	std::string entry;

	for( std::getline(filein,header); filein && ! utility::startswith(header,"$$$$"); std::getline(filein,header) ) {
		if( header.size() == 0 || header[0] != '>' ) continue;
		utility::trim( header, " <>\t" ); // Remove tabs, spaces and angle brackets from front and back.
		if( TR.Trace.visible() ) {
			TR.Trace << "Parsing header " << header << std::endl;
		}
		// The data entry ends with a single newline (which getline() will strip out)
		for( std::getline(filein,line); filein && line.size() && ! utility::startswith(line,"$$$$"); std::getline(filein,line) ) {
			if( TR.Trace.visible() ) {
				TR.Trace << "Adding data '" << line << "' (" << line.size() << " chars)" << std::endl;
			}
			if( entry.size() ) {
				entry.push_back('\n'); // Append the newline that getline stripped off.
			}
			entry.append(line);
		}
		molecule.add_str_str_data(header, entry);
		header.clear();
		entry.clear();
		////////// Done processing single data entry.
	}
}

/// @brief Ignore everything up to and including the next '$$$$' entry delimeter.
void
SDFParser::eat_until_delimiter( std::istream & filein ) const {
	std::string line;

	std::getline(filein,line);
	while( filein && !utility::startswith(line,"$$$$") ) {
		std::getline(filein,line);
	}
}


} // sdf
} // io
} // core
