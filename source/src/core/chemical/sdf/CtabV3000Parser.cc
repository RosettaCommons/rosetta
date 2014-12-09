// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/chemical/sdf/CtabV3000Parser.cc
///
/// @brief
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <core/chemical/sdf/CtabV3000Parser.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/numbers.hh>

#include <string>
#include <sstream>
#include <istream>

//// Boost Headers
//#include <boost/foreach.hpp>
//#define foreach BOOST_FOREACH

namespace core {
namespace chemical {
namespace sdf {

static thread_local basic::Tracer TR( "core.chemical.sdf.CtabV3000Parser" );


bool CtabV3000Parser::parse(std::istream & tablein, std::string const & headerline, MolFileIOMolecule & molecule) {

	assert( headerline.compare(33,5,"V3000") == 0 );

	std::string line;
	std::string M, V30, entry, tag;
	core::Size natoms(0),nbonds(0); //,nsg(0),n3d(0),chiral(0);
	for( getline(tablein,line); tablein; getline(tablein,line) ) {
		std::stringstream lstream( line );
		lstream >> M >> V30;
		if( V30 == "END" ) { break; } // The "M  END" line ends the Ctab
		lstream >> entry;
		if( entry == "COUNTS" ) {\
			// Parse counts line
			lstream >> natoms >> nbonds; // >> nsg >> n3d >> chiral;
			// Should we be checking these values? We read them in but never use them.
			continue;
		}
		// The blocks we're currently interested all start with BEGIN - ignore lines until we get there
		if( entry != "BEGIN") { continue; }
		lstream >> tag;
		if( tag == "ATOM" ) {
			std::string line2;
			for( getline(tablein,line2); tablein; getline(tablein,line2) ) {
				if( line2.compare(7,3,"END") == 0 ) {
					break;
				}
				MolFileIOAtomOP atom( new MolFileIOAtom );
				if( !parse_atom_line( line, *atom) ) { return false; }
				molecule.add_atom(atom);
			}
		} else if( tag == "BOND" ) {
			std::string line2;
			for( getline(tablein,line2); tablein; getline(tablein,line2) ) {
				if( line2.compare(7,3,"END") == 0 ) {
					break;
				}
				MolFileIOBondOP bond( new MolFileIOBond );
				if( !parse_bond_line( line, *bond) ) { return false; }
				molecule.add_bond(bond);
			}

		} else {;} // Ignore the other blocks - including the CTAB BEGIN statement.
	}
	/// The SDF file data block (with '>' headers) is handled in SDFParser,
	/// as it is not V2000/V3000 specific.
	return true;
}

bool CtabV3000Parser::parse_atom_line( std::string line, MolFileIOAtom & atom) {
	std::stringstream lstream( line );
	std::string M, V30, type, kvpair;
	core::Size index, aamap;
	core::Real x,y,z;

	lstream >> M >> V30 >> index >> type >> x >> y >> z >> aamap;
	if( !lstream ) {
		TR.Error << "Error reading Atom line in sdf/mol file." << std::endl;
		return false;
	}
	atom.index( index );
	atom.position( Vector(x,y,z) );
	atom.element( type );
	//Now we read the key/value pairs.
	for( lstream >> kvpair; lstream; lstream >> kvpair ) {
		std::string key,value;
		splitkv(kvpair,key,value);
		if( key == "CHG") {
			core::Real charge( utility::string2Real(value));
			if( utility::is_undefined(charge) ) {
				TR.Warning << "Problem parsing formal charge in sdf/mol file, skipping." << std::endl;
				continue;
			}
			atom.formal_charge( (int)charge );
		}
		// Put other properties here -- we don't care about them at the moment, though.
	}
	return true;
}

bool CtabV3000Parser::parse_bond_line( std::string line, MolFileIOBond & bond) {
	std::stringstream lstream( line );
	std::string M, V30, kvpair;
	core::Size type, atom1, atom2, index;

	lstream >> M >> V30 >> index >> type >> atom1 >> atom2;
	if( !lstream ) {
		TR.Error << "Error reading Bond line in sdf/mol file." << std::endl;
		return false;
	}
	bond.index( index );
	bond.atom1( atom1 );
	bond.atom2( atom2 );
	bond.sdf_type( type );
	//Now we read the key/value pairs -- if we were to care about any of them
//	for( lstream >> kvpair; lstream; lstream >> kvpair ) {
//		std::string key,value;
//		splitkv(kvpair,key,value);
//		if( key == "CFG") {
//			core::Size config( utility::string2Real(value));
//			if( utility::is_undefined(config) ) {
//				TR.Warning << "Problem parsing bod configuration in sdf/mol file, skipping." << std::endl;
//				continue;
//			}
//			bond.configuration( config );
//		}
//		// Put other properties here -- we don't care about them at the moment, though.
//	}
	return true;
}

/// @brief Having a '-' at the end of the line is a continuation character for V3000 sdf files.
void
CtabV3000Parser::getline(std::istream & istream, std::string line) const {
	std::getline(istream,line);
	char lastchar;
	while ( istream ) {
		while( lastchar = line[ line.size()-1 ], lastchar == '\n' || lastchar == '\r' || lastchar == ' ' ) {
			line.erase(line.size()-1); //Erase last character to end of string.
		}
		lastchar = line[ line.size()-1 ];
		if( lastchar == '-' ) {
			line.erase(line.size()-1); //Erase '-' from end of string.
			std::string newline;
			std::getline(istream,newline);
			line.append(newline.begin()+6,newline.end()); // append everything but the "M  V30"
		} else {
			// Done!
			return;
		}
	}
}

void
CtabV3000Parser::splitkv(std::string const & kvpair, std::string & key, std::string & value) const {
	core::Size location( kvpair.find('=') );
	if( location == std::string::npos ) { // Not found
		key = kvpair;
		value.clear();
	}
	key.clear();
	value.clear();
	key.append(kvpair,0,location);
	value.append(kvpair,location+1,std::string::npos); //To end of string
}

} // sdf
} // io
} // core
