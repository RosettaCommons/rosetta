// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/chemical/sdf/CtabV2000Parser.cc
///
/// @brief
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <core/chemical/sdf/CtabV2000Parser.hh>


#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/numbers.hh>

#include <string>
#include <istream>

//// Boost Headers
//#include <boost/foreach.hpp>
//#define foreach BOOST_FOREACH

namespace core {
namespace chemical {
namespace sdf {

static THREAD_LOCAL basic::Tracer TR( "core.chemical.sdf.CtabV2000Parser" );


bool CtabV2000Parser::parse(std::istream & tablein, std::string const & headerline, MolFileIOMolecule & molecule) {

	///////////////// HEADER/Counts line
	// aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
	// 0123456789012345678901234567890123456789
	// 0         1         2         3
	if ( headerline.size() < 39 ) {
		TR.Error << "Counts line for Ctab in mol/sdf file too short." << std::endl;
		return false;
	}

	core::Size const natoms( utility::string2Size( headerline.substr(0,3) ) );
	core::Size const nbonds( utility::string2Size( headerline.substr(3,6) ) );
	//std::string const natomlists(line,6,3);
	//std::string const chiral(line,12,3);
	//std::string const  nstext(line,15,3);
	//std::string const nmlines(line,30,3); // Typically ignored.
	//std::string const version(line,33,6);

	//Sometimes V2000 don't have appropriate version numbers.
	//debug_assert( std::string(line,33,6) == "V2000");

	if ( utility::is_undefined(natoms) || utility::is_undefined(nbonds) ) {
		TR.Error << "Could not read the number of atoms and bonds from header of mol/sdf Ctab." << std::endl;
		return false;
	}

	std::string line;

	///////////////// ATOMS
	for ( core::Size aa(1); aa <= natoms; ++aa ) {
		MolFileIOAtomOP atom( new MolFileIOAtom );
		std::getline(tablein, line);
		if ( !parse_atom_line( line, *atom) ) { return false; }
		atom->index(aa);
		molecule.add_atom(atom);
	}
	///////////////// BONDS
	for ( core::Size bb(1); bb <= nbonds; ++bb ) {
		MolFileIOBondOP bond( new MolFileIOBond );
		std::getline(tablein, line);
		if ( !parse_bond_line( line, *bond) ) { return false; }
		bond->index(bb);
		molecule.add_bond(bond);
	}
	// ATOM LIST - We'd parse the atom list block here, if we had any use for it
	// STEXT - We'd parse the stext block here, if we had any use for it.
	//////////////// PROPERTIES
	for ( std::getline(tablein, line);
			tablein && ! utility::startswith( line, "M  END" );
			std::getline(tablein, line) ) {
		if ( ! parse_property_line( line, molecule ) ) { return false; }
	}
	/// The SDF file data block (with '>' headers) is handled in SDFParser,
	/// as it is not V2000/V3000 specific.
	return true;
}

bool CtabV2000Parser::parse_atom_line( std::string line, MolFileIOAtom & atom) {
	///////////////// ATOM line
	// xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
	// 0123456789012345678901234567890123456789012345678901234567890123456789
	// 0         1         2         3         4         5         6
	if ( line.size() < 34 ) { // 34 as we're only currently interested in the coordinates and the symbol
		TR.Error << "Atom line in mol/sdf Ctab too short." << std::endl;
		return false;
	}

	core::Real const x( utility::string2Real( line.substr(0,10) ) );
	core::Real const y( utility::string2Real( line.substr(10,10) ) );
	core::Real const z( utility::string2Real( line.substr(20,10) ) );
	std::string const symbol(line,31,3);
	//std::string const mass_difference(line,34,2);
	//std::string const charge(line,36,3);
	//std::string const stereo(line,39,3);
	//std::string const hydrogens(line,42,3);
	//std::string const bond_stereo(line,45,3);
	//std::string const valence(line,48,3);
	//std::string const rxn_mapping(line,60,3);
	//std::string const inversion_flag(line,63,3);
	//std::string const inversion_flag(line,66,3);

	if ( utility::is_undefined(x) || utility::is_undefined(y) || utility::is_undefined(z) ) {
		TR.Error << "Cannot read coordinates for atom in mol/sdf Ctab" << std::endl;
		return false;
	}
	atom.position( Vector(x,y,z) );
	atom.element( utility::trim(symbol) );
	return true;
}

bool CtabV2000Parser::parse_bond_line( std::string line, MolFileIOBond & bond) {
	///////////////// Bond line
	// 111222tttsssxxxrrrccc
	// 012345678901234567890
	// 0         1         2
	if ( line.size() < 9 ) { // 34 as we're only currently interested in the atoms and type
		TR.Error << "Bond line in mol/sdf Ctab too short." << std::endl;
		return false;
	}

	core::Size const a1( utility::string2Size( line.substr(0,3) ) );
	core::Size const a2( utility::string2Size( line.substr(3,3) ) );
	core::Size const type( (core::Size)(utility::string2Real( line.substr(6,3) )) );
	//std::string const stereo(line,9,3);
	//std::string const topology(line,15,3);
	//std::string const reacting_center(line,18,3);
	if ( utility::is_undefined(a1) || utility::is_undefined(a2) || utility::is_undefined(type) ) {
		TR.Error << "Cannot read bond in mol/sdf Ctab" << std::endl;
		return false;
	}
	bond.atom1(a1);
	bond.atom2(a2);
	bond.sdf_type(type);
	return true;
}

bool CtabV2000Parser::parse_property_line( std::string line, MolFileIOMolecule & mol) {
	// Because we're ignoring other blocks, we should just silent ignore lines we don't understand,
	// rather than choking on malformated ones.
	if ( utility::startswith(line,"M  CHG") ) {
		/// M  CHGnn8 aaa vvv aaa vvv aaa vvv
		// 0123456789012345678901234567890123456789
		// 0         1         2         3
		core::Size nrecords( utility::string2Size( line.substr(6,3) ) );
		if ( utility::is_undefined(nrecords) ) { // We could check if this is greater than 8, but we don't need to.
			TR.Error << "Malformed CHG record in mol/sdf file" << std::endl;
			return false;
		}
		for ( core::Size ii(1); ii <= nrecords; ++ii ) {
			core::Size const atomi( utility::string2Size( line.substr(8*ii+2,3) ) );
			core::Real const charge( utility::string2Real( line.substr(8*ii+6,3) ) );
			if ( utility::is_undefined(atomi) || utility::is_undefined(charge) ) {
				TR.Error << "Malformed CHG record in mol/sdf file" << std::endl;
				return false;
			}
			MolFileIOAtomOP atom( mol.atom_index(atomi) );
			if ( !atom ) {
				TR.Error << "CHG record in mol/sdf file refers to non-existant atom." << std::endl;
				return false;
			} else {
				atom->formal_charge((int)charge);
			}
		}
	}
	return true;
}

} // sdf
} // io
} // core
