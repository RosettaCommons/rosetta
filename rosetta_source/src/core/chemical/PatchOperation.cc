// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/chemical/Residue_.cc
/// @brief implementation class for abstract class Residue if you want to see what it actually does, look at comments in Patch.cc
/// @author Phil Bradley
// Unit headers
#include <core/chemical/PatchOperation.hh>

// Package Headers


// ObjexxFCL headers

// Numeric headers
#include <numeric/conversions.hh>

// Utility headers
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


// Commented by inclean daemon #include <ObjexxFCL/string.functions.hh>

// C++ headers
// Commented by inclean daemon #include <sstream>

namespace core {
namespace chemical {

static basic::Tracer tr("core.chemical");

/// helper function
std::string
expand_icoor_atom_name( std::string name, ResidueType const & rsd )
{
	std::string const nconn_tag( "%NCONN" );
	Size pos( name.find( nconn_tag ) );
	if ( pos < name.size() ) {
		//std::cout << "name before replace: " << name << std::endl;
		name.replace( pos, nconn_tag.size(), ObjexxFCL::string_of( rsd.n_residue_connections() ) );
		//std::cout << "name after replace: " << name << std::endl;
	}
	return name;

}

bool
SetICoor::apply( ResidueType & rsd ) const
{
	//std::cout << "SetICoor::apply: " << atom_ << ' ' << stub1_ << ' ' << stub2_ << ' ' << stub3_ <<
	//	std::endl;
	std::string const atom ( expand_icoor_atom_name( atom_ , rsd ) );
	std::string const stub1( expand_icoor_atom_name( stub1_, rsd ) );
	std::string const stub2( expand_icoor_atom_name( stub2_, rsd ) );
	std::string const stub3( expand_icoor_atom_name( stub3_, rsd ) );
	// 	bool const rebuild_icoor_xyz( ICoorAtomID( stub1, rsd ).is_internal() &&
	// 																ICoorAtomID( stub2, rsd ).is_internal() &&
	// 																ICoorAtomID( stub3, rsd ).is_internal() );
	bool const rebuild_icoor_xyz( true );
	rsd.set_icoor( atom, phi_, theta_, d_, stub1, stub2, stub3, rebuild_icoor_xyz );
	return false;
}


PatchOperationOP
patch_operation_from_patch_file_line( std::string const & line ) {
	using numeric::conversions::radians;
	std::istringstream l( line );
	std::string tag, atom1, atom2, atom3, atom4, atom_name, atom_type_name, mm_atom_type_name, property;
	Real charge, mean, sdev;
	Size chino;
	l >> tag;
	if ( l.fail() || tag[0] == '#' ) return 0;
	if ( tag == "ADD_ATOM" ) {
		if ( line.size() < 25 ) return 0;
		atom_name = line.substr( 9,4); l >> tag;
		l >> atom_type_name; // = line.substr( 14,4);
		l >> mm_atom_type_name; // = line.substr( 19,4);
		l >> charge;
		if ( l.fail() ) return 0;
		return new AddAtom( atom_name, atom_type_name, mm_atom_type_name, charge );
	} else if ( tag == "DELETE_ATOM" ) {
		l >> atom_name;
		if ( l.fail() ) return 0;
		return new DeleteAtom( atom_name );
	} else if ( tag == "SET_BACKBONE_HEAVYATOM" ) {
		l >> atom_name;
		if ( l.fail() ) return 0;
		return new SetBackboneHeavyatom( atom_name );
	} else if ( tag == "SET_IO_STRING" ) { // 13 character tag
		// NOTE - USE FIXED WIDTH IO SINCE NAME3 CAN CONTAIN INTERNAL WHITESPACE (EG DNA,RNA)
		if ( line.size() < 19 ) return 0;
		std::string const three_letter_code( line.substr(14,3) ), one_letter_code( line.substr(18,1) );
		return new SetIO_String( three_letter_code, one_letter_code[0] );
	} else if ( tag == "ADD_PROPERTY" ) {
		l >> property;
		if ( l.fail() ) return 0;
		return new AddProperty( property );
	} else if ( tag == "DELETE_PROPERTY" ) {
		l >> property;
		if ( l.fail() ) return 0;
		return new DeleteProperty( property );
		//Added by Andy M. Chen in June 2009
		//  This is needed for deleting properties, which occurs in certain PTM's
	} else if ( tag == "ADD_CHI" ) {
		l >> chino >> atom1 >> atom2 >> atom3 >> atom4;
		if ( l.fail() ) return 0;
		return new AddChi( chino, atom1, atom2, atom3, atom4 );
		//Added by Andy M. Chen in June 2009
		//    This is needed for PTM's, which often result in one or more extra chi angles
	} else if ( tag == "REDEFINE_CHI" ) {
		l >> chino >> atom1 >> atom2 >> atom3 >> atom4;
		if ( l.fail() ) return 0;
		return new RedefineChi( chino, atom1, atom2, atom3, atom4 );
		//Added by Andy M. Chen in June 2009
		//    This is needed for PTM's
	} else if ( tag == "ADD_CHI_ROTAMER" ) {
		l >> chino >> mean >> sdev;
		if ( l.fail() ) return 0;
		return new AddChiRotamer( chino, mean, sdev );
		//Added by Andy M. Chen in June 2009
		//    This is needed for PTM's
	} else if ( tag == "ADD_BOND" ) {
		l >> atom1 >> atom2;
		if ( l.fail() ) return 0;
		return new AddBond( atom1, atom2 );
	} else if ( tag == "ADD_CONNECT" ) {
		std::string connect_atom;
		l >> connect_atom;
		if ( l.fail() ) return 0;
		l >> tag;
		if ( l.fail() ) {
			return new AddConnect( connect_atom, 0.0, 0.0, 0.0, connect_atom, connect_atom, connect_atom );
		} else {
			Real phi, theta, d;
			std::string parent_atom, angle_atom, torsion_atom;
			l >> phi >> theta >> d >> parent_atom >> angle_atom >> torsion_atom;
			if ( l.fail() || tag != "ICOOR" ) {
				utility_exit_with_message( "bad line in patchfile: "+line );
			}
			return new AddConnect( connect_atom, radians(phi), radians(theta), d, parent_atom, angle_atom, torsion_atom );
		}
	} else if ( tag == "SET_ATOM_TYPE" ) {
		l >> atom_name >> atom_type_name;
		if ( l.fail() ) return 0;
		return new SetAtomType( atom_name, atom_type_name );
	} else if ( tag == "SET_MM_ATOM_TYPE" ) {
		l >> atom_name >> mm_atom_type_name;
		if ( l.fail() ) return 0;
		return new SetMMAtomType( atom_name, mm_atom_type_name );
	} else if ( tag == "SET_ATOMIC_CHARGE" ) {
		l >> atom_name >> charge;
		if ( l.fail() ) return 0;
		return new SetAtomicCharge( atom_name, charge );
	} else if ( tag == "SET_POLYMER_CONNECT" ) {
		l >> tag >> atom_name; // tag should be "UPPER" or "LOWER"
		if ( l.fail() ) return 0;
		return new SetPolymerConnectAtom( atom_name, tag );
	} else if ( tag == "SET_ICOOR" ) {
		Real phi,theta,d;
		std::string stub1, stub2, stub3;
		l >> atom_name >> phi >> theta >> d >> stub1 >> stub2 >> stub3;
		if ( l.fail() ) return 0;
		return new SetICoor( atom_name, radians(phi), radians(theta), d, stub1, stub2, stub3 );
	} else if ( tag == "PREPEND_MAINCHAIN_ATOM" ) {
		l >> atom_name;
		return new PrependMainchainAtom( atom_name );
	} else if ( tag == "APPEND_MAINCHAIN_ATOM" ) {
		l >> atom_name;
		return new AppendMainchainAtom( atom_name );
	} else if ( tag == "NCAA_ROTLIB_PATH" ) {
		std::string path;
		l >> path;
		return new NCAARotLibPath( path );
	}
	tr.Warning << "patch_operation_from_patch_file_line: bad line: " << line << std::endl;

	return 0;
}


} // chemical
} // core
