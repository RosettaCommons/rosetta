// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief  Class to write kinemage-formatted output for Residue and chemical
/// @file   core/chemical/ResidueTypeKinWriter.cc
/// @author Rocco Moretti (rmorettiase@gmail.com)


// Unit headers
#include <core/chemical/ResidueTypeKinWriter.hh>

// Package headers
#include <core/chemical/AtomType.hh>
#include <core/chemical/MMAtomType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/MutableResidueType.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/Atom.hh>

#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

// C++ headers
#include <iostream>
#include <string>

namespace core {
namespace chemical {

std::string
get_element_color(std::string const & element) {
	std::string ELE( element );
	ObjexxFCL::uppercase(ELE);
	if ( ELE == "C" ) return "green";
	if ( ELE == "N" ) return "sky";
	if ( ELE == "O" ) return "red";
	if ( ELE == "H" ) return "gray";
	if ( ELE == "S" ) return "yellow";
	if ( ELE == "P" ) return "peach";
	if ( ELE == "F" ) return "bluetint";
	if ( ELE == "CL" ) return "cyan";
	if ( ELE == "BR" ) return "sea";
	if ( ELE == "I" ) return "lilac";
	return "hotpink";
}


ResidueTypeKinWriter::ResidueTypeKinWriter() = default;

ResidueTypeKinWriter::~ResidueTypeKinWriter() = default;

void ResidueTypeKinWriter::write_kin_header(
	std::ostream & ostr,
	core::chemical::ResidueType const & restype,
	core::Size which_kin
) const
{
	ostr << "@text\n";
	ostr << "View this file with KiNG or Mage from http://kinemage.biochem.duke.edu\n";
	ostr << "@kinemage {" << which_kin << "}\n";
	ostr << "@title { " << restype.name() << "(" <<restype.name3() << ":" << restype.name1() << ") }\n";
	ostr << "@onewidth\n";
	//ostr << "@1center " << ctr.x() << " " << ctr.y() << " " << ctr.z() << "\n";
	//ostr << "@1span 25\n";
}


void
ResidueTypeKinWriter::write_restype(
	std::ostream & ostr,
	core::chemical::ResidueType const & restype,
	core::Size which_kin /*= 1*/ // multiple kinemages can be put in a single kinemage file.
) const
{
	core::Size natoms( restype.natoms() );

	write_kin_header(ostr, restype, which_kin );

	ostr << "@balllist {element balls} color= gray radius= 0.1\n";
	for ( core::Size ii(1); ii <= natoms; ++ii ) {
		core::Vector const & pos( restype.ideal_xyz(ii) );
		core::chemical::AtomType const & at( restype.atom_type(ii) );

		ostr << "{" << at.element() << "} " << get_element_color(at.element()) << " "
			<< pos.x() << " " << pos.y() << " " << pos.z() << "\n";
	}

	ostr << "@vectorlist {all bonds} color= gray width= 1\n";
	// To get bonds, we iterate over atoms, looking for the bonds they're attached to (with ordering so we don't double count.
	// It would be nice to be able to annotate these bonds with the bond type.
	// (Would involve drawing multiple lines, offsetting, etc.)
	for ( core::Size ii(1); ii <= natoms; ++ii ) {
		core::Vector const & pos1( restype.ideal_xyz(ii) );
		core::chemical::AtomIndices const & neighbors( restype.bonded_neighbor( ii ) );
		for ( core::Size jj: neighbors ) {
			if ( ii >= jj ) { continue; }
			core::Vector const & pos2( restype.ideal_xyz(jj) );
			ostr << "{" << restype.atom_name(ii) << "}P "
				<< pos1.x() << " " << pos1.y() << " " << pos1.z() << "\n";
			ostr << "{" << restype.atom_name(jj) << "}L "
				<< pos2.x() << " " << pos2.y() << " " << pos2.z() << "\n";
		}
	}

	ostr << "@vectorlist {rotatable bonds} color= white width= 4\n";
	for ( core::Size ii(1); ii <= restype.nchi(); ++ii ) {
		AtomIndices const & chi( restype.chi_atoms( ii ) );
		debug_assert( chi.size() == 4 );
		core::Vector const & pos2( restype.ideal_xyz(chi[2]) );
		core::Vector const & pos3( restype.ideal_xyz(chi[3]) );
		ostr << "{" << restype.atom_name(chi[2]) << "}P "
			<< pos2.x() << " " << pos2.y() << " " << pos2.z() << "\n";
		ostr << "{" << restype.atom_name(chi[3]) << "}L "
			<< pos3.x() << " " << pos3.y() << " " << pos3.z() << "\n";
	}

	ostr << "@labellist {atom indices} color= white off\n";
	for ( core::Size ii(1); ii <= natoms; ++ii ) {
		core::Vector const & pos( restype.ideal_xyz(ii) );

		ostr << "{" << ii << "} "
			<< pos.x() << " " << pos.y() << " " << pos.z() << "\n";
	}

	ostr << "@labellist {atom names} color= white off\n";
	for ( core::Size ii(1); ii <= natoms; ++ii ) {
		core::Vector const & pos( restype.ideal_xyz(ii) );

		ostr << "{" << restype.atom_name(ii) << "} "
			<< pos.x() << " " << pos.y() << " " << pos.z() << "\n";
	}

	ostr << "@labellist {Rosetta types} color= white\n";
	for ( core::Size ii(1); ii <= natoms; ++ii ) {
		core::Vector const & pos( restype.ideal_xyz(ii) );
		core::chemical::AtomType const & at( restype.atom_type( ii ) );

		ostr << "{" << at.name() << "} "
			<< pos.x() << " " << pos.y() << " " << pos.z() << "\n";
	}

	ostr << "@labellist {MM types} color= white off\n";
	for ( core::Size ii(1); ii <= natoms; ++ii ) {
		core::Vector const & pos( restype.ideal_xyz(ii) );
		MMAtomType const & mm( restype.mm_name( ii ) );

		ostr << "{" << mm.name() << "} "
			<< pos.x() << " " << pos.y() << " " << pos.z() << "\n";
	}

	ostr << "@labellist {partial charges} color= white off\n";
	std::ios::fmtflags oldflags( ostr.flags() );
	std::streamsize oldprec( ostr.precision(2) );
	ostr.setf( std::ios::fixed );
	for ( core::Size ii(1); ii <= natoms; ++ii ) {
		core::Vector const & pos( restype.ideal_xyz(ii) );

		ostr << "{" << restype.atom_charge(ii) << "} "
			<< pos.x() << " " << pos.y() << " " << pos.z() << "\n";
	}
	ostr.flags( oldflags );
	ostr.precision( oldprec );

	// Fragment sets

	ostr << "@balllist {root atom} color= purple radius= 0.2 \n";
	// The default root atom for jumps is always 1 for ligands.
	// For residues with backbones it's a little more complicated.
	// See core/conformation/util.cc:904 for details
	core::Size root_atom(1);
	if ( restype.mainchain_atoms().size() > 0 ) {
		root_atom = restype.mainchain_atoms()[ (restype.mainchain_atoms().size()-1)/2 + 1 ];
	}

	{ // For scope limitation
		core::Vector const & pos( restype.ideal_xyz(root_atom) );
		ostr << "{" << restype.atom_name(root_atom) << "} "
			<< pos.x() << " " << pos.y() << " " << pos.z() << "\n";
	}

	{ // For scope limitation
		ostr << "@ringlist {nbr atom} color= purple radius= " << restype.nbr_radius() << " off\n";
		core::Size nbr( restype.nbr_atom() );
		core::Vector const & pos( restype.ideal_xyz(nbr) );
		ostr << "{" << restype.atom_name(nbr) << "} "
			<< pos.x() << " " << pos.y() << " " << pos.z() << "\n";
		ostr << "{" << restype.atom_name(nbr) << "} r=0.3 "
			<< pos.x() << " " << pos.y() << " " << pos.z() << "\n";
	}


	ostr << "@arrowlist {atom tree} color= purple \n";
	for ( core::Size ii(1); ii <= natoms; ++ii ) {
		if ( ii == root_atom ) { continue; }

		core::Vector const & pos1( restype.ideal_xyz(ii) );
		core::Vector const & pos2( restype.ideal_xyz( restype.atom_base(ii) ) );

		ostr << "{" << restype.atom_name( restype.atom_base(ii) ) << "}P "
			<< pos2.x() << " " << pos2.y() << " " << pos2.z() << "\n";
		ostr << "{" << restype.atom_name(ii) << "}L "
			<< pos1.x() << " " << pos1.y() << " " << pos1.z() << "\n";
	}
}

} // chemical
} // core

