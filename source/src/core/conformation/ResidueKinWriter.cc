// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief  Class to write kinemage-formatted output for Residue and Conformation
/// @file   core/conformation/ResidueKinWriter.cc
/// @author Andrew Leaver-Fay


// Unit headers
#include <core/conformation/ResidueKinWriter.hh>

// Package headers
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>

// C++ headers
#include <iostream>

#include <utility/vector1.hh>


namespace core {
namespace conformation {

void write_kinemage_header(
	std::ostream & ostr,
	Size const kin_number,
	std::string const & title,
	Vector const & ctr
)
{
	ostr << "@kinemage {" << kin_number << "}\n";
	ostr << "@title { " << title << " }\n";
	ostr << "@1center " << ctr.x() << " " << ctr.y() << " " << ctr.z() << "\n";
	ostr << "@1span 25\n";
}

/// TWO FUNCTIONS STOLEN FROM IAN: and slightly modified.
void print_node(
	std::ostream & out,
	int residue_num,
	int atom_num,
	core::Vector const & atom_xyz,
	core::chemical::ResidueType const & res,
	std::string extras = "" //< P for points, color, width/radius, etc.
)
{
	// atom_num is often 0 in fold tree, means no specific atom.
	// might as well use the first() one:
	if ( atom_num == 0 ) atom_num = 1;
	core::chemical::AtomType const & atom_type = res.atom_type(atom_num);
	// This info appears when you click on the point
	out << "{" << res.name3() << " " << residue_num
		<< " " << res.atom_name(atom_num) << " (" << atom_type.name() << ")"
		<< "}";
	// Color, width, etc. followed by coordinates
	//out << "col" << residue_num << " ";
	out << extras;
	out << " " << atom_xyz.x() << " " << atom_xyz.y() << " " << atom_xyz.z() << "\n";
}

void print_node(
	std::ostream & out,
	int residue_num,
	int atom_num,
	core::conformation::Residue const & res,
	std::string extras = "" //< P for points, color, width/radius, etc.
)
{
	print_node( out, residue_num, atom_num, res.xyz( atom_num ), res.type(), extras );
}


void print_node(
	std::ostream & out,
	int residue_num,
	std::string atom_name,
	core::conformation::Residue const & res,
	std::string extras = "" //< P for points, color, width/radius, etc.
)
{
	// atom_num is often 0 in fold tree, means no specific atom.
	// might as well use the first one:

	int atom_num;
	if ( atom_name == "" ) {
		atom_num = 1;
	} else {
		atom_num = res.atom_index( atom_name );
	}
	print_node( out, residue_num, atom_num, res, extras );
}

ResidueKinWriter::ResidueKinWriter() :
	master_( "" ),
	dominant_( true ),
	animate_( true ),
	group_( true ),
	write_apolar_hydrogens_( false ),
	write_polar_hydrogens_( false ),
	write_backbone_hydrogens_( true ),
	write_virtual_atoms_( false )
{}

ResidueKinWriter::~ResidueKinWriter() {}

void ResidueKinWriter::write_kin_header(
	std::ostream & ostr,
	core::conformation::Residue const & rsd,
	core::Size atom_to_center_on,
	core::Size which_kin
) const
{
	Vector center_point( 0 );
	if ( atom_to_center_on == 0 ) {
		center_point = rsd.xyz( rsd.nbr_atom() );
	} else {
		center_point = rsd.xyz( atom_to_center_on );
	}
	write_kinemage_header( ostr, which_kin, "Residue Kinemage", center_point );
}


/// @details If you're drawing multiple instances of a single rotamer, use the
/// "instance" flag to point to the coordinates already written out
/// in this file.  This creates a smaller output file.
void
ResidueKinWriter::write_rsd_coords(
	std::ostream & ostr,
	core::conformation::Residue const & rsd,
	bool is_instance /* = false */
) const
{

	// intra-residue connections
	// do residues in different (~random) colors to help distinguish them
	//int const num_colors = 6;
	//std::string colors[num_colors] = {"pinktint", "peachtint", "yellowtint", "greentint", "bluetint", "lilactint"};
	//std::string color = colors[ rsd.seqpos() % num_colors ];
	std::string color = "bluetint";

	std::string tag = "";

	ostr << "@" << ( group_ ? "group" : "subgroup" ) << " { Res " << rsd.seqpos() << " }";
	if ( animate_ ) ostr << " animate";
	if ( dominant_ ) ostr << " dominant";
	ostr << "\n";
	ostr << "@vectorlist {";
	if ( is_instance ) {
		ostr << rsd.name() << " " << rsd.seqpos() << "i";
	} else {
		ostr << rsd.name() << " " << rsd.seqpos();
	}
	ostr << "} color= " << color;
	if ( master_ != "" ) ostr << " master= {" << master_ << "}";
	if ( is_instance ) {
		ostr << " instance= {" << rsd.name() << " " << rsd.seqpos() << "}";
	}
	ostr << "\n";

	if ( ! is_instance ) {
		for ( core::Size atom_i = 1; atom_i <= rsd.natoms(); ++atom_i ) {
			core::conformation::Residue::AtomIndices const & nbrs = rsd.nbrs(atom_i);
			for ( core::conformation::Residue::AtomIndices::const_iterator j = nbrs.begin(), end_j = nbrs.end(); j != end_j; ++j ) {
				core::Size atom_j = *j;
				if ( atom_j <= atom_i ) continue; // so we draw each bond just once, not twice
				bool const is_H = rsd.atom_is_hydrogen(atom_j) || rsd.atom_is_hydrogen(atom_i);

				if ( is_H && ! write_apolar_hydrogens_ && ! write_polar_hydrogens_ ) continue;

				if ( ! write_virtual_atoms_ && ( rsd.atom_type( atom_i ).element() == "X" || rsd.atom_type( atom_j ).element() == "X" ) ) continue;

				/// backbone hydrogens?
				if ( ! write_backbone_hydrogens_ && rsd.atom_is_backbone( atom_i ) && rsd.atom_is_hydrogen( atom_i ) ) continue;
				if ( ! write_backbone_hydrogens_ && rsd.atom_is_backbone( atom_j ) && rsd.atom_is_hydrogen( atom_j ) ) continue;
				std::string const ptmaster = ( is_H ? " 'h'" : "" );
				print_node( ostr, rsd.seqpos(), atom_i, rsd, tag+"P"+ptmaster);
				print_node( ostr, rsd.seqpos(), atom_j, rsd, tag+ptmaster);
			}
		}
	}
}


void ResidueKinWriter::dominant( bool setting ) { dominant_ = setting; }
void ResidueKinWriter::animate(  bool setting ) { animate_  = setting; }
void ResidueKinWriter::group(    bool setting ) { group_    = setting; }
void ResidueKinWriter::master( std::string const & setting ) { master_ = setting; }
void ResidueKinWriter::write_virtual_atoms( bool setting ) { write_virtual_atoms_ = setting; }
void ResidueKinWriter::write_hydrogens(          bool setting ) {
	write_apolar_hydrogens( setting );
	write_polar_hydrogens( setting );
	write_backbone_hydrogens( setting );
}

void ResidueKinWriter::write_apolar_hydrogens(   bool setting ) { write_apolar_hydrogens_ = setting; }
void ResidueKinWriter::write_polar_hydrogens(    bool setting ) { write_polar_hydrogens_ = setting; }
void ResidueKinWriter::write_backbone_hydrogens( bool setting ) { write_backbone_hydrogens_ = setting; }

ConformationKinWriter::~ConformationKinWriter() {}

/// @brief Write out the coordinates for an entire conformation; this includes
/// inter-residue bonds that would be missed by the ResidueKinWriter.
void
ConformationKinWriter::write_coords(
	std::ostream & ostr,
	core::conformation::Conformation const & conf,
	bool is_instance
) const
{
	using namespace core::conformation;

	ResidueKinWriter rsd_writer;
	rsd_writer.animate( false );
	rsd_writer.group( false );
	ostr << "@group { conformation  } dominant on\n";

	for ( Size ii = 1; ii <= conf.size(); ++ii ) {
		Residue const & ii_rsd( conf.residue( ii ) );
		rsd_writer.write_rsd_coords( ostr, ii_rsd, is_instance );
		/// insert code here to write inter-residue bonds.

		for ( Size jj = 1; jj <= ii_rsd.n_possible_residue_connections(); ++jj ) {
			Size const jj_conn_residue = ii_rsd.connected_residue_at_resconn( jj );
			if ( jj_conn_residue < ii ) continue; // we've already output this connection
			Residue const & jj_rsd( conf.residue( jj_conn_residue ) );
			Size const ii_conn_atom = ii_rsd.residue_connect_atom_index( jj );
			Size const jj_conn_id   = ii_rsd.connect_map( jj ).connid();
			Size const jj_conn_atom = jj_rsd.residue_connect_atom_index( jj_conn_id );

			print_node( ostr, ii_rsd.seqpos(), ii_conn_atom, ii_rsd, "P");
			print_node( ostr, jj_rsd.seqpos(), jj_conn_atom, jj_rsd, "");

		}
	}
}

void ConformationKinWriter::write_virtual_atoms( bool setting ) { write_virtual_atoms_ = setting; }
void ConformationKinWriter::master( std::string const & setting ) { master_ = setting; }


} // conformation
} // core

