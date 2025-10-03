// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/chemical/sdf/mmCIFWriter.cc
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <core/chemical/mmCIF/mmCIFWriter.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/Element.hh>

#include <core/types.hh>

//Utility functions
#include <basic/Tracer.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
#include <utility/stream_util.hh>
#include <utility/gemmi_util.hh>

//external CIF includes
#include <gemmi/cif.hpp>
#include <gemmi/to_cif.hpp>
#include <gemmi/numb.hpp> // for as_number

namespace core {
namespace chemical {
namespace mmCIF {


//Load up the tracer for this class
static basic::Tracer TR( "core.io.mmCIF.mmCIFWriter" );

void
mmCIFWriter::write_file(std::string const & file_name, core::chemical::ResidueType const & restype) {
	utility::io::ozstream outfile;
	outfile.open(file_name.c_str(), std::ios::out);
	if ( !outfile ) {
		throw CREATE_EXCEPTION(utility::excn::FileNotFound, "Cannot open file"+file_name);
	}
	write_stream(outfile,restype);
	outfile.close();
}

void
mmCIFWriter::write_stream(std::ostream & output_stream, core::chemical::ResidueType const & restype) {
	gemmi::cif::WriteOptions options;
	options.prefer_pairs = true;
	options.misuse_hash = true;

	gemmi::cif::write_cif_block_to_stream( output_stream, generate_block(restype), options );
}

gemmi::cif::Block
mmCIFWriter::generate_block( core::chemical::ResidueType const & restype ) {
	using utility::gemmi_add_table;
	using utility::gemmi_add_row;

	gemmi::cif::Block block;
	block.name = restype.name();

	block.set_pair( "_chem_comp.id", restype.name() );
	block.set_pair( "_chem_comp.name", restype.name() );
	block.set_pair( "_chem_comp.one_letter_code", std::string(1, restype.name1() ) );
	block.set_pair( "_chem_comp.three_letter_code", restype.name3() );
	if ( restype.is_polymer() ) {
		if ( restype.is_l_aa() ) {
			block.set_pair( "_chem_comp.type", "L-PEPTIDE LINKING" );
		} else if ( restype.is_d_aa() ) {
			block.set_pair( "_chem_comp.type", "D-PEPTIDE LINKING" );
		} else if ( restype.is_RNA() ) {
			block.set_pair( "_chem_comp.type", "RNA LINKING" );
		} else if ( restype.is_DNA() ) {
			block.set_pair( "_chem_comp.type", "DNA LINKING" );
		}
	} else {
		block.set_pair( "_chem_comp.type", "NON-POLYMER" );
	}


	gemmi::cif::Loop & atom_comp = gemmi_add_table(block, "_chem_comp_atom", {
		"comp_id",
		"atom_id",
		"type_symbol",
		"charge",
		"model_Cartn_x",
		"model_Cartn_y",
		"model_Cartn_z",
		"pdbx_ordinal",
		} );

	for ( core::Size ii(1); ii <= restype.natoms(); ++ii ) {
		std::vector< std::string > vec;
		vec.push_back( restype.name() );
		vec.push_back( restype.atom_name(ii) );
		vec.push_back( restype.element_type(ii)->get_chemical_symbol() );
		vec.push_back( std::to_string( restype.formal_charge(ii) ) );
		core::Vector pos = restype.ideal_xyz(ii);
		vec.push_back( std::to_string(pos.x()) );
		vec.push_back( std::to_string(pos.y()) );
		vec.push_back( std::to_string(pos.z()) );
		vec.push_back( std::to_string(ii) );

		gemmi_add_row( atom_comp, vec );
	}

	gemmi::cif::Loop & bond_comp = gemmi_add_table(block, "_chem_comp_bond", {
		"comp_id",
		"atom_id_1",
		"atom_id_2",
		"value_order",
		"pdbx_aromatic_flag",
		"pdbx_ordinal"
		} );

	utility::vector1< std::pair< core::Size, core::Size > > const & all_bonds = restype.bonds();

	for ( core::Size ii(1); ii <= restype.nbonds(); ++ii ) {
		std::vector< std::string > vec;
		vec.push_back( restype.name() );
		core::Size atm1 = all_bonds[ii].first;
		core::Size atm2 = all_bonds[ii].second;
		vec.push_back( restype.atom_name( atm1 ) );
		vec.push_back( restype.atom_name( atm2 ) );
		std::string bond_type = "UNK";
		switch ( restype.bond_type( atm1, atm2 ) ) {
			case SingleBond:
				bond_type = "SING"; break;
			case DoubleBond:
				bond_type = "DOUB"; break;
			case TripleBond:
				bond_type = "TRIP"; break;
			case AromaticBond:
				bond_type = "AROM"; break;
			default:
				bond_type = "UNK"; break;
		}
		vec.push_back( bond_type );
		vec.push_back( (restype.bond_type( atm1, atm2 ) == AromaticBond) ? "Y" : "N" );
		vec.push_back( std::to_string(ii) );

		gemmi_add_row( bond_comp, vec );
	}

	return block;
}


}
}
}
