// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/rmoretti/restype_converter.cc
/// @brief A utility for converting restype restype files (params, sdf, pdb)
/// @author Rocco Moretti (rmorettiase@gmail.com)

// devel headers
#include <devel/init.hh>

// core headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/residue_io.hh>
#include <core/chemical/residue_support.hh>
#include <core/chemical/sdf/MolFileIOReader.hh>
#include <core/chemical/sdf/mol_writer.hh>
#include <core/chemical/mmCIF/mmCIFParser.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pose/Pose.hh>


// utility headers
#include <utility/excn/Exceptions.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <utility/options/OptionCollection.hh>
#include <basic/options/option_macros.hh>

OPT_1GRP_KEY( StringVector, restype_convert, types )
OPT_1GRP_KEY( StringVector, restype_convert, name3 )

OPT_1GRP_KEY( Boolean, restype_convert, params_out )
OPT_1GRP_KEY( Boolean, restype_convert, sdf_out )

static basic::Tracer TR("restype_converter");

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	option.add_relevant( in::file::centroid );
	option.add_relevant( restype_convert::types );
	option.add_relevant( restype_convert::name3 );
	option.add_relevant( in::file::extra_res_cen );
	option.add_relevant( in::file::extra_res_fa );
	option.add_relevant( in::file::extra_res_mol );
	option.add_relevant( in::file::extra_res_mmCIF );

	option.add_relevant( out::path::all );
	option.add_relevant( out::prefix );
	option.add_relevant( out::suffix );

	option.add_relevant( out::pdb );
	option.add_relevant( restype_convert::params_out );
	option.add_relevant( restype_convert::sdf_out );

}

/// Some filetypes need to be loaded as full atom, and then converted appropriately.
utility::vector1< core::chemical::ResidueTypeCOP >
load_as_fullatom() {
	using namespace core::chemical;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::vector1< core::chemical::ResidueTypeCOP > fullatom;

	if ( ! option[OptionKeys::in::file::extra_res_mol].active() && ! option[OptionKeys::in::file::extra_res_mmCIF].active() ) {
		return fullatom;
	}

	core::chemical::AtomTypeSetCOP atom_types = ChemicalManager::get_instance()->atom_type_set("fa_standard");
	core::chemical::ElementSetCOP elements = ChemicalManager::get_instance()->element_set("default");
	core::chemical::MMAtomTypeSetCOP mm_atom_types = ChemicalManager::get_instance()->mm_atom_type_set("fa_standard");

	if ( option[OptionKeys::in::file::extra_res_mol].active() ) {
		sdf::MolFileIOReader molfile_reader;
		for ( auto const & filename: option[OptionKeys::in::file::extra_res_mol]() ) {
			utility::vector1< sdf::MolFileIOMoleculeOP > data( molfile_reader.parse_file( filename ) );
			fullatom.append( sdf::convert_to_ResidueTypes( data, /* load_rotamers= */ true, atom_types, elements, mm_atom_types ) );
		}
	}

	if ( option[OptionKeys::in::file::extra_res_mmCIF].active() ) {
		mmCIF::mmCIFParser mmCIF_parser;
		for ( auto const & filename : option[OptionKeys::in::file::extra_res_mmCIF]() ) {
			utility::vector1< sdf::MolFileIOMoleculeOP> molecules( mmCIF_parser.parse( filename ) );
			fullatom.append( sdf::convert_to_ResidueTypes( molecules, true, atom_types, elements, mm_atom_types ) );
		}
	}

	return fullatom;
}

utility::vector1< core::chemical::ResidueTypeCOP >
collect_residue_types() {
	using namespace core::chemical;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::vector1< core::chemical::ResidueTypeCOP > restypes;

	ResidueTypeSetCOP type_set;
	std::string type_set_name; // For diagnostic output
	if ( option[ in::file::centroid ]() ) {
		type_set_name = "Centroid";
		type_set = ChemicalManager::get_instance()->residue_type_set( CENTROID_t );
	} else {
		type_set_name = "Full atom";
		type_set = ChemicalManager::get_instance()->residue_type_set( FULL_ATOM_t );
	}

	for ( auto const & name: option[ OptionKeys::restype_convert::types ]() ) {
		if ( ! type_set->has_name( name ) ) {
			TR.Error << type_set_name << " type set does not have residue '" << name << "'" << std::endl;
			continue;
		}
		restypes.push_back( type_set->name_mapOP( name ) );
	}

	for ( auto const & name3: option[ restype_convert::name3 ]() ) {
		if ( ! type_set->has_name3( name3 ) ) {
			TR.Error << type_set_name << " type set does not have any residue with the three letter code '" << name3 << "'" << std::endl;
			continue;
		}
		ResidueTypeFinder rtf( *type_set );
		rtf.name3( name3 );
		restypes.append( rtf.get_all_possible_residue_types( true /* allow_extra_variants */ ) ); // Match how PDB loading does it.
	}

	if ( type_set->mode() == FULL_ATOM_t ) {
		for ( auto const & filename: option[OptionKeys::in::file::extra_res_fa]() ) {
			restypes.push_back( read_topology_file( filename, type_set ) );
		}
	}

	if ( type_set->mode() == CENTROID_t ) {
		for ( auto const & filename: option[OptionKeys::in::file::extra_res_cen]() ) {
			restypes.push_back( read_topology_file( filename, type_set ) );
		}
	}

	for ( core::chemical::ResidueTypeCOP fa_type: load_as_fullatom() ) {
		if ( type_set->mode() == CENTROID_t ) {
			ResidueTypeOP centroid_type( make_centroid( *fa_type ) );
			restypes.push_back( centroid_type );
		} else {
			restypes.push_back( fa_type );
		}
	}

	return restypes;
}

std::string
determine_output_name( std::string const & resname, std::string const & extension ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::string outname;
	if ( option[ out::path::pdb ].user() && extension == "pdb" ) {
		outname = option[ out::path::pdb ]().path();
	} else if ( option[ out::path::all ].user() ) {
		outname = option[ out::path::all ]().path();
	}
	outname += option[ out::prefix ]();
	outname += resname;
	outname += option[ out::suffix ]();
	outname += '.';
	outname += extension;
	return outname;
}


void
output_residue_types( utility::vector1< core::chemical::ResidueTypeCOP > const & restypes ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;

	for ( ResidueTypeCOP restype: restypes ) {
		TR << "Outputting residue type '" << restype->name() << "' with three letter code '" << restype->name3() << "'" << std::endl;

		std::string const & name = restype->name();

		if ( option[ out::pdb ]() ) {
			core::conformation::ResidueOP res = core::conformation::ResidueFactory::create_residue( *restype );
			core::pose::Pose pose;
			pose.append_residue_by_bond( *res );

			pose.dump_pdb( determine_output_name( name, "pdb" ) );
		}

		if ( option[ restype_convert::params_out ]() ) {
			write_topology_file( *restype, determine_output_name( name, "params" ) );
		}

		if ( option[ restype_convert::sdf_out ]() ) {
			sdf::MolWriter writer;
			writer.output_residue( determine_output_name( name, "sdf" ), restype );
		}

	}
}

int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		NEW_OPT( restype_convert::types, "Which types from the standard database to convert, full name", utility::vector1<std::string>{} );
		NEW_OPT( restype_convert::name3, "Which types from the standard database to output, three letter codes", utility::vector1<std::string>{} );
		NEW_OPT( restype_convert::params_out, "Output params files as output", false );
		NEW_OPT( restype_convert::sdf_out, "Output sdf files as output", false );

		devel::init( argc, argv );
		register_options();

		utility::vector1< core::chemical::ResidueTypeCOP > restypes = collect_residue_types();
		output_residue_types( restypes );

	} catch ( utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;
}
