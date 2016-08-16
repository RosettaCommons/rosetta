// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/sdf/single_mol_loading.cxxtest.hh
/// @brief  test suite for loading a single ResidueType from a large Molfile
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Program headers
#include <core/chemical/residue_io.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/sdf/MolFileIOReader.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/database/open.hh>

// Utility headers
#include <utility/io/izstream.hh>

// ObjexxFCL headers

// C++ headers
#include <string>

static basic::Tracer TR("core.io.pdb.file_data_fixup.cxxtest");

using namespace core;


class single_mol_loading_Tests : public CxxTest::TestSuite
{

public:
	// Shared initialization goes here.
	void setUp() {
		core_init();

		using namespace core::chemical;

		ChemicalManager * cm(ChemicalManager::get_instance());
		std::string const tag(FA_STANDARD);
		AtomTypeSetCAP atom_types = cm->atom_type_set(tag);
		ElementSetCAP element_types = cm->element_set("default");
		MMAtomTypeSetCAP mm_atom_types = cm->mm_atom_type_set(tag);
		orbitals::OrbitalTypeSetCAP orbital_types = cm->orbital_type_set(tag);

	}

	/// @brief
	void test_file_positioning() {
		std::string components_file( "core/chemical/sdf/Components_trimmed.sdf" );
		utility::io::izstream filestream( components_file );

		std::map< std::string, core::Size > location_map; // Number of lines to discard before the item
		std::string line;
		bool name_next_line = true;
		core::Size n_lines = 0;
		while ( getline( filestream, line ) ) {
			if ( name_next_line ) {
				std::string name( utility::strip_whitespace( line ) );
				if ( name.size() <= 3 ) { // The database file has some modified standard residues with more than 3 letters
					TR << name << " : " << n_lines << std::endl;
					location_map[ name ] = n_lines;
				}
				name_next_line = false;
			}
			++n_lines;
			if ( utility::startswith(line, "$$$$") ) {
				name_next_line = true;
			}
		}

		utility::io::izstream filestream2( components_file );
		core::chemical::sdf::MolFileIOReader sdf_reader;
		TR << "Position for WOW: " << location_map["WOW"] << std::endl;
		// "Quickly" advance to the appropriate location
		for ( core::Size ii(1); ii <= location_map["WOW"]; ++ii ) {
			getline( filestream2, line );
		}
		TR << ">>>>WOW<<<" << std::endl;
		utility::vector1< core::chemical::sdf::MolFileIOMoleculeOP > entries( sdf_reader.parse_file( filestream2, "sdf", 1 ) );
		runtime_assert( entries.size() == 1 );
		core::chemical::ResidueTypeOP azt_res( core::chemical::sdf::convert_to_ResidueType( entries ) );
		TR << azt_res->name() << std::endl;
		TR << azt_res->name3() << std::endl;
		TR << azt_res->name1() << std::endl;
		TR << "Natoms " << azt_res->natoms() << " nheavy " << azt_res->nheavyatoms() << std::endl;
		azt_res->show_all_atom_names( TR );
	}

};

