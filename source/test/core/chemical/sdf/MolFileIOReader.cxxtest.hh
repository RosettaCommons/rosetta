// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/chemical/sdf/MolfileIOReader.cxxtest.hh
/// @brief unit tests for the MolfileIOReader file
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

#include <core/chemical/sdf/MolFileIOData.hh>
#include <core/chemical/sdf/MolFileIOReader.hh>

// Unit Headers
#include <core/chemical/ResidueType.hh>

#include <core/chemical/rotamers/RotamerLibrarySpecification.hh>
#include <core/chemical/rotamers/StoredRotamerLibrarySpecification.hh>

#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>
#include <core/pack/rotamers/SingleLigandRotamerLibrary.hh>

#include <core/chemical/AtomICoor.hh>

// Project Headers
#include <core/types.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <basic/Tracer.hh>


// C++ Headers
#include <string>

static basic::Tracer TR("core.chemical.sdf.MolfileIOReader.cxxtest");

using namespace core::chemical;
using namespace core::chemical::sdf;

class MolFileIOReaderTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	void test_read_rotamers() {
		TR << "Starting test " << std::endl;
		sdf::MolFileIOReader molfile_reader;

		utility::vector1< MolFileIOMoleculeOP > molfile_data = molfile_reader.parse_file( "core/chemical/sdf/multi_rotamers.sdf" );

		TR << "File parsed" << std::endl;

		utility::vector1< ResidueTypeOP > no_rot = convert_to_ResidueTypes( molfile_data, false );

		TR << "SDF file contains " << no_rot.size() << " entries " << std::endl;

		TS_ASSERT_EQUALS( no_rot.size(), 9 ); // 9 entries, each a different residue type

		utility::vector1< ResidueTypeOP > rotamers = convert_to_ResidueTypes( molfile_data, true );

		TR << "These are divided into  " << rotamers.size() << " rotamer groups " << std::endl;
		TS_ASSERT_EQUALS( rotamers.size(), 2 ); // 2 entries for different residue types

		for ( core::Size ii(1); ii<= rotamers.size(); ++ii ) {
			core::chemical::rotamers::RotamerLibrarySpecificationCOP rotspec = rotamers[ii]->rotamer_library_specification();
			TS_ASSERT( rotspec.get() );
			core::chemical::rotamers::StoredRotamerLibrarySpecificationCOP stored_rotspec(
				utility::pointer::dynamic_pointer_cast< core::chemical::rotamers::StoredRotamerLibrarySpecification const >(rotspec) );
			TS_ASSERT( stored_rotspec.get() );
			TR << "Entry " << ii << " ('" << rotamers[ii]->name() << "') contains " << stored_rotspec->coordinates().size() << " entries. " << std::endl;
			if ( ii == 1 ) {
				TS_ASSERT_EQUALS( rotamers[1]->name(), "SBT" );
				TS_ASSERT_EQUALS( stored_rotspec->coordinates().size(), 7 );
			}
			if ( ii == 2 ) {
				TS_ASSERT_EQUALS( rotamers[2]->name(), "SBF" );
				TS_ASSERT_EQUALS( stored_rotspec->coordinates().size(), 2 );
			}
		}
	}

	void test_build_rotamers() {
		TR << "Starting building test " << std::endl;

		sdf::MolFileIOReader molfile_reader;
		utility::vector1< MolFileIOMoleculeOP > molfile_data = molfile_reader.parse_file( "core/chemical/sdf/multi_rotamers.sdf" );
		utility::vector1< ResidueTypeOP > restypes = convert_to_ResidueTypes( molfile_data, true );

		core::pack::rotamers::SingleResidueRotamerLibraryCOP rotlib( core::pack::rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( *restypes[1] ) );
		TS_ASSERT( rotlib.get() );
		core::pack::rotamers::SingleLigandRotamerLibraryCOP slrotlib( utility::pointer::dynamic_pointer_cast< core::pack::rotamers::SingleLigandRotamerLibrary const >(rotlib) );
		TS_ASSERT( slrotlib.get() );

		// core::chemical::pretty_print_atomicoor( TR, *restypes[1] );

		core::pack::rotamers::RotamerVector rotamers;
		slrotlib->build_base_rotamers( *restypes[1], rotamers );

		TR << "Built " << rotamers.size() << " rotamers. " << std::endl;
		TS_ASSERT_EQUALS( rotamers.size(), 7 );

		// We're doing long int here to avoid issues with Real rounding
		// the * 3000 + 0.5 is to make sure that input precision rounds to the same value
		std::set< long int > coordinates;
		coordinates.insert( 10.8330 * 3000 + 0.5);
		coordinates.insert( 16.1490 * 3000 + 0.5);
		coordinates.insert(-11.2080 * 3000 + 0.5);
		coordinates.insert(-25.6300 * 3000 + 0.5);
		coordinates.insert(  7.3540 * 3000 + 0.5);
		coordinates.insert(  6.8840 * 3000 + 0.5);
		coordinates.insert(107.7510 * 3000 + 0.5);

		for ( core::Size ii(1); ii <= rotamers.size(); ++ii ) {
			core::Vector o_pos = rotamers[ii]->xyz("O1");
			if ( coordinates.find( o_pos.y() * 3000 + 0.5 ) == coordinates.end() ) {
				TR << "Coordinate value not found: " << o_pos.y() << std::endl;
				TS_FAIL("Cannot find coordinate in approved map.");
			} else {
				TR << "Found " << o_pos.y() << std::endl;
			}
		}

	}
};
