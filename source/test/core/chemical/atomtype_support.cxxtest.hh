// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/chemical/atomtype_support.cxxtest.hh
/// @brief unit tests for ResidueType
/// @author Rocco Moretti (rmorettiase@gmail.com)


// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/atomtype_support.hh>

// Project Headers
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/MMAtomType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/PoseResidueTypeSet.hh>
#include <core/chemical/residue_io.hh>

// Platform Headers
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <string>
#include <ostream>


using std::endl;
using std::string;
using basic::Tracer;
using core::chemical::AtomTypeSetCAP;
using core::chemical::AtomType;
using core::chemical::ChemicalManager;
using core::chemical::ElementSetCAP;
using core::chemical::MMAtomTypeSetCAP;
using core::chemical::orbitals::OrbitalTypeSetCAP;
using core::chemical::ResidueType;
using core::chemical::ResidueTypeSet;
using core::chemical::FA_STANDARD;
using utility::vector1;

static Tracer TR("core.chemical.atomtype_support.cxxtest");

class atomtype_support_Tests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	void test_fullatom_retyping() {
		using namespace core::chemical;

		ChemicalManager * cm(ChemicalManager::get_instance());
		ResidueTypeSetCOP rsd_types( cm->residue_type_set(FA_STANDARD) );

		utility::io::izstream paramslist("core/chemical/params/retype_list.txt");
		std::string filename;
		paramslist >> filename;
		while ( paramslist ) {
			TR << "Retyping " << filename << std::endl;
			core::chemical::ResidueTypeOP rsd = read_topology_file("core/chemical/"+filename, rsd_types);
			core::chemical::ResidueTypeOP ref( new core::chemical::ResidueType(*rsd) );

			for ( core::Size ii(1); ii <= rsd->natoms(); ++ii ) {
				rsd->atom(ii).atom_type_index(0);
			}
			core::chemical::rosetta_retype_fullatom(*rsd);

			for ( core::Size ii(1); ii <= rsd->natoms(); ++ii ) {
				TS_ASSERT( rsd->atom(ii).atom_type_index() != 0 );
				TS_ASSERT_EQUALS( rsd->atom_type(ii).name(), ref->atom_type(ii).name() );
			}

			paramslist >> filename;
		}

		// TODO: Test partial reassignment.
	}
};
