// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rotamers/SingleLigandRotamerLibrary.cxxtest.hh
/// @brief
/// @author Steven Lewis smlewi@gmail.com

// Test headers
#include <test/core/pack/rotamers/SingleLigandRotamerLibrary_and_StoredRotamerLibrarySpecification.util.hh>
#include <core/pack/rotamers/SingleLigandRotamerLibrary.hh>

#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/residue_io.hh>

#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>


#include <core/types.hh>

#include <basic/Tracer.hh>

#include <sstream>

static basic::Tracer TR("core.pack.rotamers.SingleLigandRotamerLibrary.cxxtest");

//typedef std::map< std::string, core::Vector > NamePosMap;
//typedef numeric::xyzVector< core::Real > Vector; //probable
using namespace core::pack::rotamers; //for NamePosMap
typedef utility::vector1< NamePosMap > AtomPositions;

class SingleLigandRotamerLibraryTests : public CxxTest::TestSuite {

private:
	AtomPositions reference_;
	core::Real const ref_energy_ = 77.0;

public:

	// Shared initialization goes here.
	void setUp() {
		core_init();

		//fill container with 'right' data
		generate_comparator( reference_ );

	}

	// Shared finalization goes here.
	void tearDown() {
	}


	//here we test that init_from_file produces the same vector1<NamePosMap> we expect
	//we call lib_from_file both directly, and implicitly via the Factory
	//(sure, could be two functions)
	void test_init_from_file_explicit_and_implicit() {

		//First make ResidueType
		using namespace core::chemical;
		ResidueTypeSetCOP rtsCOP(ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD));

		core::chemical::ResidueTypeOP RTOP(core::chemical::read_topology_file( "core/pack/rotamers/4tim.params", rtsCOP ));


		//Run SingleLigandRotamerLibrary machinery to read same data from file

		SingleLigandRotamerLibrary lib_from_file;
		lib_from_file.init_from_file("core/pack/rotamers/4tim_confs.pdb", *RTOP);

		AtomPositions const & data(lib_from_file.get_atom_positions());

		//I assume this is zero since it's effectively unset; ctor should take care
		//of it.  Exact comparison should be ok since 0 is exactly represented, but
		//feel free to make this delta if needed.
		TS_ASSERT_EQUALS(lib_from_file.get_reference_energy(), ref_energy_);

		//This DOES compile and work, but I don't really trust it; the coordinates should be delta'd
		//TS_ASSERT_EQUALS(reference_, data);
		compare_AtomPositions( reference_, data ); //do the TS_ASSERTS manually instead

		//Also use the implicit call to init_from_file from the Factory
		//This is using the actual PDB_ROTAMERS string
		//unfortunately the Factory returns a parent class pointer, so we have to downcast
		SingleResidueRotamerLibraryCOP SRRL(SingleResidueRotamerLibraryFactory::get_instance()->get( *RTOP ));
		TS_ASSERT(SRRL);

		SingleLigandRotamerLibraryCOP lib_from_factory( utility::pointer::dynamic_pointer_cast< SingleLigandRotamerLibrary const >(SRRL) );
		TS_ASSERT(lib_from_factory);
		assert(lib_from_factory);

		TS_ASSERT_EQUALS(lib_from_factory->get_reference_energy(), ref_energy_);
		compare_AtomPositions( reference_, lib_from_factory->get_atom_positions() );

	}

	//Now we test that the init_from_stream utility function returns the correct object
	void test_init_from_stream() {

		//create ResidueType to test with
		//not just using the class's version because it had PDB_ROTAMERS set
		//(which we can clear, but let's be safe)
		std::istringstream param_stream;
		param_stream.str(params_string());

		core::chemical::ResidueTypeSetCOP rtsCOP(core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FULL_ATOM_t ));

		core::chemical::ResidueTypeOP RTOP(core::chemical::read_topology_file( param_stream, "dummy_filename", rtsCOP ));

		//create stringstream to init from
		std::istringstream conformers_stream;
		conformers_stream.str(conformers_pdb_string());

		//empty data to fill
		AtomPositions data;
		core::Real e_ref(0.0);

		//test utility function - refactoring this function from SLRL::init_from_file was the point of this exercise
		core::pack::rotamers::rotamer_information_from_PDB_stream(conformers_stream, *RTOP, data, e_ref);

		TS_ASSERT_EQUALS(ref_energy_, e_ref);
		compare_AtomPositions( reference_, data );

	}

};
