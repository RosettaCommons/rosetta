// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rotamers/StoredRotamerLibrarySpecification.cxxtest.hh
/// @brief
/// @author Steven Lewis smlewi@gmail.com

// Test headers
#include <test/core/pack/rotamers/SingleLigandRotamerLibrary_and_StoredRotamerLibrarySpecification.util.hh>
#include <core/chemical/rotamers/StoredRotamerLibrarySpecification.hh>

#include <core/pack/rotamers/StoredRotamerLibraryCreator.hh>

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

static basic::Tracer TR("core.pack.rotamers.StoredRotamerLibrarySpecification.cxxtest");

//typedef std::map< std::string, core::Vector > NamePosMap;
//typedef numeric::xyzVector< core::Real > Vector; //probable
using namespace core::pack::rotamers; //for NamePosMap
typedef utility::vector1< NamePosMap > AtomPositions;

using namespace core::chemical::rotamers;

class StoredRotamerLibrarySpecificationTests : public CxxTest::TestSuite {

private:
	AtomPositions reference_;
	AtomPositions conformers_;
	core::Real const ref_energy_ = 77.0;
	core::Real e_ref_;
	core::chemical::ResidueTypeOP RTOP_;

public:

	// Shared initialization goes here.
	void setUp() {
		core_init();

		//fill container with 'right' data
		generate_comparator( reference_ );

		//create ResidueType to test with
		//this is same setup as SingleLigandRotamerLibraryTests::test_init_from_stream
		std::istringstream param_stream;
		param_stream.str(params_string());

		core::chemical::ResidueTypeSetCOP rtsCOP(core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FULL_ATOM_t ));

		RTOP_ = core::chemical::read_topology_file( param_stream, "dummy_filename", rtsCOP );

		//create stringstream to init from
		std::istringstream conformers_stream;
		conformers_stream.str(conformers_pdb_string());

		//fill conformers
		e_ref_ = 0.0;
		core::pack::rotamers::rotamer_information_from_PDB_stream(conformers_stream, *RTOP_, conformers_, e_ref_);

	}

	// Shared finalization goes here.
	void tearDown() {
	}


	void compare_SLRL_to_reference(SingleLigandRotamerLibraryCOP SLRL) {

		TS_ASSERT(SLRL);
		assert(SLRL);

		TS_ASSERT_EQUALS(SLRL->get_reference_energy(), ref_energy_);
		compare_AtomPositions( reference_, SLRL->get_atom_positions() );

	}

	//utility function shared between several tests
	void get_and_compare_SLRL_to_reference(StoredRotamerLibrarySpecificationOP SRLS) {

		//load SRLS into RTOP
		RTOP_->strip_rotamer_library_specification();
		RTOP_->rotamer_library_specification(SRLS);

		//get SingleLigandRotamerLibrary from StoredRotamerLibraryCreator
		StoredRotamerLibraryCreator SRLC;
		SingleResidueRotamerLibraryCOP SRRL(SRLC.create(*RTOP_));
		SingleLigandRotamerLibraryCOP SLRL( utility::pointer::dynamic_pointer_cast< SingleLigandRotamerLibrary const >(SRRL) );

		compare_SLRL_to_reference(SLRL);

		//get SingleLigandRotamerLibrary from SingleResidueRotamerLibraryFactory instead to make sure that works
		SRRL = SingleResidueRotamerLibraryFactory::get_instance()->get( *RTOP_, false );
		SLRL = utility::pointer::dynamic_pointer_cast< SingleLigandRotamerLibrary const >(SRRL);
		compare_SLRL_to_reference(SLRL);


	}

	void test_add_rotamer() {

		//the setUp function has prepared the conformers_ data for us

		//create our StoredRotamerLibrarySpecification, and test its add_rotamer function
		StoredRotamerLibrarySpecificationOP SRLS_onebyone(new StoredRotamerLibrarySpecification);
		for ( auto const & item : conformers_ ) SRLS_onebyone->add_rotamer(item);

		SRLS_onebyone->set_reference_energy(e_ref_);

		//extracts a SingleLigandRotamerLibrary and compares it to the reference data
		get_and_compare_SLRL_to_reference(SRLS_onebyone);

	}

	void test_add_conformers() {

		//the setUp function has prepared the conformers_ data for us

		//create our StoredRotamerLibrarySpecification, and test its add_rotamer function
		StoredRotamerLibrarySpecificationOP SRLS_all(new StoredRotamerLibrarySpecification);
		SRLS_all->add_rotamers(conformers_);

		SRLS_all->set_reference_energy(e_ref_);

		//extracts a SingleLigandRotamerLibrary and compares it to the reference data
		get_and_compare_SLRL_to_reference(SRLS_all);

	}

};
