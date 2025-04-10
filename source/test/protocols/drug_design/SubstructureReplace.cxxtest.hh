// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/drug_design/SubstructureReplace.cxxtest.hh
/// @brief  test for SubstructureReplace Chemistry
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Project Headers
#include <core/types.hh>

#include <protocols/drug_design/SubstructureReplace.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/sdf/MolFileIOReader.hh>

// Utility Headers
#include <utility/file/FileName.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR("protocols.drug_design.SubstructureReplace.cxxtest.hh");

// --------------- Test Class --------------- //

class SubstructureReplaceTests : public CxxTest::TestSuite {

private:
	core::chemical::ResidueTypeOP restype_;

public:

	void setUp() {
		core_init();
	}

	void tearDown() {
	}

	void test_remove_hydro() {
		protocols::drug_design::SubstructureReplace replace;

		replace.substructure_database( "protocols/drug_design/benzene_stub.sdf", /*append=*/ true);
		replace.substructure_database( "protocols/drug_design/pyrimidine_stub.sdf", /*append=*/ true);
		replace.distance_threshold(1.5);

		utility::vector1< core::chemical::MutableResidueTypeOP > restypes( core::chemical::sdf::convert_to_ResidueTypes("protocols/drug_design/IClBrbenzene.sdf") );
		TS_ASSERT( restypes.size() > 0 );
		core::chemical::MutableResidueTypeOP restype( restypes[1] );

		core::chemical::NameVDMapping orig_map( *restype );

		TS_ASSERT_EQUALS( restype->natoms(), 12 );
		TS_ASSERT_EQUALS( restype->nheavyatoms(), 9 );

		//orig_map.show(TR);

		TS_ASSERT( restype->has("CL1") );
		TS_ASSERT( restype->has("I1") );
		TS_ASSERT( restype->has("BR1") );
		TS_ASSERT( restype->has("C4") );

		replace.apply(*restype);
		TS_ASSERT_EQUALS( replace.get_last_status(), core::chemical::modifications::SUCCESS );

		TS_ASSERT_EQUALS( restype->natoms(), 10 );
		TS_ASSERT_EQUALS( restype->nheavyatoms(), 9 );

		// Need to have the new nitrogens
		TS_ASSERT( restype->has("N1") );
		TS_ASSERT_EQUALS( restype->atom_type( restype->atom_vertex("N1") ).atom_type_name(), "Nhis" ); // No hydrogen
		TS_ASSERT( restype->has("N2") );
		TS_ASSERT_EQUALS( restype->atom_type( restype->atom_vertex("N2") ).atom_type_name(), "Nhis" ); // No hydrogen

		core::chemical::VDVDMapping map( replace.get_mapping() );

		//core::chemical::VDNameMapping n_map( *restype );
		//map.show(TR);
		//TR << "---------------------------------------------------------------" << std::endl;
		//n_map.show(TR);
		//TR << "---------------------------------------------------------------" << std::endl;
		//orig_map.downstream_combine( map ).downstream_combine( n_map ).show( TR );

		TS_ASSERT_DIFFERS( map[ orig_map[ " CL1" ] ], map.invalid_entry() );
		TS_ASSERT_DIFFERS( map[ orig_map[ " I1 " ] ], map.invalid_entry() );
		TS_ASSERT_DIFFERS( map[ orig_map[ " BR1" ] ], map.invalid_entry() );
		TS_ASSERT_EQUALS( map[ orig_map[ " C4 " ] ], map.invalid_entry() ); // The rings shouldn't be mapped.

		TS_ASSERT_EQUALS( restype->atom_name( map[ orig_map[ " CL1" ] ] ), " CL1" );
		TS_ASSERT_EQUALS( restype->atom_name( map[ orig_map[ " I1 " ] ] ),  " I1 " );
		TS_ASSERT_EQUALS( restype->atom_name( map[ orig_map[ " BR1" ] ] ), " BR1" );

		// Shouldn't map anything on the nitrogens.
		TS_ASSERT_EQUALS( map.reverse_lookup( restype->atom_vertex("N1") ), map.invalid_key() );
		TS_ASSERT_EQUALS( map.reverse_lookup( restype->atom_vertex("N2") ), map.invalid_key() );

	}

	void test_no_stub_match() {
		protocols::drug_design::SubstructureReplace replace;

		replace.substructure_database( "protocols/drug_design/benzene_stub.sdf", /*append=*/ true);
		replace.substructure_database( "protocols/drug_design/pyridazine_stub.sdf", /*append=*/ true);
		replace.distance_threshold(1.5);

		utility::vector1< core::chemical::MutableResidueTypeOP > restypes( core::chemical::sdf::convert_to_ResidueTypes("protocols/drug_design/IClBrbenzene.sdf") );
		TS_ASSERT( restypes.size() > 0 );
		core::chemical::MutableResidueTypeOP restype( restypes[1] );

		replace.apply(*restype);

		// Shouldn't be able to find a match.
		TS_ASSERT_EQUALS( replace.get_last_status(), core::chemical::modifications::FAIL_RETRY );
	}

	void test_add_hydro() {
		protocols::drug_design::SubstructureReplace replace;

		replace.substructure_database( "protocols/drug_design/benzene_stub.sdf", /*append=*/ true);
		replace.substructure_database( "protocols/drug_design/pyrimidine_stub.sdf", /*append=*/ true);
		replace.distance_threshold(1.5);

		utility::vector1< core::chemical::MutableResidueTypeOP > restypes( core::chemical::sdf::convert_to_ResidueTypes("protocols/drug_design/IClBrpyrimidine.sdf") );
		TS_ASSERT( restypes.size() > 0 );
		core::chemical::MutableResidueTypeOP restype( restypes[1] );

		core::chemical::NameVDMapping orig_map( *restype );

		TS_ASSERT_EQUALS( restype->natoms(), 10 );
		TS_ASSERT_EQUALS( restype->nheavyatoms(), 9 );

		TS_ASSERT( restype->has("CL1") );
		TS_ASSERT( restype->has("I1") );
		TS_ASSERT( restype->has("BR1") );
		TS_ASSERT( restype->has("C4") );
		TS_ASSERT( restype->has("H1") );

		replace.apply(*restype);
		TS_ASSERT_EQUALS( replace.get_last_status(), core::chemical::modifications::SUCCESS );

		TS_ASSERT_EQUALS( restype->natoms(), 12 );
		TS_ASSERT_EQUALS( restype->nheavyatoms(), 9 );

		// Need to have removed the nitrogens
		TS_ASSERT( ! restype->has("N1") );
		TS_ASSERT( ! restype->has("N2") );

		core::chemical::VDVDMapping map( replace.get_mapping() );

		//core::chemical::VDNameMapping n_map( *restype );
		//map.show(TR);
		//TR << "---------------------------------------------------------------" << std::endl;
		//n_map.show(TR);
		//TR << "---------------------------------------------------------------" << std::endl;
		//orig_map.downstream_combine( map ).downstream_combine( n_map ).show( TR );

		TS_ASSERT_DIFFERS( map[ orig_map[ " CL1" ] ], map.invalid_entry() );
		TS_ASSERT_DIFFERS( map[ orig_map[ " BR1" ] ], map.invalid_entry() );
		TS_ASSERT_DIFFERS( map[ orig_map[ " I1 " ] ], map.invalid_entry() );
		TS_ASSERT_DIFFERS( map[ orig_map[ " H1 " ] ], map.invalid_entry() );
		TS_ASSERT_EQUALS( map[ orig_map[ " C4 " ] ], map.invalid_entry() ); // The rings shouldn't be mapped.

		TS_ASSERT_EQUALS( restype->atom_name( map[ orig_map[ " CL1" ] ] ), " CL1" );
		TS_ASSERT_EQUALS( restype->atom_name( map[ orig_map[ " BR1" ] ] ), " BR1" );
		TS_ASSERT_EQUALS( restype->atom_name( map[ orig_map[ " I1 " ] ] ),  " I1 " );
		TS_ASSERT_EQUALS( restype->atom_name( map[ orig_map[ " H1 " ] ] ),  " H1 " );

		// Should have three hydrogens, each connected to a carbon.
		TS_ASSERT( restype->has("H1") );
		TS_ASSERT( restype->has("H2") );
		TS_ASSERT( restype->has("H3") );

		TS_ASSERT_EQUALS(restype->atom( restype->atom_base( restype->atom_vertex("H1") ) ).element(), core::chemical::element::C);
		TS_ASSERT_EQUALS(restype->atom( restype->atom_base( restype->atom_vertex("H2") ) ).element(), core::chemical::element::C);
		TS_ASSERT_EQUALS(restype->atom( restype->atom_base( restype->atom_vertex("H3") ) ).element(), core::chemical::element::C);

		// Should have six carbons, each connected to a carbon. They shouldn't map backwards
		TS_ASSERT( restype->has("C1") );
		TS_ASSERT( restype->has("C2") );
		TS_ASSERT( restype->has("C3") );
		TS_ASSERT( restype->has("C4") );
		TS_ASSERT( restype->has("C5") );
		TS_ASSERT( restype->has("C6") );
		TS_ASSERT_EQUALS( map.reverse_lookup( restype->atom_vertex("C1") ), map.invalid_key() );
		TS_ASSERT_EQUALS( map.reverse_lookup( restype->atom_vertex("C2") ), map.invalid_key() );
		TS_ASSERT_EQUALS( map.reverse_lookup( restype->atom_vertex("C3") ), map.invalid_key() );
		TS_ASSERT_EQUALS( map.reverse_lookup( restype->atom_vertex("C4") ), map.invalid_key() );
		TS_ASSERT_EQUALS( map.reverse_lookup( restype->atom_vertex("C5") ), map.invalid_key() );
		TS_ASSERT_EQUALS( map.reverse_lookup( restype->atom_vertex("C6") ), map.invalid_key() );
	}

	void test_double_bonds() {
		protocols::drug_design::SubstructureReplace replace;
		replace.distance_threshold(1.5);

		replace.substructure_database( "protocols/drug_design/ethylbenzene_stub.sdf", /*append=*/ true);
		replace.substructure_database( "protocols/drug_design/ethylpyrimidine_stub.sdf", /*append=*/ true);

		utility::vector1< core::chemical::MutableResidueTypeOP > restypes( core::chemical::sdf::convert_to_ResidueTypes("protocols/drug_design/methylclstyrene.sdf") );
		TS_ASSERT( restypes.size() > 0 );
		core::chemical::MutableResidueTypeOP restype( restypes[1] );

		core::chemical::NameVDMapping orig_map( *restype );

		TS_ASSERT_EQUALS( restype->natoms(), 19 );
		TS_ASSERT_EQUALS( restype->nheavyatoms(), 10 );

		TS_ASSERT( restype->has("CL1") );

		replace.apply(*restype);
		TS_ASSERT_EQUALS( replace.get_last_status(), core::chemical::modifications::FAIL_DO_NOT_RETRY ); // Can't find template

		replace.substructure_database( "protocols/drug_design/styrene_stub.sdf", /*append=*/ true);

		replace.apply(*restype);
		TS_ASSERT_EQUALS( replace.get_last_status(), core::chemical::modifications::FAIL_RETRY ); // Can find template, can't find replacement

		replace.substructure_database( "protocols/drug_design/vinylpyrimidine_stub.sdf", /*append=*/ true);

		replace.apply(*restype);
		TS_ASSERT_EQUALS( replace.get_last_status(), core::chemical::modifications::SUCCESS );

		TS_ASSERT_EQUALS( restype->natoms(), 17 );
		TS_ASSERT_EQUALS( restype->nheavyatoms(), 10 );

		core::chemical::VDVDMapping map( replace.get_mapping() );

		TS_ASSERT_DIFFERS( map[ orig_map[ " CL1" ] ], map.invalid_entry() );
		TS_ASSERT_EQUALS( restype->atom_name( map[ orig_map[ " CL1" ] ] ), " CL1" );

		// Need to have the new nitrogens
		TS_ASSERT( restype->has("N1") );
		TS_ASSERT_EQUALS( restype->atom_type( restype->atom_vertex("N1")).atom_type_name(), "Nhis" ); // No hydrogen
		TS_ASSERT( restype->has("N2") );
		TS_ASSERT_EQUALS( restype->atom_type( restype->atom_vertex("N2")).atom_type_name(), "Nhis" ); // No hydrogen

		// Shouldn't map anything on the nitrogens.
		TS_ASSERT_EQUALS( map.reverse_lookup( restype->atom_vertex("N1") ), map.invalid_key() );
		TS_ASSERT_EQUALS( map.reverse_lookup( restype->atom_vertex("N2") ), map.invalid_key() );
	}

	void test_HV_items() {
		protocols::drug_design::SubstructureReplace replace;

		replace.substructure_database( "protocols/drug_design/benzene.sdf", /*append=*/ true);
		replace.substructure_database( "protocols/drug_design/pyrimidine_Vstub.sdf", /*append=*/ true);
		replace.distance_threshold(1.5);

		utility::vector1< core::chemical::MutableResidueTypeOP > restypes( core::chemical::sdf::convert_to_ResidueTypes("protocols/drug_design/IClBrbenzene.sdf") );
		TS_ASSERT( restypes.size() > 0 );
		core::chemical::MutableResidueTypeOP restype( restypes[1] );

		core::chemical::NameVDMapping orig_map( *restype );

		TS_ASSERT_EQUALS( restype->natoms(), 12 );
		TS_ASSERT_EQUALS( restype->nheavyatoms(), 9 );

		//orig_map.show(TR);

		TS_ASSERT( restype->has("CL1") );
		TS_ASSERT( restype->has("I1") );
		TS_ASSERT( restype->has("BR1") );
		TS_ASSERT( restype->has("C4") );

		replace.apply(*restype);
		TS_ASSERT_EQUALS( replace.get_last_status(), core::chemical::modifications::FAIL_DO_NOT_RETRY ); // Template can't be found

		replace.H_as_dummy( true );
		replace.substructure_database( "protocols/drug_design/benzene.sdf", /*append=*/ true);

		replace.apply(*restype);
		TS_ASSERT_EQUALS( replace.get_last_status(), core::chemical::modifications::FAIL_RETRY ); // replace can't be found

		replace.V_as_dummy( true );
		replace.substructure_database( "protocols/drug_design/pyrimidine_Vstub.sdf", /*append=*/ true);

		replace.apply(*restype);
		TS_ASSERT_EQUALS( replace.get_last_status(), core::chemical::modifications::SUCCESS );

		TS_ASSERT_EQUALS( restype->natoms(), 10 );
		TS_ASSERT_EQUALS( restype->nheavyatoms(), 9 );

		// Need to have the new nitrogens
		TS_ASSERT( restype->has("N1") );
		TS_ASSERT_EQUALS( restype->atom_type( restype->atom_vertex("N1")).atom_type_name(), "Nhis" ); // No hydrogen
		TS_ASSERT( restype->has("N2") );
		TS_ASSERT_EQUALS( restype->atom_type( restype->atom_vertex("N2")).atom_type_name(), "Nhis" ); // No hydrogen

		core::chemical::VDVDMapping map( replace.get_mapping() );

		//core::chemical::VDNameMapping n_map( *restype );
		//map.show(TR);
		//TR << "---------------------------------------------------------------" << std::endl;
		//n_map.show(TR);
		//TR << "---------------------------------------------------------------" << std::endl;
		//orig_map.downstream_combine( map ).downstream_combine( n_map ).show( TR );

		TS_ASSERT_DIFFERS( map[ orig_map[ " CL1" ] ], map.invalid_entry() );
		TS_ASSERT_DIFFERS( map[ orig_map[ " I1 " ] ], map.invalid_entry() );
		TS_ASSERT_DIFFERS( map[ orig_map[ " BR1" ] ], map.invalid_entry() );
		TS_ASSERT_EQUALS( map[ orig_map[ " C4 " ] ], map.invalid_entry() ); // The rings shouldn't be mapped.

		TS_ASSERT_EQUALS( restype->atom_name( map[ orig_map[ " CL1" ] ] ), " CL1" );
		TS_ASSERT_EQUALS( restype->atom_name( map[ orig_map[ " I1 " ] ] ),  " I1 " );
		TS_ASSERT_EQUALS( restype->atom_name( map[ orig_map[ " BR1" ] ] ), " BR1" );

		// Shouldn't map anything on the nitrogens.
		TS_ASSERT_EQUALS( map.reverse_lookup( restype->atom_vertex("N1") ), map.invalid_key() );
		TS_ASSERT_EQUALS( map.reverse_lookup( restype->atom_vertex("N2") ), map.invalid_key() );
	}

	// @brief Test if a structure which has extra hydrogens can gracefully delete them.
	void test_dropped_hydrogens() {
		protocols::drug_design::SubstructureReplace replace;

		replace.H_as_dummy( true );
		replace.substructure_database( "protocols/drug_design/benzene.sdf", /*append=*/ true);
		replace.substructure_database( "protocols/drug_design/piperazine_stub.sdf", /*append=*/ true);
		replace.distance_threshold(1.5);

		utility::vector1< core::chemical::MutableResidueTypeOP > restypes( core::chemical::sdf::convert_to_ResidueTypes("protocols/drug_design/piperazine_prot.sdf") );
		TS_ASSERT( restypes.size() > 0 );
		core::chemical::MutableResidueTypeOP restype( restypes[1] );

		TS_ASSERT_EQUALS( restype->natoms(), 18 );
		TS_ASSERT_EQUALS( restype->nheavyatoms(), 6 );

		replace.apply(*restype);
		TS_ASSERT_EQUALS( replace.get_last_status(), core::chemical::modifications::SUCCESS );

		TS_ASSERT_EQUALS( restype->natoms(), 12 );
		TS_ASSERT_EQUALS( restype->nheavyatoms(), 6 );

		// But it shouldn't work if the items are not hydrogens
		utility::vector1< core::chemical::MutableResidueTypeOP > restypes2( core::chemical::sdf::convert_to_ResidueTypes("protocols/drug_design/piperazine_quat.sdf") );
		TS_ASSERT( restypes2.size() > 0 );
		core::chemical::MutableResidueTypeOP restype2( restypes2[1] );

		TS_ASSERT_EQUALS( restype2->natoms(), 30 );
		TS_ASSERT_EQUALS( restype2->nheavyatoms(), 10 );

		replace.apply(*restype2);
		TS_ASSERT_EQUALS( replace.get_last_status(), core::chemical::modifications::FAIL_DO_NOT_RETRY );
	}

};
