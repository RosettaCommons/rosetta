// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/cyclic_peptide/PeptideStubMoverTests.cxxtest.hh
/// @brief  Unit tests for the PeptideStubMover, used to build peptide geometry.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/util/rosettascripts.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/cyclic_peptide/PeptideStubMover.hh>
#include <protocols/simple_moves/ModifyVariantTypeMover.hh>

// Protocols Headers
#include <protocols/simple_moves/MutateResidue.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/variant_util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

// Numeric Headers
#include <numeric/conversions.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>
#include <utility/pointer/memory.hh>

static basic::Tracer TR("PeptideStubMoverTests");

#define TEST_DELTA 0.00001

class PeptideStubMoverTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options( "-write_all_connect_info true" );
		core::pose::PoseOP pose( utility::pointer::make_shared< core::pose::Pose >() );
		core::import_pose::pose_from_file( *pose, "protocols/cyclic_peptide/peptidestubmover_test_pose.pdb", false, core::import_pose::PDB_file );
		for ( core::Size i(1), imax(pose->total_residue()); i<=imax; ++i ) {
			protocols::simple_moves::MutateResidue mutres( i, "VAL" ); //Mutate the whole pose to valine.
			mutres.apply(*pose);
		}
		core::pose::add_lower_terminus_type_to_pose_residue(*pose, 1);
		core::pose::add_upper_terminus_type_to_pose_residue(*pose, pose->size());

		pose->update_residue_neighbors();
		pose_ = pose; //Nonconst to const.
	}

	void tearDown(){

	}

	void test_append() {
		core::pose::Pose pose( *(pose_) ); //Local copy

		core::Size const original_pose_size( pose.total_residue() );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Append", "ALA", 0, false, "", 0, pose.total_residue(), "" );
		stubmover.apply(pose);
		TS_ASSERT_EQUALS(pose.total_residue(), original_pose_size+1 );
		TS_ASSERT_EQUALS( pose.residue_type( pose.total_residue() ).name3(), "ALA" );
		TS_ASSERT_EQUALS( pose.residue_type( original_pose_size ).name3(), "VAL" );
		TS_ASSERT_EQUALS( pose.residue(pose.total_residue()).connected_residue_at_lower(), original_pose_size );
		TS_ASSERT_EQUALS( pose.residue(original_pose_size).connected_residue_at_upper(), pose.total_residue() );
		TS_ASSERT_DELTA( pose.residue( original_pose_size ).xyz("C").distance( pose.residue( original_pose_size + 1 ).xyz("N") ), 1.328685, 0.01 ); //Check peptide bond length.
	}

	void test_append2() {
		core::pose::Pose pose( *(pose_) ); //Local copy

		core::Size const original_pose_size( pose.total_residue() );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Append", "ALA", 0, false, "", 0, 0, "" ); //Should automatically select last
		stubmover.apply(pose);
		TS_ASSERT_EQUALS(pose.total_residue(), original_pose_size+1 );
		TS_ASSERT_EQUALS( pose.residue_type( pose.total_residue() ).name3(), "ALA" );
		TS_ASSERT_EQUALS( pose.residue_type( original_pose_size ).name3(), "VAL" );
		TS_ASSERT_EQUALS( pose.residue(pose.total_residue()).connected_residue_at_lower(), original_pose_size );
		TS_ASSERT_EQUALS( pose.residue(original_pose_size).connected_residue_at_upper(), pose.total_residue() );
		TS_ASSERT_DELTA( pose.residue( original_pose_size ).xyz("C").distance( pose.residue( original_pose_size + 1 ).xyz("N") ), 1.328685, 0.01 ); //Check peptide bond length.
	}

	/// @brief Test appending to an N-methylamidated C-terminus and prepending on an
	/// N-acetylated N-terminus.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
	void test_append_NME_and_prepend_ACE() {
		using namespace protocols::cyclic_peptide;

		utility::vector1< core::Real > const phivals{ 21.395, -73.554, 129.445 };
		utility::vector1< core::Real > const psivals{ 166.44, -144.31, 60.32 };
		utility::vector1< core::Real > const omgvals{ 178.33, -174.41, 179.12 };
		core::Real const omg_pre_val( -177.773 );

		core::pose::Pose pose;

		// Set up the peptide (three-residue pose with acetylated N-terminus and N-methylamidated C-terminus):
		{
			PeptideStubMover mover1;
			mover1.set_reset_mode( true );
			mover1.add_residue( "Append", "ALA", 0, true, "", 0, 0, "" );
			mover1.add_residue( "Append", "CYS", 1, false, "", 0, 0, "" );
			mover1.add_residue( "Append", "SER", 2, false, "", 0, 0, "" );
			mover1.apply( pose );

			core::select::residue_selector::ResidueIndexSelectorCOP select_one(
				utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >( 1 )
			);
			core::select::residue_selector::ResidueIndexSelectorCOP select_three(
				utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >( 3 )
			);

			protocols::simple_moves::ModifyVariantTypeMover add_ace;
			add_ace.set_residue_selector( select_one );
			add_ace.set_additional_type_to_add("N_ACETYLATION");
			add_ace.set_update_polymer_bond_dependent_atoms(true);
			add_ace.apply(pose);

			protocols::simple_moves::ModifyVariantTypeMover add_nme;
			add_nme.set_residue_selector( select_three );
			add_nme.set_additional_type_to_add("C_METHYLAMIDATION");
			add_nme.set_update_polymer_bond_dependent_atoms(true);
			add_nme.apply(pose);

			for ( core::Size i(1); i<=3; ++i ) {
				pose.set_phi( i, phivals[i] );
				pose.set_psi( i, psivals[i] );
				pose.set_omega( i, omgvals[i] );
			}
			core::chemical::ResidueType const & restype( pose.residue_type(1) );
			pose.conformation().set_torsion_angle(
				core::id::AtomID( restype.atom_index( "CQ" ), 1 ),
				core::id::AtomID( restype.atom_index( "CP" ), 1 ),
				core::id::AtomID( restype.atom_index( "N" ), 1 ),
				core::id::AtomID( restype.atom_index( "CA" ), 1 ),
				numeric::conversions::radians( omg_pre_val )
			);
			pose.update_residue_neighbors();
		}

		// Check the conformation:
		numeric::xyzVector< core::Real > old_ca_position( pose.xyz( core::id::NamedAtomID( "CA", 1 ) ) );
		numeric::xyzVector< core::Real > old_cb_position( pose.xyz( core::id::NamedAtomID( "CB", 1 ) ) );
		{
			TS_ASSERT_EQUALS( pose.total_residue(), 3 );

			core::chemical::ResidueType const & restype1( pose.residue_type(1) );
			core::chemical::ResidueType const & restype2( pose.residue_type(2) );
			core::chemical::ResidueType const & restype3( pose.residue_type(3) );

			TS_ASSERT_EQUALS( restype1.name(), "ALA:N_acetylated" );
			TS_ASSERT_EQUALS( restype2.name(), "CYS" );
			TS_ASSERT_EQUALS( restype3.name(), "SER:C_methylamidated" );

			TS_ASSERT_DELTA( numeric::conversions::degrees( pose.conformation().torsion_angle(
				core::id::AtomID( restype1.atom_index( "CQ" ), 1 ),
				core::id::AtomID( restype1.atom_index( "CP" ), 1 ),
				core::id::AtomID( restype1.atom_index( "N" ), 1 ),
				core::id::AtomID( restype1.atom_index( "CA" ), 1 )
				) ),
				omg_pre_val, TEST_DELTA
			);

			TS_ASSERT_DELTA( numeric::conversions::degrees( pose.conformation().torsion_angle(
				core::id::AtomID( restype1.atom_index( "CP" ), 1 ),
				core::id::AtomID( restype1.atom_index( "N" ), 1 ),
				core::id::AtomID( restype1.atom_index( "CA" ), 1 ),
				core::id::AtomID( restype1.atom_index( "C" ), 1 )
				) ),
				phivals[1], TEST_DELTA
			);

			TS_ASSERT_DELTA( numeric::conversions::degrees( pose.conformation().torsion_angle(
				core::id::AtomID( restype1.atom_index( "N" ), 1 ),
				core::id::AtomID( restype1.atom_index( "CA" ), 1 ),
				core::id::AtomID( restype1.atom_index( "C" ), 1 ),
				core::id::AtomID( restype2.atom_index( "N" ), 2 )
				) ),
				psivals[1], TEST_DELTA
			);

			TS_ASSERT_DELTA( numeric::conversions::degrees( pose.conformation().torsion_angle(
				core::id::AtomID( restype1.atom_index( "CA" ), 1 ),
				core::id::AtomID( restype1.atom_index( "C" ), 1 ),
				core::id::AtomID( restype2.atom_index( "N" ), 2 ),
				core::id::AtomID( restype2.atom_index( "CA" ), 2 )
				) ),
				omgvals[1], TEST_DELTA
			);

			TS_ASSERT_DELTA( numeric::conversions::degrees( pose.conformation().torsion_angle(
				core::id::AtomID( restype1.atom_index( "C" ), 1 ),
				core::id::AtomID( restype2.atom_index( "N" ), 2 ),
				core::id::AtomID( restype2.atom_index( "CA" ), 2 ),
				core::id::AtomID( restype2.atom_index( "C" ), 2)
				) ),
				phivals[2], TEST_DELTA
			);

			TS_ASSERT_DELTA( numeric::conversions::degrees( pose.conformation().torsion_angle(
				core::id::AtomID( restype2.atom_index( "N" ), 2 ),
				core::id::AtomID( restype2.atom_index( "CA" ), 2 ),
				core::id::AtomID( restype2.atom_index( "C" ), 2 ),
				core::id::AtomID( restype3.atom_index( "N" ), 3 )
				) ),
				psivals[2], TEST_DELTA
			);

			TS_ASSERT_DELTA( numeric::conversions::degrees( pose.conformation().torsion_angle(
				core::id::AtomID( restype2.atom_index( "CA" ), 2 ),
				core::id::AtomID( restype2.atom_index( "C" ), 2 ),
				core::id::AtomID( restype3.atom_index( "N" ), 3 ),
				core::id::AtomID( restype3.atom_index( "CA" ), 3 )
				) ),
				omgvals[2], TEST_DELTA
			);

			TS_ASSERT_DELTA( numeric::conversions::degrees( pose.conformation().torsion_angle(
				core::id::AtomID( restype2.atom_index( "C" ), 2 ),
				core::id::AtomID( restype3.atom_index( "N" ), 3 ),
				core::id::AtomID( restype3.atom_index( "CA" ), 3 ),
				core::id::AtomID( restype3.atom_index( "C" ), 3)
				) ),
				phivals[3], TEST_DELTA
			);

			TS_ASSERT_DELTA( numeric::conversions::degrees( pose.conformation().torsion_angle(
				core::id::AtomID( restype3.atom_index( "N" ), 3 ),
				core::id::AtomID( restype3.atom_index( "CA" ), 3 ),
				core::id::AtomID( restype3.atom_index( "C" ), 3 ),
				core::id::AtomID( restype3.atom_index( "NR" ), 3 )
				) ),
				psivals[3], TEST_DELTA
			);

			TS_ASSERT_DELTA( numeric::conversions::degrees( pose.conformation().torsion_angle(
				core::id::AtomID( restype3.atom_index( "CA" ), 3 ),
				core::id::AtomID( restype3.atom_index( "C" ), 3 ),
				core::id::AtomID( restype3.atom_index( "NR" ), 3 ),
				core::id::AtomID( restype3.atom_index( "CS" ), 3 )
				) ),
				omgvals[3], TEST_DELTA
			);
		}

		// Prepend and append residues:
		{
			//pose.dump_pdb( "ACE_NME_TEST_BEFORE.pdb" );
			PeptideStubMover mover2;
			mover2.set_reset_mode( false );
			mover2.add_residue( "Append", "VAL", 3, false, "", 0, 0, "" );
			mover2.add_residue( "Prepend", "GLY", 1, false, "", 0, 0, "" );
			mover2.apply( pose );
			//pose.dump_pdb( "ACE_NME_TEST_AFTER.pdb" );
		}

		// Check the new conformation:
		numeric::xyzVector< core::Real > new_ca_position( pose.xyz( core::id::NamedAtomID( "CA", 2 ) ) );
		numeric::xyzVector< core::Real > new_cb_position( pose.xyz( core::id::NamedAtomID( "CB", 2 ) ) );
		{
			TS_ASSERT_EQUALS( pose.total_residue(), 5 );

			TS_ASSERT_DELTA( new_ca_position.x(), old_ca_position.x(), TEST_DELTA );
			TS_ASSERT_DELTA( new_ca_position.y(), old_ca_position.y(), TEST_DELTA );
			TS_ASSERT_DELTA( new_ca_position.z(), old_ca_position.z(), TEST_DELTA );
			TS_ASSERT_DELTA( new_cb_position.x(), old_cb_position.x(), TEST_DELTA );
			TS_ASSERT_DELTA( new_cb_position.y(), old_cb_position.y(), TEST_DELTA );
			TS_ASSERT_DELTA( new_cb_position.z(), old_cb_position.z(), TEST_DELTA );

			core::chemical::ResidueType const & restype1( pose.residue_type(1) );
			core::chemical::ResidueType const & restype2( pose.residue_type(2) );
			core::chemical::ResidueType const & restype3( pose.residue_type(3) );
			core::chemical::ResidueType const & restype4( pose.residue_type(4) );
			core::chemical::ResidueType const & restype5( pose.residue_type(5) );

			TS_ASSERT_EQUALS( restype1.name(), "GLY" );
			TS_ASSERT_EQUALS( restype2.name(), "ALA" );
			TS_ASSERT_EQUALS( restype3.name(), "CYS" );
			TS_ASSERT_EQUALS( restype4.name(), "SER" );
			TS_ASSERT_EQUALS( restype5.name(), "VAL" );

			TS_ASSERT_DELTA( numeric::conversions::degrees( pose.conformation().torsion_angle(
				core::id::AtomID( restype1.atom_index( "CA" ), 1 ),
				core::id::AtomID( restype1.atom_index( "C" ), 1 ),
				core::id::AtomID( restype2.atom_index( "N" ), 2 ),
				core::id::AtomID( restype2.atom_index( "CA" ), 2 )
				) ),
				omg_pre_val, TEST_DELTA
			);

			TS_ASSERT_DELTA( numeric::conversions::degrees( pose.conformation().torsion_angle(
				core::id::AtomID( restype1.atom_index( "C" ), 1 ),
				core::id::AtomID( restype2.atom_index( "N" ), 2 ),
				core::id::AtomID( restype2.atom_index( "CA" ), 2 ),
				core::id::AtomID( restype2.atom_index( "C" ), 2 )
				) ),
				phivals[1], TEST_DELTA
			);

			TS_ASSERT_DELTA( numeric::conversions::degrees( pose.conformation().torsion_angle(
				core::id::AtomID( restype2.atom_index( "N" ), 2 ),
				core::id::AtomID( restype2.atom_index( "CA" ), 2 ),
				core::id::AtomID( restype2.atom_index( "C" ), 2 ),
				core::id::AtomID( restype3.atom_index( "N" ), 3 )
				) ),
				psivals[1], TEST_DELTA
			);

			TS_ASSERT_DELTA( numeric::conversions::degrees( pose.conformation().torsion_angle(
				core::id::AtomID( restype2.atom_index( "CA" ), 2 ),
				core::id::AtomID( restype2.atom_index( "C" ), 2 ),
				core::id::AtomID( restype3.atom_index( "N" ), 3 ),
				core::id::AtomID( restype3.atom_index( "CA" ), 3 )
				) ),
				omgvals[1], TEST_DELTA
			);

			TS_ASSERT_DELTA( numeric::conversions::degrees( pose.conformation().torsion_angle(
				core::id::AtomID( restype2.atom_index( "C" ), 2 ),
				core::id::AtomID( restype3.atom_index( "N" ), 3 ),
				core::id::AtomID( restype3.atom_index( "CA" ), 3 ),
				core::id::AtomID( restype3.atom_index( "C" ), 3)
				) ),
				phivals[2], TEST_DELTA
			);

			TS_ASSERT_DELTA( numeric::conversions::degrees( pose.conformation().torsion_angle(
				core::id::AtomID( restype3.atom_index( "N" ), 3 ),
				core::id::AtomID( restype3.atom_index( "CA" ), 3 ),
				core::id::AtomID( restype3.atom_index( "C" ), 3 ),
				core::id::AtomID( restype4.atom_index( "N" ), 4 )
				) ),
				psivals[2], TEST_DELTA
			);

			TS_ASSERT_DELTA( numeric::conversions::degrees( pose.conformation().torsion_angle(
				core::id::AtomID( restype3.atom_index( "CA" ), 3 ),
				core::id::AtomID( restype3.atom_index( "C" ), 3 ),
				core::id::AtomID( restype4.atom_index( "N" ), 4 ),
				core::id::AtomID( restype4.atom_index( "CA" ), 4 )
				) ),
				omgvals[2], TEST_DELTA
			);

			TS_ASSERT_DELTA( numeric::conversions::degrees( pose.conformation().torsion_angle(
				core::id::AtomID( restype3.atom_index( "C" ), 3 ),
				core::id::AtomID( restype4.atom_index( "N" ), 4 ),
				core::id::AtomID( restype4.atom_index( "CA" ), 4 ),
				core::id::AtomID( restype4.atom_index( "C" ), 4)
				) ),
				phivals[3], TEST_DELTA
			);

			TS_ASSERT_DELTA( numeric::conversions::degrees( pose.conformation().torsion_angle(
				core::id::AtomID( restype4.atom_index( "N" ), 4 ),
				core::id::AtomID( restype4.atom_index( "CA" ), 4 ),
				core::id::AtomID( restype4.atom_index( "C" ), 4 ),
				core::id::AtomID( restype5.atom_index( "N" ), 5 )
				) ),
				psivals[3], TEST_DELTA
			);

			TS_ASSERT_DELTA( numeric::conversions::degrees( pose.conformation().torsion_angle(
				core::id::AtomID( restype4.atom_index( "CA" ), 4 ),
				core::id::AtomID( restype4.atom_index( "C" ), 4 ),
				core::id::AtomID( restype5.atom_index( "N" ), 5 ),
				core::id::AtomID( restype5.atom_index( "CA" ), 5 )
				) ),
				omgvals[3], TEST_DELTA
			);
		}
	}

	void test_multi_append() {
		core::pose::Pose pose( *(pose_) ); //Local copy

		core::Size const original_pose_size( pose.total_residue() );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		for ( core::Size i(0); i<5; ++i ) {
			stubmover.add_residue( "Append", "ALA", 0, false, "", 0, pose.total_residue() + i, "" );
		}
		stubmover.apply(pose);
		TS_ASSERT_EQUALS(pose.total_residue(), original_pose_size+5 );
		TS_ASSERT_EQUALS( pose.residue_type( original_pose_size ).name3(), "VAL" );
		for ( core::Size i(1); i<=5; ++i ) {
			TS_ASSERT_EQUALS( pose.residue_type( original_pose_size + i ).name3(), "ALA" );
			TS_ASSERT_EQUALS( pose.residue( original_pose_size + i ).connected_residue_at_lower(), original_pose_size + i - 1 );
			TS_ASSERT_EQUALS( pose.residue( original_pose_size + i - 1 ).connected_residue_at_upper(), original_pose_size + i );
			TS_ASSERT_DELTA( pose.residue( original_pose_size + i - 1 ).xyz("C").distance( pose.residue( original_pose_size + i ).xyz("N") ), 1.328685, 0.01 ); //Check peptide bond length.
		}
	}

	void test_append_repeat() {
		core::pose::Pose pose( *(pose_) ); //Local copy

		core::Size const original_pose_size( pose.total_residue() );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Append", "ALA", 0, false, "", 5, pose.total_residue(), "" );
		stubmover.apply(pose);
		TS_ASSERT_EQUALS(pose.total_residue(), original_pose_size+5 );
		TS_ASSERT_EQUALS( pose.residue_type( original_pose_size ).name3(), "VAL" );
		for ( core::Size i(1); i<=5; ++i ) {
			TS_ASSERT_EQUALS( pose.residue_type( original_pose_size + i ).name3(), "ALA" );
			TS_ASSERT_EQUALS( pose.residue( original_pose_size + i ).connected_residue_at_lower(), original_pose_size + i - 1 );
			TS_ASSERT_EQUALS( pose.residue( original_pose_size + i - 1 ).connected_residue_at_upper(), original_pose_size + i );
			TS_ASSERT_DELTA( pose.residue( original_pose_size + i - 1 ).xyz("C").distance( pose.residue( original_pose_size + i ).xyz("N") ), 1.328685, 0.01 ); //Check peptide bond length.
		}
	}

	void test_prepend() {
		core::pose::Pose pose( *(pose_) ); //Local copy

		core::Size const original_pose_size( pose.total_residue() );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Prepend", "ALA", 0, false, "", 0, 1, "" );
		stubmover.apply(pose);
		TS_ASSERT_EQUALS(pose.total_residue(), original_pose_size+1 );
		TS_ASSERT_EQUALS( pose.residue_type( pose.total_residue() ).name3(), "VAL" );
		TS_ASSERT_EQUALS( pose.residue_type( original_pose_size ).name3(), "VAL" );
		TS_ASSERT_EQUALS( pose.residue_type( 2 ).name3(), "VAL" );
		TS_ASSERT_EQUALS( pose.residue_type( 1 ).name3(), "ALA" );
		TS_ASSERT_EQUALS( pose.residue(2).connected_residue_at_lower(), 1 );
		TS_ASSERT_EQUALS( pose.residue(1).connected_residue_at_upper(), 2 );
		TS_ASSERT_DELTA( pose.residue(1).xyz("C").distance( pose.residue(2).xyz("N") ), 1.328685, 0.01 ); //Check peptide bond length.
	}

	void test_prepend2() {
		core::pose::Pose pose( *(pose_) ); //Local copy

		core::Size const original_pose_size( pose.total_residue() );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Prepend", "ALA", 0, false, "", 0, 0, "" ); //Should auto-detect?
		stubmover.apply(pose);
		TS_ASSERT_EQUALS(pose.total_residue(), original_pose_size+1 );
		TS_ASSERT_EQUALS( pose.residue_type( pose.total_residue() ).name3(), "VAL" );
		TS_ASSERT_EQUALS( pose.residue_type( original_pose_size ).name3(), "VAL" );
		TS_ASSERT_EQUALS( pose.residue_type( 2 ).name3(), "VAL" );
		TS_ASSERT_EQUALS( pose.residue_type( 1 ).name3(), "ALA" );
		TS_ASSERT_EQUALS( pose.residue(2).connected_residue_at_lower(), 1 );
		TS_ASSERT_EQUALS( pose.residue(1).connected_residue_at_upper(), 2 );
		TS_ASSERT_DELTA( pose.residue(1).xyz("C").distance( pose.residue(2).xyz("N") ), 1.328685, 0.01 ); //Check peptide bond length.
	}

	void test_multi_prepend() {
		core::pose::Pose pose( *(pose_) ); //Local copy

		core::Size const original_pose_size( pose.total_residue() );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		for ( core::Size i(0); i<5; ++i ) {
			stubmover.add_residue( "Prepend", "ALA", 0, false, "", 0, 1, "" );
		}
		stubmover.apply(pose);
		TS_ASSERT_EQUALS(pose.total_residue(), original_pose_size+5 );
		TS_ASSERT_EQUALS( pose.residue_type( pose.total_residue() ).name3(), "VAL" );
		TS_ASSERT_EQUALS( pose.residue_type( 6 ).name3(), "VAL" );
		for ( core::Size i(1); i<=5; ++i ) {
			TS_ASSERT_EQUALS( pose.residue_type( i ).name3(), "ALA" );
			TS_ASSERT_EQUALS( pose.residue( i+1 ).connected_residue_at_lower(), i );
			TS_ASSERT_EQUALS( pose.residue( i ).connected_residue_at_upper(), i+1 );
			TS_ASSERT_DELTA( pose.residue(i).xyz("C").distance( pose.residue(i+1).xyz("N") ), 1.328685, 0.01 ); //Check peptide bond length.
		}
	}

	void test_prepend_repeat() {
		core::pose::Pose pose( *(pose_) ); //Local copy

		core::Size const original_pose_size( pose.total_residue() );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Prepend", "ALA", 0, false, "", 5, 1, "" );
		stubmover.apply(pose);
		TS_ASSERT_EQUALS(pose.total_residue(), original_pose_size+5 );
		TS_ASSERT_EQUALS( pose.residue_type( pose.total_residue() ).name3(), "VAL" );
		TS_ASSERT_EQUALS( pose.residue_type( 6 ).name3(), "VAL" );
		for ( core::Size i(1); i<=5; ++i ) {
			TS_ASSERT_EQUALS( pose.residue_type( i ).name3(), "ALA" );
			TS_ASSERT_EQUALS( pose.residue( i+1 ).connected_residue_at_lower(), i );
			TS_ASSERT_EQUALS( pose.residue( i ).connected_residue_at_upper(), i+1 );

			TS_ASSERT_DELTA( pose.residue(i).xyz("C").distance( pose.residue(i+1).xyz("N") ), 1.328685, 0.01 ); //Check peptide bond length.
		}
	}

	void test_prepend_fail() {
		core::pose::Pose pose( *(pose_) ); //Local copy
		core::Size const original_pose_size( pose.total_residue() );

		basic::datacache::DataMap data;
		Filters_map filters;
		Movers_map movers;
		TagCOP tag = tagptr_from_string("<PeptideStubMover name=\"psm\" reset=\"false\">\n<Append resname=\"ALA\" anchor_rsd=\"72\"/>\n <Prepend resname=\"ALA\" anchor_rsd=\"1\" repeat=\"3\"/></PeptideStubMover>\n");

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.parse_my_tag( tag, data, filters, movers, pose );
		stubmover.apply(pose);
		TS_ASSERT_EQUALS( pose.total_residue(), original_pose_size + 4 );
		TS_ASSERT_EQUALS( pose.residue_type( 1 ).name3(), "ALA" );

	}


private:

	core::pose::PoseCOP pose_;



};
