// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/core/pose/Pose.cxxtest.hh
/// @brief  unit tests for core::pose::Pose
/// @author

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/variant_util.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/signals/ConformationEvent.hh>
#include <core/pose/signals/DestructionEvent.hh>
#include <core/pose/signals/EnergyEvent.hh>
#include <core/pose/signals/GeneralEvent.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <test/util/pose_funcs.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/angle.functions.hh>
#include <numeric/constants.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>

// C++ headers
#include <iostream>

#include <basic/Tracer.hh>

#ifdef SERIALIZATION
// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR("core.pose.Pose.cxxtest");

namespace test_pose {

// faux observer for Pose observer interface
struct Obs {
	typedef core::Size Size;
	typedef core::pose::signals::DestructionEvent DestructionEvent;
	typedef core::pose::signals::GeneralEvent GeneralEvent;
	typedef core::pose::signals::ConformationEvent ConformationEvent;
	typedef core::pose::signals::EnergyEvent EnergyEvent;

	Obs() : count( 0 ), g_count( 0 ) {}

	void on_destruction_change( DestructionEvent const & ) { ++count; }
	void on_general_change( GeneralEvent const & ) { ++g_count; }
	void on_conformation_change( ConformationEvent const & ) { ++count; }
	void on_energy_change( EnergyEvent const & ) { ++count; }

	int count;
	int g_count;
};

} // namespace test_pose


class PoseTests : public CxxTest::TestSuite
{


public: //setup


	PoseTests() {};


	// Shared initialization.
	void setUp() {
		core_init_with_additional_options( "-output_virtual -write_all_connect_info" );
	}


	// Shared finalization.
	void tearDown() {
	}


public: // tests

	/// @brief Test whether mainchain torsions are measured correctly
	/// for single- and multi-chain poses (in the latter case, with covalent
	/// bonds between the chains).
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_Pose_mainchain_torsions() {
		core::Real const deltathresh_weak(0.06);
		core::Real const deltathresh_strong(0.0001);

		core::pose::Pose onechainpose;
		core::pose::Pose multichainpose;
		core::import_pose::pose_from_file( onechainpose, "core/pose/onechain.pdb", core::import_pose::PDB_file);
		core::import_pose::pose_from_file( multichainpose, "core/pose/multichain.pdb", core::import_pose::PDB_file);

		core::pose::remove_variant_type_from_pose_residue( multichainpose, core::chemical::UPPER_TERMINUS_VARIANT, 7 );
		core::pose::remove_variant_type_from_pose_residue( multichainpose, core::chemical::CUTPOINT_LOWER, 7 );
		core::pose::remove_variant_type_from_pose_residue( multichainpose, core::chemical::CUTPOINT_UPPER, 7 );
		core::pose::remove_variant_type_from_pose_residue( multichainpose, core::chemical::LOWER_TERMINUS_VARIANT, 8 );
		core::pose::remove_variant_type_from_pose_residue( multichainpose, core::chemical::CUTPOINT_LOWER, 7 );
		core::pose::remove_variant_type_from_pose_residue( multichainpose, core::chemical::CUTPOINT_UPPER, 7 );

		core::pose::remove_variant_type_from_pose_residue( multichainpose, core::chemical::UPPER_TERMINUS_VARIANT, 14 );
		core::pose::remove_variant_type_from_pose_residue( multichainpose, core::chemical::CUTPOINT_LOWER, 14 );
		core::pose::remove_variant_type_from_pose_residue( multichainpose, core::chemical::CUTPOINT_UPPER, 14 );
		core::pose::remove_variant_type_from_pose_residue( multichainpose, core::chemical::LOWER_TERMINUS_VARIANT, 15 );
		core::pose::remove_variant_type_from_pose_residue( multichainpose, core::chemical::CUTPOINT_LOWER, 15 );
		core::pose::remove_variant_type_from_pose_residue( multichainpose, core::chemical::CUTPOINT_UPPER, 15 );

		core::pose::remove_variant_type_from_pose_residue( multichainpose, core::chemical::UPPER_TERMINUS_VARIANT, 20 );
		core::pose::remove_variant_type_from_pose_residue( multichainpose, core::chemical::CUTPOINT_LOWER, 20 );
		core::pose::remove_variant_type_from_pose_residue( multichainpose, core::chemical::CUTPOINT_UPPER, 20 );
		core::pose::remove_variant_type_from_pose_residue( multichainpose, core::chemical::LOWER_TERMINUS_VARIANT, 21 );
		core::pose::remove_variant_type_from_pose_residue( multichainpose, core::chemical::CUTPOINT_LOWER, 21 );
		core::pose::remove_variant_type_from_pose_residue( multichainpose, core::chemical::CUTPOINT_UPPER, 21 );

		multichainpose.conformation().declare_chemical_bond(7, "C", 8, "N");
		multichainpose.conformation().declare_chemical_bond(14, "C", 15, "N");
		multichainpose.conformation().declare_chemical_bond(20, "C", 21, "N");
		multichainpose.conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(7);
		multichainpose.conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(8);
		multichainpose.conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(14);
		multichainpose.conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(15);
		multichainpose.conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(20);
		multichainpose.conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(21);

		utility::vector1< core::Real > expected_phi;
		expected_phi.reserve(31);
		expected_phi.push_back( 0.0 );
		expected_phi.push_back( -135.0 );
		expected_phi.push_back( -135.1 );
		expected_phi.push_back( -134.9 );
		expected_phi.push_back( -135.1 );
		expected_phi.push_back( 38.6 );
		expected_phi.push_back( 52.9 );
		expected_phi.push_back( -43.0 );
		expected_phi.push_back( 1.5 );
		expected_phi.push_back( -135.0 );
		expected_phi.push_back( -134.9 );
		expected_phi.push_back( -135.1 );
		expected_phi.push_back( -122.9 );
		expected_phi.push_back( -36.1 );
		expected_phi.push_back( -167.8 );
		expected_phi.push_back( -96.1 );
		expected_phi.push_back( -135.1 );
		expected_phi.push_back( -134.9 );
		expected_phi.push_back( -135.1 );
		expected_phi.push_back( -86.2 );
		expected_phi.push_back( -77.3 );
		expected_phi.push_back( -34.6 );
		expected_phi.push_back( 115.1 );
		expected_phi.push_back( -83.7 );
		expected_phi.push_back( -165.6 );
		expected_phi.push_back( -58.7 );
		expected_phi.push_back( 107.9 );
		expected_phi.push_back( -135.1 );
		expected_phi.push_back( -135.0 );
		expected_phi.push_back( -135.1 );
		expected_phi.push_back( -134.9 );

		TR << "RES\tPHI_EXP\tPHI_1chain\tPHI_nchain" << std::endl;
		for ( core::Size ir=1; ir<=31; ++ir ) {
			TR << ir << "\t" << expected_phi[ir] << "\t" << onechainpose.phi(ir) << "\t" << multichainpose.phi(ir) << std::endl;
			TS_ASSERT_DELTA(expected_phi[ir], onechainpose.phi(ir), deltathresh_weak);
			TS_ASSERT_DELTA(expected_phi[ir], multichainpose.phi(ir), deltathresh_weak);
			TS_ASSERT_DELTA(onechainpose.phi(ir), multichainpose.phi(ir), deltathresh_strong);
		}
		TR << std::endl;

		utility::vector1< core::Real > expected_psi;
		expected_psi.reserve(31);
		expected_psi.push_back( 135.1 );
		expected_psi.push_back( 134.9 );
		expected_psi.push_back( 135.1 );
		expected_psi.push_back( 135.0 );
		expected_psi.push_back( 135.1 );
		expected_psi.push_back( 70.2 );
		expected_psi.push_back( -155.2 );
		expected_psi.push_back( -33.0 );
		expected_psi.push_back( -75.2 );
		expected_psi.push_back( 135.1 );
		expected_psi.push_back( 134.9 );
		expected_psi.push_back( 135.1 );
		expected_psi.push_back( -149.5 );
		expected_psi.push_back( -39.3 );
		expected_psi.push_back( -8.4 );
		expected_psi.push_back( 139.2 );
		expected_psi.push_back( 135.1 );
		expected_psi.push_back( 135.0 );
		expected_psi.push_back( 135.0 );
		expected_psi.push_back( 165.8 );
		expected_psi.push_back( 9.1 );
		expected_psi.push_back( -38.9 );
		expected_psi.push_back( 163.3 );
		expected_psi.push_back( 138.0 );
		expected_psi.push_back( 156.0 );
		expected_psi.push_back( 159.5 );
		expected_psi.push_back(  71.3 );
		expected_psi.push_back( 135.1 );
		expected_psi.push_back( 134.9 );
		expected_psi.push_back( 135.1 );
		expected_psi.push_back( 0.0 );

		TR << "RES\tPSI_EXP\tPSI_1chain\tPSI_nchain" << std::endl;
		for ( core::Size ir=1; ir<=31; ++ir ) {
			TR << ir << "\t" << expected_psi[ir] << "\t" << onechainpose.psi(ir) << "\t" << multichainpose.psi(ir) << std::endl;
			TS_ASSERT_DELTA(expected_psi[ir], onechainpose.psi(ir), deltathresh_weak);
			TS_ASSERT_DELTA(expected_psi[ir], multichainpose.psi(ir), deltathresh_weak);
			TS_ASSERT_DELTA(onechainpose.psi(ir), multichainpose.psi(ir), deltathresh_strong);
		}
		TR << std::endl;
	}

	/// @brief test Pose observer interface
	void test_Pose_observer() {
		using namespace core::scoring;
		using core::pose::Pose;

		using namespace test_pose;

		ScoreFunctionOP scorefxn = get_score_function();

		Pose * pose = new Pose;
		core::import_pose::pose_from_file( *pose, "core/pose/pdbinfo_test_in.pdb" , core::import_pose::PDB_file);

		// attach observer
		Obs obs;
		pose->attach_destruction_obs( &Obs::on_destruction_change, &obs );
		pose->attach_general_obs( &Obs::on_general_change, &obs );
		pose->attach_conformation_obs( &Obs::on_conformation_change, &obs );
		pose->attach_energy_obs( &Obs::on_energy_change, &obs );

		// test ConformationEvent & GeneralEvent
		pose->set_phi( 3, 45 );
		pose->residue( 2 ); // force refold
		TS_ASSERT_EQUALS( obs.count, 1 );
		TS_ASSERT_EQUALS( obs.g_count, 1 );
		obs.count = 0;
		obs.g_count = 0;

		// test EnergyEvent & GeneralEvent
		scorefxn->score( *pose );
		TS_ASSERT_EQUALS( obs.count, 1 );
		TS_ASSERT_EQUALS( obs.g_count, 1 );
		obs.count = 0;
		obs.g_count = 0;

		// test DestructionEvent
		delete pose;
		TS_ASSERT_EQUALS( obs.count, 1 );
		obs.count = 0;
	}

	void test_append_pose_by_jump()
	{
		using namespace core::pose;

		Pose pose;
		make_pose_from_sequence(pose, "TE", "fa_standard");
		pose.pdb_info(PDBInfoOP( new PDBInfo(pose) ));

		Pose pose2;
		make_pose_from_sequence(pose2, "ST", "fa_standard");
		pose2.pdb_info(PDBInfoOP( new PDBInfo(pose2) ));

		char pose2_chain = 'B';
		pose2.pdb_info()->set_chains(pose2_chain);

		Pose work_pose(pose);
		work_pose.append_pose_by_jump(pose2, 1);

		TS_ASSERT_EQUALS(work_pose.size(), pose.size() + pose2.size());
		TS_ASSERT_EQUALS(work_pose.sequence(), pose.sequence() + pose2.sequence());
		TS_ASSERT_EQUALS(work_pose.sequence(2, 3), "ES");

		for ( core::Size i = 1; i <= pose.size(); i++ ) {
			TS_ASSERT_DELTA(
				work_pose.residue(i).xyz(1),
				pose.residue(i).xyz(1),
				1e-6);
			TS_ASSERT_EQUALS(
				work_pose.pdb_info()->chain(i),
				pose.pdb_info()->chain(i));
		}

		for ( core::Size i = 1; i <= pose2.size(); i++ ) {
			TS_ASSERT_DELTA(
				work_pose.residue(i + pose.size()).xyz(1),
				pose2.residue(i).xyz(1),
				1e-6);
			TS_ASSERT_EQUALS(
				work_pose.pdb_info()->chain(i + pose.size()),
				pose2.pdb_info()->chain(i));
		}
	}

	void test_update_residue_neighbors() {
		// These are tests mainly just to make sure that the calls don't crash.
		using namespace core::pose;

		Pose empty;
		empty.update_residue_neighbors();

		Pose pose;
		core::import_pose::pose_from_file( pose, "core/pose/pdbinfo_test_in.pdb" , core::import_pose::PDB_file);
		pose.update_residue_neighbors();
	}

	/// @brief Test that the cutpoint virtuals are updated properly when the cutpoint variant type is added.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void test_update_cutpoint_virtual_atoms_when_cutpoint_variant_added() {
		core::pose::Pose pose;
		make_pose_from_sequence( pose, "AAAAAAAAAA", "fa_standard" );
		for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
			pose.set_phi(i,-135);
			pose.set_psi(i,135);
			pose.set_omega(i,180);
		}
		core::kinematics::FoldTree ft;
		ft.add_edge( 1, 5, core::kinematics::Edge::PEPTIDE );
		ft.add_edge( 1, 10, 1  );
		ft.add_edge( 10, 6, core::kinematics::Edge::PEPTIDE );
		pose.fold_tree(ft);

		//Translate part of the pose
		numeric::xyzVector< core::Real > transvect(2, 7, 12);
		for ( core::Size ir(6), irmax(pose.total_residue()); ir<=irmax; ++ir ) {
			for ( core::Size ia(1), iamax(pose.residue_type(ir).natoms()); ia<=iamax; ++ia ) {
				core::id::AtomID curat( ia, ir );
				pose.set_xyz( curat, transvect + pose.xyz( curat ) );
			}
		}
		//Arbitrarily move alpha carbon to make non-ideal geometry.
		numeric::xyzVector< core::Real > randvect1(0.5,0.7,0.25);
		pose.set_xyz(core::id::NamedAtomID("CA", 6), pose.xyz(core::id::NamedAtomID("CA",6)) + randvect1);
		pose.update_residue_neighbors();

		core::Real const N_CA_dist( pose.residue(6).xyz("N").distance( pose.residue(6).xyz("CA") ) );
		core::Real const LOWER_N_dist( pose.residue(6).lower_connect().icoor().d() );
		core::Real const LOWER_N_angle( pose.residue(6).lower_connect().icoor().theta() );
		core::Real const C_UPPER_dist( pose.residue(5).upper_connect().icoor().d() );

		//Add the modified variant types.  Should call update_cutpoint_virtual_atoms_if_connected():
		core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_LOWER, 5 );
		//pose.dump_pdb("vtemp5.pdb"); //DELETE ME
		core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_UPPER, 6 );

		//Check that virtuals have been updated properly:
		core::Real const C_OVL1_dist( pose.residue(5).xyz("OVL1").distance( pose.residue(5).xyz("C") ) );
		core::Real const OVL1_OVL2_dist( pose.residue(5).xyz("OVL2").distance( pose.residue(5).xyz("OVL1") ) );
		core::Real const OVL2_angle( numeric::constants::d::pi - numeric::angle_radians( pose.residue(5).xyz("C"), pose.residue(5).xyz("OVL1"), pose.residue(5).xyz("OVL2") ) );
		core::Real const OVU1_N_dist( pose.residue(6).xyz("OVU1").distance( pose.residue(6).xyz("N") ) );

		TR << "\nN_CA_dist:\t" << N_CA_dist << "\n";
		TR << "OVL1_OVL2_dist:\t" << OVL1_OVL2_dist << "\n";
		TR << "LOWER_N_dist:\t" << LOWER_N_dist << "\n";
		TR << "C_OVL1_dist:\t" << C_OVL1_dist << "\n";
		TR << "LOWER_N_angle:\t" << LOWER_N_angle << "\n";
		TR << "OVL2_angle:\t" << OVL2_angle << "\n";
		TR << "C_UPPER_dist:\t" << C_UPPER_dist << "\n";
		TR << "OVU1_N_dist:\t" << OVU1_N_dist << "\n";
		TR.flush();

		TS_ASSERT_DELTA( N_CA_dist, OVL1_OVL2_dist, 0.001 );
		TS_ASSERT_DELTA( LOWER_N_dist, C_OVL1_dist, 0.001 );
		TS_ASSERT_DELTA( LOWER_N_angle, OVL2_angle, 0.001 );
		TS_ASSERT_DELTA( C_UPPER_dist, OVU1_N_dist, 0.001 );

		//pose.dump_pdb("vtemp6.pdb"); //DELETE ME
	}

	/// @brief Test that cutpoint virtuals are updated properly when a bond is declared.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void test_update_cutpoint_virtual_atoms_when_bond_declared() {
		core::pose::Pose pose;
		make_pose_from_sequence( pose, "AAAAAAAAAA", "fa_standard" );
		for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
			pose.set_phi(i,-135);
			pose.set_psi(i,135);
			pose.set_omega(i,180);
		}
		core::kinematics::FoldTree ft;
		ft.add_edge( 1, 5, core::kinematics::Edge::PEPTIDE );
		ft.add_edge( 1, 10, 1  );
		ft.add_edge( 10, 6, core::kinematics::Edge::PEPTIDE );
		pose.fold_tree(ft);

		//Translate part of the pose
		numeric::xyzVector< core::Real > transvect(2, 7, 12);
		for ( core::Size ir(6), irmax(pose.total_residue()); ir<=irmax; ++ir ) {
			for ( core::Size ia(1), iamax(pose.residue_type(ir).natoms()); ia<=iamax; ++ia ) {
				core::id::AtomID curat( ia, ir );
				pose.set_xyz( curat, transvect + pose.xyz( curat ) );
			}
		}
		//Arbitrarily move alpha carbon to make non-ideal geometry.
		numeric::xyzVector< core::Real > randvect1(0.5,0.7,0.25);
		pose.set_xyz(core::id::NamedAtomID("CA", 6), pose.xyz(core::id::NamedAtomID("CA",6)) + randvect1);
		pose.update_residue_neighbors();

		core::Real const N_CA_dist( pose.residue(6).xyz("N").distance( pose.residue(6).xyz("CA") ) );
		core::Real const LOWER_N_dist( pose.residue(6).lower_connect().icoor().d() );
		core::Real const LOWER_N_angle( pose.residue(6).lower_connect().icoor().theta() );
		core::Real const C_UPPER_dist( pose.residue(5).upper_connect().icoor().d() );

		//Make sure that there's no bond there:
		pose.conformation().sever_chemical_bond( 5, pose.residue_type(5).upper_connect_id(), 6, pose.residue_type(6).lower_connect_id() );

		//Add the modified variant types WITHOUT adding a chemical bond:
		core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_LOWER, 5 );
		core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_UPPER, 6 );
		//pose.dump_pdb("vtemp3.pdb"); //DELETE ME

		//Re-declare the bond:
		pose.conformation().declare_chemical_bond( 5, "C", 6, "N" );

		//Check that virtuals have been updated properly:
		core::Real const C_OVL1_dist( pose.residue(5).xyz("OVL1").distance( pose.residue(5).xyz("C") ) );
		core::Real const OVL1_OVL2_dist( pose.residue(5).xyz("OVL2").distance( pose.residue(5).xyz("OVL1") ) );
		core::Real const OVL2_angle( numeric::constants::d::pi - numeric::angle_radians( pose.residue(5).xyz("C"), pose.residue(5).xyz("OVL1"), pose.residue(5).xyz("OVL2") ) );
		core::Real const OVU1_N_dist( pose.residue(6).xyz("OVU1").distance( pose.residue(6).xyz("N") ) );

		TR << "\nN_CA_dist:\t" << N_CA_dist << "\n";
		TR << "OVL1_OVL2_dist:\t" << OVL1_OVL2_dist << "\n";
		TR << "LOWER_N_dist:\t" << LOWER_N_dist << "\n";
		TR << "C_OVL1_dist:\t" << C_OVL1_dist << "\n";
		TR << "LOWER_N_angle:\t" << LOWER_N_angle << "\n";
		TR << "OVL2_angle:\t" << OVL2_angle << "\n";
		TR << "C_UPPER_dist:\t" << C_UPPER_dist << "\n";
		TR << "OVU1_N_dist:\t" << OVU1_N_dist << "\n";
		TR.flush();

		TS_ASSERT_DELTA( N_CA_dist, OVL1_OVL2_dist, 0.001 );
		TS_ASSERT_DELTA( LOWER_N_dist, C_OVL1_dist, 0.001 );
		TS_ASSERT_DELTA( LOWER_N_angle, OVL2_angle, 0.001 );
		TS_ASSERT_DELTA( C_UPPER_dist, OVU1_N_dist, 0.001 );

		//pose.dump_pdb("vtemp4.pdb"); //DELETE ME

	}

	/// @brief Test the function in core/pose/variant_util.hh that updates the coordinates of virtual
	/// atoms in cutpoint residues.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void test_update_cutpoint_virtual_atoms_if_connected() {
		core::pose::Pose pose;
		make_pose_from_sequence( pose, "AAAAAAAAAA", "fa_standard" );
		for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
			pose.set_phi(i,-135);
			pose.set_psi(i,135);
			pose.set_omega(i,180);
		}
		core::kinematics::FoldTree ft;
		ft.add_edge( 1, 5, core::kinematics::Edge::PEPTIDE );
		ft.add_edge( 1, 10, 1  );
		ft.add_edge( 10, 6, core::kinematics::Edge::PEPTIDE );
		pose.fold_tree(ft);

		//Translate part of the pose
		numeric::xyzVector< core::Real > transvect(2, 7, 12);
		for ( core::Size ir(6), irmax(pose.total_residue()); ir<=irmax; ++ir ) {
			for ( core::Size ia(1), iamax(pose.residue_type(ir).natoms()); ia<=iamax; ++ia ) {
				core::id::AtomID curat( ia, ir );
				pose.set_xyz( curat, transvect + pose.xyz( curat ) );
			}
		}
		//Arbitrarily move alpha carbon to make non-ideal geometry.
		numeric::xyzVector< core::Real > randvect1(0.5,0.7,0.25);
		pose.set_xyz(core::id::NamedAtomID("CA", 6), pose.xyz(core::id::NamedAtomID("CA",6)) + randvect1);
		pose.update_residue_neighbors();

		core::Real const N_CA_dist( pose.residue(6).xyz("N").distance( pose.residue(6).xyz("CA") ) );
		core::Real const LOWER_N_dist( pose.residue(6).lower_connect().icoor().d() );
		core::Real const LOWER_N_angle( pose.residue(6).lower_connect().icoor().theta() );
		core::Real const C_UPPER_dist( pose.residue(5).upper_connect().icoor().d() );

		core::pose::correctly_add_cutpoint_variants( pose, 5 ); //Calls the update_cutpoint_virtual_atoms_if_connected() function.

		core::Real const C_OVL1_dist( pose.residue(5).xyz("OVL1").distance( pose.residue(5).xyz("C") ) );
		core::Real const OVL1_OVL2_dist( pose.residue(5).xyz("OVL2").distance( pose.residue(5).xyz("OVL1") ) );
		core::Real const OVL2_angle( numeric::constants::d::pi - numeric::angle_radians( pose.residue(5).xyz("C"), pose.residue(5).xyz("OVL1"), pose.residue(5).xyz("OVL2") ) );
		core::Real const OVU1_N_dist( pose.residue(6).xyz("OVU1").distance( pose.residue(6).xyz("N") ) );

		TR << "\nN_CA_dist:\t" << N_CA_dist << "\n";
		TR << "OVL1_OVL2_dist:\t" << OVL1_OVL2_dist << "\n";
		TR << "LOWER_N_dist:\t" << LOWER_N_dist << "\n";
		TR << "C_OVL1_dist:\t" << C_OVL1_dist << "\n";
		TR << "LOWER_N_angle:\t" << LOWER_N_angle << "\n";
		TR << "OVL2_angle:\t" << OVL2_angle << "\n";
		TR << "C_UPPER_dist:\t" << C_UPPER_dist << "\n";
		TR << "OVU1_N_dist:\t" << OVU1_N_dist << "\n";
		TR.flush();

		TS_ASSERT_DELTA( N_CA_dist, OVL1_OVL2_dist, 0.001 );
		TS_ASSERT_DELTA( LOWER_N_dist, C_OVL1_dist, 0.001 );
		TS_ASSERT_DELTA( LOWER_N_angle, OVL2_angle, 0.001 );
		TS_ASSERT_DELTA( C_UPPER_dist, OVU1_N_dist, 0.001 );

		//pose.dump_pdb("vtemp2.pdb"); //DELETE ME
	}

	/// @brief Test that cutpoints allow free rotation of the virtual atoms.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void test_cutpoint_bond_rotation() {
		core::pose::Pose pose;
		make_pose_from_sequence( pose, "AAAAAAAAAA", "fa_standard" );
		for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
			pose.set_phi(i,-135);
			pose.set_psi(i,135);
			pose.set_omega(i,180);
		}
		core::kinematics::FoldTree ft;
		ft.add_edge( 1, 5, core::kinematics::Edge::PEPTIDE );
		ft.add_edge( 1, 10, 1  );
		ft.add_edge( 10, 6, core::kinematics::Edge::PEPTIDE );
		pose.fold_tree(ft);
		core::pose::correctly_add_cutpoint_variants( pose, 5 );

		TR << "\nOMEGA\tEXPECTED\n";
		for ( core::Size i(1); i<=36; ++i ) {
			pose.set_omega( 5, static_cast<core::Real>(i)*10 - 5 );
			/*char outfile [256];
			sprintf( outfile, "out1_%04lu.pdb", i );
			pose.dump_pdb( std::string(outfile) );*/
			core::Real const measured_angle( numeric::dihedral_degrees(
				pose.residue(5).xyz("CA"),
				pose.residue(5).xyz("C"),
				pose.residue(5).xyz("OVL1"),
				pose.residue(5).xyz("OVL2")
				) );
			TR << numeric::nonnegative_principal_angle_degrees(measured_angle) << "\t" << numeric::nonnegative_principal_angle_degrees( static_cast<core::Real>(i)*10 - 5 ) << "\n";
			TS_ASSERT_DELTA( numeric::nonnegative_principal_angle_degrees( measured_angle ), numeric::nonnegative_principal_angle_degrees( static_cast<core::Real>(i)*10 - 5 ), 0.01 );
		}

		TR << "\nPHI\tEXPECTED\n";
		for ( core::Size i(1); i<=36; ++i ) {
			pose.set_phi( 6, static_cast<core::Real>(i)*10 - 5 );
			/*char outfile [256];
			sprintf( outfile, "out2_%04lu.pdb", i );
			pose.dump_pdb( std::string(outfile) );*/
			core::Real const measured_angle( numeric::dihedral_degrees(
				pose.residue(6).xyz("OVU1"),
				pose.residue(6).xyz("N"),
				pose.residue(6).xyz("CA"),
				pose.residue(6).xyz("C")
				) );
			TR << numeric::nonnegative_principal_angle_degrees(measured_angle) << "\t" << numeric::nonnegative_principal_angle_degrees( static_cast<core::Real>(i)*10 - 5 ) << "\n";
			TS_ASSERT_DELTA( numeric::nonnegative_principal_angle_degrees( measured_angle ), numeric::nonnegative_principal_angle_degrees( static_cast<core::Real>(i)*10 - 5 ), 0.01 );
		}
		TR.flush();
	}

	void test_serialize_pose() {
		TS_ASSERT( true );
#ifdef SERIALIZATION
		using namespace core::conformation;
		using namespace core::pose;
		using namespace core::id;
		using namespace core::kinematics;
		using namespace core::scoring::constraints;
		using namespace core::scoring::func;

		// TO DO: Test constraint set propertly observes the deserialized pose (i.e. following residue insertion)
		// TO DO: Test CacheableObservers properly observe the deserialized pose
		// TO DO: Test that a symmetric conformation is properly deserialized inside a pose

		PoseOP trpcage = create_trpcage_ideal_poseop();
		HarmonicFuncOP harm_0_1( new HarmonicFunc( 0, 1 ) );
		ConstraintOP apc_1_10( new AtomPairConstraint( AtomID(3,1), AtomID(3,10), harm_0_1 ));
		ConstraintOP apc_1_11( new AtomPairConstraint( AtomID(3,1), AtomID(3,11), harm_0_1 ));
		trpcage->add_constraint( apc_1_10 );

		// Now serialize the Pose
		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arch( oss );
			arch( trpcage );
		}

		PoseOP trpcage_copy;
		std::istringstream iss( oss.str() );
		{
			cereal::BinaryInputArchive arch( iss );
			arch( trpcage_copy );
		}

		TS_ASSERT( trpcage->total_residue() == trpcage_copy->total_residue() );

		// Make sure constraint set was serialized properly
		TS_ASSERT( ! trpcage_copy->constraint_set()->is_empty() );
		TS_ASSERT( trpcage_copy->constraint_set()->get_all_constraints().size() == 1 );
		TS_ASSERT( *trpcage_copy->constraint_set()->get_all_constraints()[1] == *apc_1_10 );
		// make sure that constraint set is observing the conformation
		// do something silly by making another copy of residue 15 and inserting it between
		// residue 5 and residue 6 (set a foldtree that allows this insertion first)
		FoldTree ft;
		ft.add_edge( 1, 5, Edge::PEPTIDE );  ft.add_edge( 5, 10, 1  );
		ft.add_edge( 10, 6, Edge::PEPTIDE ); ft.add_edge( 10, 20, Edge::PEPTIDE );
		trpcage_copy->fold_tree( ft );
		trpcage_copy->insert_residue_by_bond( trpcage_copy->residue(15), 6, 5 );
		TS_ASSERT( ! ( *trpcage_copy->constraint_set()->get_all_constraints()[1] == *apc_1_10 ));
		TS_ASSERT(     *trpcage_copy->constraint_set()->get_all_constraints()[1] == *apc_1_11 );

		// Make sure that the deserialized PDBInfo object is reestablished as observing the
		// Pose's Conformation
		core::conformation::ConformationCOP obs_by_pdb_info = trpcage_copy->pdb_info()->is_observing().lock();
		TS_ASSERT_EQUALS( obs_by_pdb_info.get(), &trpcage_copy->conformation() );


#endif
	}

};
