// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/Conformation.cxxtest.hh
/// @brief  test suite for Minimizer
/// @author Phil Bradley
/// @author Sergey Lyskov

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/signals/ConnectionEvent.hh>
#include <core/conformation/signals/GeneralEvent.hh>
#include <core/conformation/signals/IdentityEvent.hh>
#include <core/conformation/signals/LengthEvent.hh>
#include <core/conformation/signals/XYZEvent.hh>
#include <core/id/TorsionID.hh>
#include <core/conformation/carbohydrates/GlycanTreeSetObserver.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/types.hh>
#include <core/pose/variant_util.hh>

#include <test/util/pose_funcs.hh>
#include <test/UTracer.hh>

//Auto Headers
#include <core/id/AtomID_Mask.hh>
#include <core/kinematics/AtomPointer.hh>
#include <core/kinematics/DomainMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <utility/vector1.hh>
#include <utility/io/all.fwd.hh>
#include <utility/io/icstream.fwd.hh>
#include <utility/io/ocstream.fwd.hh>
#include <numeric/BodyPosition.fwd.hh>
#include <numeric/Quaternion.fwd.hh>
#include <numeric/all.fwd.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.conformation.Conformation.cxxtest");

using namespace core;


namespace test_conf {

// faux observer for Conformation observer interface
struct Obs {
	typedef core::Size Size;
	typedef core::conformation::signals::ConnectionEvent ConnectionEvent;
	typedef core::conformation::signals::GeneralEvent GeneralEvent;
	typedef core::conformation::signals::IdentityEvent IdentityEvent;
	typedef core::conformation::signals::LengthEvent LengthEvent;
	typedef core::conformation::signals::XYZEvent XYZEvent;

	Obs() {
		reset();
	}

	void on_connection_change( ConnectionEvent const & ) { c_last = ++count; }
	void on_general_change( GeneralEvent const & ) { ++g_count; }
	void on_identity_change( IdentityEvent const & e ) { i_event = e; i_last = ++count; }
	void on_length_change( LengthEvent const & e ) { l_event = e; l_last = ++count; }
	void on_xyz_change( XYZEvent const & ) { x_last = ++count; }
	void reset() {
		count = 0;
		g_count = 0;
		c_last = 0;
		i_last = 0;
		l_last = 0;
		x_last = 0;
	}
	int count;
	int g_count;
	IdentityEvent i_event;
	LengthEvent l_event;
	// The count at the last respective event.
	int c_last, i_last, l_last, x_last;
};

} // test_conf


class ConformationTests : public CxxTest::TestSuite
{
	chemical::ResidueTypeSetCAP residue_set;

public:
	ConformationTests() {};

	// Shared initialization goes here.
	void setUp() {
		//core_init();
		core_init_with_additional_options("-extra_res_fa core/chemical/ASX.params core/chemical/LYX.params"); // Remove if database versions get permanently enabled

		residue_set = chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	///////////////////////////////////////////////////////////////////////////////
	// ------------------------------------------ //
	/// @brief simple Conformation test
	void test_simple_conformation()
	{
		using namespace io::pdb;
		using namespace conformation;
		using namespace chemical;

		test::UTracer UT("core/conformation/test_simple_conformation.u");

		//pose::Pose pose;
		//pose_from_file( pose, "core/conformation/test_in.pdb" , core::import_pose::PDB_file);
		pose::Pose pose( create_test_in_pdb_pose());

		kinematics::FoldTree f( pose.size() );
		f.new_jump( 8, 26, 18 );
		f.reorder( 8 );
		pose.fold_tree( f );

		dump_pdb( pose, "core/conformation/start.pdb.new" );
		UTRACE << pose.fold_tree() << "\n";

		pose.conformation().delete_polymer_residue( 1 );

		dump_pdb( pose, "core/conformation/del_0.pdb.new" );
		UTRACE << pose.fold_tree() << "\n";

		// prepend an alanine residue
		ResidueTypeSetCOP residue_set
			( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
		ResidueOP ala_rsd( ResidueFactory::create_residue( residue_set->name_map( "ALA" ) ) );

		pose.conformation().prepend_polymer_residue_before_seqpos( *ala_rsd, 1, true );
		pose.set_omega(1,180);

		dump_pdb( pose, "core/conformation/ala_1.pdb.new" );
		UTRACE << pose.fold_tree() << "\n";

		pose.conformation().prepend_polymer_residue_before_seqpos( *ala_rsd, 1, true );
		pose.set_omega(1,180);

		dump_pdb( pose, "core/conformation/ala_2.pdb.new" );
		UTRACE << pose.fold_tree() << "\n";

		pose.conformation().delete_polymer_residue( pose.fold_tree().cutpoint(1) );

		dump_pdb( pose, "core/conformation/del_1.pdb.new" );
		UTRACE << pose.fold_tree() << "\n";

		pose.conformation().delete_polymer_residue( pose.fold_tree().cutpoint(1) );

		dump_pdb( pose, "core/conformation/del_2.pdb.new" );
		UTRACE << pose.fold_tree() << "\n";

		pose.conformation().append_polymer_residue_after_seqpos( *ala_rsd, pose.fold_tree().cutpoint(1), true );

		dump_pdb( pose, "core/conformation/ala_3.pdb.new" );
		UTRACE << pose.fold_tree() << "\n";
	} // test_simple_conformation


	/// @brief test Conformation observer interface
	void test_Conformation_observer() {
		// project namespaces
		using namespace core::io::pdb;
		using namespace core::conformation;
		using namespace core::chemical;
		using core::pose::Pose;

		// test namespaces
		using test_conf::Obs;

		// grab alanine residue
		ResidueTypeSetCOP residue_set( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
		ResidueOP ala_rsd( ResidueFactory::create_residue( residue_set->name_map( "ALA" ) ) );

		//Pose pose;
		//pose_from_file( pose, "core/conformation/test_in.pdb" , core::import_pose::PDB_file);
		pose::Pose pose( create_test_in_pdb_pose());
		ConformationOP conf_op = pose.conformation().clone();
		Conformation & conf = *conf_op;
		Obs obs;

		// connection observers
		{ // begin scope
			ConformationOP c2_op = conf.clone();
			Conformation & c2 = *c2_op;
			TS_ASSERT( c2.attach_connection_obs( &Obs::on_connection_change, &obs ).valid() );
			TS_ASSERT( c2.detach_connection_obs( &Obs::on_connection_change, &obs ) );
			TS_ASSERT( c2.attach_connection_obs( &Obs::on_connection_change, &obs ).valid() );
		} // end scope
		TS_ASSERT_EQUALS( obs.count, 1 ); // ConnectionEvent::DISCONNECT sent on ~Conformation()
		obs.reset();

		// identity observers
		TS_ASSERT( conf.attach_identity_obs( &Obs::on_identity_change, &obs ).valid() );

		conf.replace_residue( 3, *ala_rsd, true );
		TS_ASSERT_EQUALS( obs.count, 1 );
		TS_ASSERT( obs.i_event.tag == signals::IdentityEvent::RESIDUE );

		TS_ASSERT( conf.detach_identity_obs( &Obs::on_identity_change, &obs ) );
		obs.reset();

		// length observers
		TS_ASSERT( conf.attach_length_obs( &Obs::on_length_change, &obs ).valid() );

		conf.append_residue_by_jump( *ala_rsd, conf.size() );
		TS_ASSERT_EQUALS( obs.count, 1 );
		TS_ASSERT( obs.l_event.tag == signals::LengthEvent::RESIDUE_APPEND );

		conf.prepend_polymer_residue_before_seqpos( *ala_rsd, 3, true );
		TS_ASSERT_EQUALS( obs.count, 2 );
		TS_ASSERT( obs.l_event.tag == signals::LengthEvent::RESIDUE_PREPEND );

		conf.delete_polymer_residue( 3 );
		TS_ASSERT_EQUALS( obs.count, 3 );
		TS_ASSERT( obs.l_event.tag == signals::LengthEvent::RESIDUE_DELETE );

		TS_ASSERT( conf.detach_length_obs( &Obs::on_length_change, &obs ) );
		obs.reset();

		// Test event ordering
		TS_ASSERT( conf.attach_length_obs( &Obs::on_length_change, &obs ).valid() );
		TS_ASSERT( conf.attach_xyz_obs( &Obs::on_xyz_change, &obs ).valid() );
		TS_ASSERT( conf.attach_general_obs( &Obs::on_general_change, &obs ).valid() );

		conf.delete_residue_slow( 4 );
		TS_ASSERT_EQUALS( obs.count, 2 ); // One length and one xyz from setup_atom_tree
		TS_ASSERT_LESS_THAN( obs.l_last, obs.x_last ); // Length event should fire before xyz event.

		obs.reset();

		conf.delete_residue_range_slow( 3, 5 );
		TS_ASSERT_EQUALS( obs.count, 2 ); // One length and one xyz from setup_atom_tree

		TS_ASSERT_LESS_THAN( obs.l_last, obs.x_last ); // Length event should fire before xyz event.

		// This is in a different order than adding on purpose.
		TS_ASSERT( conf.detach_xyz_obs( &Obs::on_xyz_change, &obs ) );
		TS_ASSERT( conf.detach_general_obs( &Obs::on_general_change, &obs ) );
		TS_ASSERT( conf.detach_length_obs( &Obs::on_length_change, &obs ) );
		obs.reset();

		// xyz and general observers
		TS_ASSERT( conf.attach_xyz_obs( &Obs::on_xyz_change, &obs ).valid() );
		TS_ASSERT( conf.attach_general_obs( &Obs::on_general_change, &obs ).valid() );

		// the following copy assignment should fire an XYZEvent and one GeneralEvent
		// initiated by the XYZEvent
		conf = Conformation();
		TS_ASSERT_EQUALS( obs.count, 1 );
		TS_ASSERT_EQUALS( obs.g_count, 1 );
		obs.reset();

		TS_ASSERT( conf.detach_xyz_obs( &Obs::on_xyz_change, &obs ) );
		TS_ASSERT( conf.detach_general_obs( &Obs::on_general_change, &obs ) );
		obs.count = 0;
		obs.g_count = 0;
	}

	void test_inter_residue_connection_partner() {
		// A sidechain connection between residues 2 (ASX) & 5 (LYX)
		core::pose::PoseOP pose( core::import_pose::pose_from_file("core/chemical/conn.pdb") );

		utility::vector1< Size > const & r2_conns( pose->residue(2).connections_to_residue( pose->residue(5) ) );
		utility::vector1< Size > const & r5_conns( pose->residue(5).connections_to_residue( pose->residue(2) ) );
		TS_ASSERT_EQUALS( r2_conns.size(), 1 );
		TS_ASSERT_EQUALS( r5_conns.size(), 1 );

		core::id::AtomID partner2( pose->conformation().inter_residue_connection_partner(2, r2_conns[1]) );
		core::id::AtomID partner5( pose->conformation().inter_residue_connection_partner(5, r5_conns[1]) );

		TS_ASSERT_EQUALS( partner2.rsd(), 5 );
		TS_ASSERT_EQUALS( partner2.atomno(), pose->residue(5).atom_index("NZ") );

		TS_ASSERT_EQUALS( partner5.rsd(), 2 );
		TS_ASSERT_EQUALS( partner5.atomno(), pose->residue(2).atom_index("CG") );
	}

	/// @brief Test how ResidueTypeSets in Conformation are cloned
	/// That is, are were doing a semi-shallow copy?
	void test_Conformation_ResidueTypeSet_clone() {
		using namespace core::io::pdb;
		using namespace core::conformation;
		using namespace core::chemical;
		using core::pose::Pose;

		pose::Pose pose( create_test_in_pdb_pose());
		Conformation & orig_conf = pose.conformation();
		// We don't start off with the residue.
		TS_ASSERT( ! pose.conformation().residue_type_set_for_conf( FULL_ATOM_t )->has_name("U11") );

		// Get the PoseRTS, and make a modification.
		PoseResidueTypeSetOP orig_rts = orig_conf.modifiable_residue_type_set_for_conf( FULL_ATOM_t );
		TS_ASSERT( ! orig_rts->has_name("U11") );
		orig_rts->add_base_residue_type( "core/chemical/params/U11.params" );
		TS_ASSERT( orig_rts->has_name("U11") );
		orig_conf.reset_residue_type_set_for_conf( orig_rts );

		// We should now have the ResidueType in the Pose.
		TS_ASSERT( pose.conformation().residue_type_set_for_conf( FULL_ATOM_t )->has_name("U11") );

		// The Conformation should have made a copy of the ResidueTypeSet we put in, so further modifications won't be reflected in the pose.
		TS_ASSERT( ! orig_rts->has_name("U12") );
		TS_ASSERT( ! pose.conformation().residue_type_set_for_conf( FULL_ATOM_t )->has_name("U12") );
		orig_rts->add_base_residue_type( "core/chemical/params/U12.params" );
		TS_ASSERT( orig_rts->has_name("U12") ); // Our RTS has it,
		TS_ASSERT( ! pose.conformation().residue_type_set_for_conf( FULL_ATOM_t )->has_name("U12") ); // but the pose's still doesn't

		// Conformation clones should make copies of the RTS
		ConformationOP clone_conf_op = pose.conformation().clone();
		Conformation & clone_conf = *clone_conf_op; //two-step so we keep the OP around and it isn't garbage collected.
		TS_ASSERT( clone_conf.residue_type_set_for_conf( FULL_ATOM_t )->has_name("U11") );

		// It should be a shallow-ish copy: The residue types shouldn't be cloned.
		ResidueTypeCOP orig_U11( orig_conf.modifiable_residue_type_set_for_conf( FULL_ATOM_t )->name_mapOP( "U11" ) );
		ResidueTypeCOP clone_U11( clone_conf.modifiable_residue_type_set_for_conf( FULL_ATOM_t )->name_mapOP( "U11" ) );
		TS_ASSERT_EQUALS( orig_U11, clone_U11 ); // Yes, comparing the pointers -- we want to make sure these are the same object.

		// The two conformations should be disconnected - alterning one shouldn't result in the other being modified
		orig_rts = orig_conf.modifiable_residue_type_set_for_conf( FULL_ATOM_t );
		orig_rts->add_base_residue_type( "core/chemical/params/U21.params" );
		orig_conf.reset_residue_type_set_for_conf( orig_rts );

		PoseResidueTypeSetOP clone_rts = clone_conf.modifiable_residue_type_set_for_conf( FULL_ATOM_t );
		clone_rts->add_base_residue_type( "core/chemical/params/U31.params" );
		clone_conf.reset_residue_type_set_for_conf( clone_rts );

		TS_ASSERT( orig_conf.residue_type_set_for_conf( FULL_ATOM_t )->has_name("U11") ); // We should still have this.
		TS_ASSERT( orig_conf.residue_type_set_for_conf( FULL_ATOM_t )->has_name("U21") );
		TS_ASSERT( ! orig_conf.residue_type_set_for_conf( FULL_ATOM_t )->has_name("U31") );
		TS_ASSERT( clone_conf.residue_type_set_for_conf( FULL_ATOM_t )->has_name("U11") ); // We should still have this.
		TS_ASSERT( ! clone_conf.residue_type_set_for_conf( FULL_ATOM_t )->has_name("U21") );
		TS_ASSERT( clone_conf.residue_type_set_for_conf( FULL_ATOM_t )->has_name("U31") );

	}

	/// @brief Check that the horrible Conformation::backbone_torsion_angle_atoms() function returns the correct
	/// atoms for a cutpoint variant.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_Conformation_backbone_torsion_angle_atoms_cutpoint() {
		pose::Pose pose( create_test_in_pdb_pose() );
		core::Size const nres(pose.total_residue());
		core::pose::correctly_add_cutpoint_variants( pose, nres, false, 1 );
		pose.conformation().declare_chemical_bond( 1, "N", nres, "C" );

		core::id::TorsionID phi_upper( 1, core::id::BB, 1 );
		core::id::TorsionID psi_upper( 1, core::id::BB, 2 );
		core::id::TorsionID omega_upper( 1, core::id::BB, 3 );
		core::id::TorsionID phi_lower( nres, core::id::BB, 1 );
		core::id::TorsionID psi_lower( nres, core::id::BB, 2 );
		core::id::TorsionID omega_lower( nres, core::id::BB, 3 );

		utility::vector1< core::id::AtomID > lower_phi_atoms(4);
		utility::vector1< core::id::AtomID > lower_psi_atoms(4);
		utility::vector1< core::id::AtomID > lower_omega_atoms(4);
		utility::vector1< core::id::AtomID > upper_phi_atoms(4);
		utility::vector1< core::id::AtomID > upper_psi_atoms(4);
		utility::vector1< core::id::AtomID > upper_omega_atoms(4);

		TS_ASSERT( ! pose.conformation().backbone_torsion_angle_atoms( phi_upper, upper_phi_atoms[1], upper_phi_atoms[2], upper_phi_atoms[3], upper_phi_atoms[4] ) );
		TS_ASSERT( ! pose.conformation().backbone_torsion_angle_atoms( psi_upper, upper_psi_atoms[1], upper_psi_atoms[2], upper_psi_atoms[3], upper_psi_atoms[4] ) );
		TS_ASSERT( ! pose.conformation().backbone_torsion_angle_atoms( omega_upper, upper_omega_atoms[1], upper_omega_atoms[2], upper_omega_atoms[3], upper_omega_atoms[4] ) );
		TS_ASSERT( ! pose.conformation().backbone_torsion_angle_atoms( phi_lower, lower_phi_atoms[1], lower_phi_atoms[2], lower_phi_atoms[3], lower_phi_atoms[4] ) );
		TS_ASSERT( ! pose.conformation().backbone_torsion_angle_atoms( psi_lower, lower_psi_atoms[1], lower_psi_atoms[2], lower_psi_atoms[3], lower_psi_atoms[4] ) );
		TS_ASSERT( ! pose.conformation().backbone_torsion_angle_atoms( omega_lower, lower_omega_atoms[1], lower_omega_atoms[2], lower_omega_atoms[3], lower_omega_atoms[4] ) );

		core::conformation::Residue const &res1(pose.residue(1));
		core::conformation::Residue const &res2(pose.residue(2));
		core::conformation::Residue const &resn_minus_1(pose.residue(nres-1));
		core::conformation::Residue const &resn(pose.residue(nres));

		TS_ASSERT_EQUALS( upper_phi_atoms[1], core::id::AtomID(res1.atom_index("OVU1"), 1) );
		TS_ASSERT_EQUALS( upper_phi_atoms[2], core::id::AtomID(res1.atom_index("N"), 1) );
		TS_ASSERT_EQUALS( upper_phi_atoms[3], core::id::AtomID(res1.atom_index("CA"), 1) );
		TS_ASSERT_EQUALS( upper_phi_atoms[4], core::id::AtomID(res1.atom_index("C"), 1) );

		TS_ASSERT_EQUALS( upper_psi_atoms[1], core::id::AtomID(res1.atom_index("N"), 1) );
		TS_ASSERT_EQUALS( upper_psi_atoms[2], core::id::AtomID(res1.atom_index("CA"), 1) );
		TS_ASSERT_EQUALS( upper_psi_atoms[3], core::id::AtomID(res1.atom_index("C"), 1) );
		TS_ASSERT_EQUALS( upper_psi_atoms[4], core::id::AtomID(res2.atom_index("N"), 2) );

		TS_ASSERT_EQUALS( upper_omega_atoms[1], core::id::AtomID(res1.atom_index("CA"), 1) );
		TS_ASSERT_EQUALS( upper_omega_atoms[2], core::id::AtomID(res1.atom_index("C"), 1) );
		TS_ASSERT_EQUALS( upper_omega_atoms[3], core::id::AtomID(res2.atom_index("N"), 2) );
		TS_ASSERT_EQUALS( upper_omega_atoms[4], core::id::AtomID(res2.atom_index("CA"), 2) );

		TS_ASSERT_EQUALS( lower_phi_atoms[1], core::id::AtomID(resn_minus_1.atom_index("C"), nres - 1) );
		TS_ASSERT_EQUALS( lower_phi_atoms[2], core::id::AtomID(resn.atom_index("N"), nres) );
		TS_ASSERT_EQUALS( lower_phi_atoms[3], core::id::AtomID(resn.atom_index("CA"), nres) );
		TS_ASSERT_EQUALS( lower_phi_atoms[4], core::id::AtomID(resn.atom_index("C"), nres) );

		TS_ASSERT_EQUALS( lower_psi_atoms[1], core::id::AtomID(resn.atom_index("N"), nres) );
		TS_ASSERT_EQUALS( lower_psi_atoms[2], core::id::AtomID(resn.atom_index("CA"), nres) );
		TS_ASSERT_EQUALS( lower_psi_atoms[3], core::id::AtomID(resn.atom_index("C"), nres) );
		TS_ASSERT_EQUALS( lower_psi_atoms[4], core::id::AtomID(resn.atom_index("OVL1"), nres) );

		TS_ASSERT_EQUALS( lower_omega_atoms[1], core::id::AtomID(resn.atom_index("CA"), nres) );
		TS_ASSERT_EQUALS( lower_omega_atoms[2], core::id::AtomID(resn.atom_index("C"), nres) );
		TS_ASSERT_EQUALS( lower_omega_atoms[3], core::id::AtomID(resn.atom_index("OVL1"), nres) );
		TS_ASSERT_EQUALS( lower_omega_atoms[4], core::id::AtomID(resn.atom_index("OVL2"), nres) );

	} //test_Conformation_backbone_torsion_angle_atoms_cutpoint

};
