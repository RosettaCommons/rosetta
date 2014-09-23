// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/conformation/Conformation.cxxtest.hh
/// @brief  test suite for Minimizer
/// @author Phil Bradley
/// @author Sergey Lyskov

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/signals/ConnectionEvent.hh>
#include <core/conformation/signals/GeneralEvent.hh>
#include <core/conformation/signals/IdentityEvent.hh>
#include <core/conformation/signals/LengthEvent.hh>
#include <core/conformation/signals/XYZEvent.hh>

// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>

// AUTO-REMOVED #include <core/optimization/AtomTreeMinimizer.hh>
// AUTO-REMOVED #include <core/optimization/MinimizerOptions.hh>

// AUTO-REMOVED #include <core/pose/Pose.hh>

// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>

#include <core/types.hh>

// AUTO-REMOVED #include <basic/Tracer.hh>

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
		core_init();

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
	//pose_from_pdb( pose, "core/conformation/test_in.pdb" );
	pose::Pose pose( create_test_in_pdb_pose());

	kinematics::FoldTree f( pose.total_residue() );
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
	//pose_from_pdb( pose, "core/conformation/test_in.pdb" );
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


};
