// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/sarel/hotspot_stub_constraint_test.cc
/// @brief For running an integration test
/// @author John Karanicolas (johnk@ku.edu)

// Project Headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <devel/init.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <basic/database/open.hh>
#include <core/kinematics/MoveMap.hh>
// AUTO-REMOVED #include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/vector1.hh>

// Unit Headers
//#include <protocols/hotspot_hashing/HotspotHashingConstraints.hh>
#include <protocols/hotspot_hashing/HotspotStubSet.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/MinMover.hh>
// AUTO-REMOVED #include <protocols/docking/DockingProtocol.hh>

#include <core/import_pose/import_pose.hh>
#include <basic/options/keys/OptionKeys.hh>



// C++ headers

static thread_local basic::Tracer TR( "apps.pilot_apps.sarel.hotspot_stub_constraint_test" );

typedef core::Size Size;
typedef core::Real Real;
typedef core::pose::Pose Pose;


void run_test() {

	using namespace core;
	using namespace basic::options;
	using namespace core::scoring;
	using namespace protocols::moves;

	pose::Pose pose;
	core::import_pose::pose_from_pdb( pose, "1a22_prot.pdb" );

	core::scoring::ScoreFunctionOP scorefxn( get_score_function() );
	scorefxn->reset();
	(*scorefxn)(pose);
	TR << "Pre-stub score is: " << pose.energies().total_energy() << std::endl;

	// jk Read the stub set here
	protocols::hotspot_hashing::HotspotStubSetOP hotspot_stub_setOP( new protocols::hotspot_hashing::HotspotStubSet );
	hotspot_stub_setOP->read_data( "jacobHGH1.stubs" );
	protocols::hotspot_hashing::HotspotStubSetCOP saved_hotspot_stub_setCOP( new protocols::hotspot_hashing::HotspotStubSet( *hotspot_stub_setOP ) );

	// jk Setup the PackerTask which will be used in setting up the constraints.
	core::Size const chain_to_redesign = 2;
//	pack::task::PackerTaskOP const hotspot_hash_packer_taskOP =
//		protocols::hotspot_hashing::prepare_hashing_packer_task(pose, chain_to_redesign);

	// Assign a fixed residue (for the constraints)
	// variables below unused ~Labonte
	//core::Size fixed_res(1);
	//if ( chain_to_redesign == 1 ) fixed_res = pose.total_residue();
	//core::id::AtomID fixed_atom_id = core::id::AtomID( pose.residue(fixed_res).atom_index("CA"), fixed_res );

	core::Real const worst_allowed_stub_bonus(-1.);
	bool const apply_self_energies(false);
	core::Real const bump_cutoff(10.);
	bool const apply_ambiguous_constraints(true);

	// Apply strong (short-range) constraints
	hotspot_stub_setOP->add_hotspot_constraints_to_pose( pose, chain_to_redesign, hotspot_stub_setOP, 5., worst_allowed_stub_bonus, apply_self_energies, bump_cutoff, apply_ambiguous_constraints );
	scorefxn->set_weight( core::scoring::backbone_stub_constraint, 1.0 );
	(*scorefxn)(pose);
	TR << "Total score with short-range constraints is: " << pose.energies().total_energies()[ core::scoring::total_score ] << std::endl;
	TR << "Stub score with short-range constraints is: " << pose.energies().total_energies()[ core::scoring::backbone_stub_constraint ] << std::endl;

	// Remove constraints
	hotspot_stub_setOP->remove_all_hotspot_constraints( pose );
//	core::scoring::constraints::ConstraintSetOP empty_constraint_set = new core::scoring::constraints::ConstraintSet;
//	pose.constraint_set( empty_constraint_set );

	(*scorefxn)(pose);
	TR << "Total score after constraint removal is: " << pose.energies().total_energies()[ core::scoring::total_score ] << std::endl;
	TR << "Stub score after constraint removal is: " << pose.energies().total_energies()[ core::scoring::backbone_stub_constraint ] << std::endl;

	// Apply weak (long-range) constraints
	hotspot_stub_setOP = protocols::hotspot_hashing::HotspotStubSetOP( new protocols::hotspot_hashing::HotspotStubSet( *saved_hotspot_stub_setCOP ) );
	hotspot_stub_setOP->add_hotspot_constraints_to_pose( pose, chain_to_redesign, hotspot_stub_setOP, 0.0001, worst_allowed_stub_bonus, apply_self_energies, bump_cutoff, apply_ambiguous_constraints );
	(*scorefxn)(pose);
	TR << "Total score with long-range constraints is: " << pose.energies().total_energies()[ core::scoring::total_score ] << std::endl;
	TR << "Stub score with long-range constraints is: " << pose.energies().total_energies()[ core::scoring::backbone_stub_constraint ] << std::endl;

	// move the chains apart
	core::Size const rb_move_jump = 1; // use the first jump as the one between partners
//	protocols::rigid::RigidBodyTransMover trans_mover( pose, rb_move_jump );
//	protocols::rigid::RigidBodyTransMover trans_mover( pose, rb_move_jump );
//	trans_mover.trans_axis( trans_mover.trans_axis() );
//	trans_mover.step_size(unbound_dist);
//	trans_mover.apply( pose );
	protocols::rigid::RigidBodyPerturbMover rb_pert( rb_move_jump, 15,  3 );
	rb_pert.apply( pose );
	(*scorefxn)(pose);
	scorefxn->show( std::cout, pose );
	TR << "Unbound total score with long-range constraints is: " << pose.energies().total_energies()[ core::scoring::total_score ] << std::endl;
	TR << "Unbound stub score with long-range constraints is: " << pose.energies().total_energies()[ core::scoring::backbone_stub_constraint ] << std::endl;

	// do a minimization, with deriv_check
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
	movemap->set_chi( false );
	movemap->set_bb( false );
	movemap->set_jump( rb_move_jump, true );
	bool const deriv_check(true);
	bool const deriv_check_verbose(true);
	protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover( movemap, scorefxn, "linmin", 0.00001, true, deriv_check, deriv_check_verbose ) );
	//	protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover( movemap, scorefxn, "dfpmin", 0.0001, true, deriv_check, deriv_check_verbose );
	(*scorefxn)(pose);
	min_mover->apply( pose );
	pose.dump_pdb("minimized.pdb");
	pose.energies().clear();
	(*scorefxn)(pose);
	TR << "Minimized total score with long-range constraints is: " << pose.energies().total_energies()[ core::scoring::total_score ] << std::endl;
	TR << "Minimized stub score with long-range constraints is: " << pose.energies().total_energies()[ core::scoring::backbone_stub_constraint ] << std::endl;

	return;
}


int
main( int argc, char * argv [] )
{

	try {

	// setup random numbers and options
	devel::init(argc, argv);

	// run the test
	run_test();

	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

