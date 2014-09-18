// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		src/apps/pilot/jkleman/mpframework_test.cc
/// @brief		testing random small stuff in MPframework
/// @author 	JKLeman (julia.koehler1982@gmail.com)

// Unit Headers
#include <devel/init.hh>

// Package Headers
#include <protocols/membrane/AddMembraneMover.hh> 
#include <protocols/membrane/MembranePositionFromTopologyMover.hh>
#include <protocols/membrane/SetMembranePositionMover.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <protocols/moves/Mover.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh> 
#include <core/kinematics/Jump.hh>
#include <core/types.hh> 

// Utility Headers
#include <utility/pointer/owning_ptr.hh> 
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/conversions.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh> 

// C++ headers
#include <iostream>
#include <cstdlib> 

static thread_local basic::Tracer TR( "apps.pilot.jkleman.mpframework_test" );

using namespace numeric;
using namespace protocols::moves;
using namespace core::kinematics;
using namespace protocols::membrane;

/// @brief Load Membrane Mover: Load and Initialize a Membrane Protein
/// using Rosetta's new Membrane Framework
class MPframeworkTestMover : public Mover {

public: 

	/// @brief Default Constructor
	MPframeworkTestMover() : Mover() {}

	/// @brief Get Mover Name
	std::string get_name() const { return "MPframeworkTestMover"; }

	/// @brief Apply Membrane Relax
	void apply( Pose & pose ) {

		// Add Membrane
		AddMembraneMoverOP add_memb = new AddMembraneMover(); 
		add_memb->apply( pose ); 

		// reorder foldtree
		pose.fold_tree().show(std::cout);
		core::kinematics::FoldTree foldtree = pose.fold_tree();
		foldtree.reorder( pose.conformation().membrane_info()->membrane_rsd_num() );
		pose.fold_tree( foldtree );
		TR << "foldtree reordered" << std::endl;
		pose.fold_tree().show(std::cout);

////////////////////////////////////////////////////////////////////////////////

//		// SetMembranePositionMover
//		// before move
//		pose.dump_pdb("before1.pdb");
//
//		Vector center(20, 20, 0);
//		Vector normal(20, 20, 10);
//
//		//Initialize Membrane
//		SetMembranePositionMoverOP rt_memb = new SetMembranePositionMover( center, normal );
//		rt_memb->apply( pose );
//
//		// after move
//		pose.dump_pdb("after1.pdb");

////////////////////////////////////////////////////////////////////////////////

		// MembranePositionFromTopologyMover
		// before move
		pose.dump_pdb("before2.pdb");

		//Initialize Membrane
		MembranePositionFromTopologyMoverOP init_memb = new MembranePositionFromTopologyMover( true );
		init_memb->apply( pose );

		// normalize normal vector to length 15
		Vector center = pose.conformation().membrane_info()->membrane_center();
		Vector normal = pose.conformation().membrane_info()->membrane_normal();
		normal.normalize( 15 );

		// Update membrane position - shift normal along center
		pose.conformation().update_membrane_position( center, normal );

		// after move
		pose.dump_pdb("after2.pdb");

////////////////////////////////////////////////////////////////////////////////
// TEST ROTATIONS AND TRANSLATIONS OUTSIDE OF MOVERS

		// Show pose membrane
//		pose.conformation().show_membrane(); 
//
//		// Grab current membrane normal from the pose
//		Vector center( pose.conformation().membrane_info()->membrane_center() );
//		Vector current_normal( pose.conformation().membrane_info()->membrane_normal() );
//		current_normal.normalize();
//		
//		TR << "current normal: " << current_normal.to_string() << std::endl;
//
//		TR << "flexible jump: " << pose.conformation().membrane_info()->membrane_jump() << std::endl;
//
//		// reorder foldtree
//		pose.fold_tree().show(std::cout);
//		core::kinematics::FoldTree foldtree = pose.fold_tree();
//		foldtree.reorder( pose.conformation().membrane_info()->membrane_rsd_num() );
//		pose.fold_tree( foldtree );
//		pose.fold_tree().show(std::cout);
//		
//		// Grab the jump from the pose
//		Jump flexible_jump = pose.jump( pose.conformation().membrane_info()->membrane_jump() );
//
//		TR << "flexible jump done" << std::endl;
//
//		// before move
//		pose.dump_pdb("before.pdb");
//
//		// set new normal
//		Vector new_normal (0, 0, 1);
//		TR << "new normal: " << new_normal.to_string() << std::endl;
//
//		// Compute Rotation (delta_rot)
//		Real const theta = numeric::conversions::degrees( angle_of( current_normal, new_normal ) );
//		TR << "theta: " << theta << std::endl;
//
//		xyzVector< Real > const axis = cross( new_normal, current_normal ).normalize();
//		TR << "axis: " << axis.to_string() << std::endl;
//
//		xyzMatrix< Real > rotmatrix = numeric::rotation_matrix_degrees( axis, theta );
//		TR << "rotmatrix: " << std::endl;
//		rotmatrix.show();
//		
//		xyzMatrix< Real > const is_rot = flexible_jump.get_rotation();
//		TR << "is_rot: " << std::endl;
//		is_rot.show();
//
//		flexible_jump.set_rotation( is_rot * rotmatrix );
//		TR << "flexible jump rotation: " << std::endl;
//		flexible_jump.get_rotation().show();
//		
//		// Set jump in the pose
//		pose.set_jump( 1, flexible_jump );
//
//		// after move
//		pose.dump_pdb("after.pdb");

		// Show new pose membrane
//		pose.conformation().show_membrane();

	}
};

typedef utility::pointer::owning_ptr< MPframeworkTestMover > MPframeworkTestMoverOP; 

/// @brief Main method
int
main( int argc, char * argv [] )
{
	try {

		devel::init(argc, argv);

		using namespace protocols::moves;
		using namespace protocols::membrane;

		protocols::jd2::register_options();

		// Create and kick off a new load membrane mover
		MPframeworkTestMoverOP load_memb = new MPframeworkTestMover(); 
		protocols::jd2::JobDistributor::get_instance()->go( load_memb );

		return 0; 

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

