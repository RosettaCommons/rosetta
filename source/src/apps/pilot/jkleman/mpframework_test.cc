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
#include <protocols/membrane/TransformIntoMembraneMover.hh>
#include <protocols/membrane/SetMembranePositionMover.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/util.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <protocols/moves/Mover.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh> 
#include <core/kinematics/Jump.hh>
#include <core/types.hh> 
#include <core/kinematics/Jump.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/RT.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <protocols/membrane/geometry/Embedding.hh>
#include <protocols/membrane/geometry/util.hh>
#include <protocols/rigid/RigidBodyMover.hh>

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
using namespace protocols::simple_moves;
using namespace protocols::membrane::geometry;

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

		// show foldtree
		pose.fold_tree().show(std::cout);

		// before move
		pose.dump_pdb("before_addmem.pdb");

////////////////////////////////////////////////////////////////////////////////

		// Add Membrane, appends MEM as jump1
		AddMembraneMoverOP add_memb( new AddMembraneMover() );
		add_memb->apply( pose );

		// before move
		pose.dump_pdb("before.pdb");

		// reorder foldtree
		pose.fold_tree().show(std::cout);
		core::kinematics::FoldTree foldtree = pose.fold_tree();
		foldtree.reorder( pose.conformation().membrane_info()->membrane_rsd_num() );
		pose.fold_tree( foldtree );
		TR << "foldtree reordered" << std::endl;
		pose.fold_tree().show(std::cout);

		// define vectors
		Vector new_center(20, 20, 20);
		Vector new_normal(1, 0, 0);
//		Vector new_center(0, 0, 0);
//		Vector new_normal(0, 0, 1);
				
		Size jumpnum = pose.conformation().membrane_info()->membrane_jump();
		TR << "jump: " << jumpnum << std::endl;
		
		// mover
//		TranslationMoverOP trans = new TranslationMover( translation, jumpnum );
//		trans->apply( pose );

//		RotationMoverOP rot = new RotationMover( old_normal, new_normal, old_center, jumpnum );
//		rot->apply( pose );

//		TranslationRotationMoverOP rt = new TranslationRotationMover( old_center, old_normal, new_center, new_normal, jumpnum );
//		rt->apply( pose );

		TransformIntoMembraneMoverOP rt( new TransformIntoMembraneMover( new_center, new_normal, spanfile_name() ) );
		rt->apply( pose );

		// after move
		pose.dump_pdb("after.pdb");


////////////////////////////////////////////////////

	}
};

typedef utility::pointer::shared_ptr< MPframeworkTestMover > MPframeworkTestMoverOP; 

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
		MPframeworkTestMoverOP load_memb( new MPframeworkTestMover() ); 
		protocols::jd2::JobDistributor::get_instance()->go( load_memb );

		return 0; 

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

