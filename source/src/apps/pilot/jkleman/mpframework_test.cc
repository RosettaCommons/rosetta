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
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>

#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/RT.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/util.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Span.hh>

#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <protocols/membrane/geometry/Embedding.hh>
#include <protocols/membrane/geometry/util.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/docking/DockingLowRes.hh>
#include <protocols/docking/DockingInitialPerturbation.hh>
#include <protocols/docking/membrane/MPDockingSetupMover.hh>
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/FlipMover.hh>
#include <protocols/membrane/TransformIntoMembraneMover.hh>
#include <protocols/membrane/SetMembranePositionMover.hh>
#include <protocols/membrane/visualize/VisualizeEmbeddingMover.hh>

// Utility Headers
#include <core/types.hh>
#include <core/pose/util.hh>
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



#include <core/conformation/membrane/types.hh>

static thread_local basic::Tracer TR( "apps.pilot.jkleman.mpframework_test" );

using namespace core;
using namespace core::scoring;
using namespace numeric;
using namespace protocols::moves;
using namespace protocols::docking;
using namespace protocols::docking::membrane;
using namespace core::kinematics;
using namespace core::conformation::membrane;
using namespace protocols::rigid;
using namespace protocols::membrane;
using namespace protocols::simple_moves;
using namespace protocols::membrane::geometry;
using namespace protocols::membrane::visualize;

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

		TR << "center: " << mem_center.to_string() << std::endl;
		TR << "normal: " << mem_normal.to_string() << std::endl;
		TR << "thickness: " << mem_thickness << std::endl;
		TR << "anchor: " << mem_anchor << std::endl;

		// show foldtree
		pose.fold_tree().show(std::cout);

		// before move
		pose.dump_pdb("before_addmem.pdb");

////////////////////////////////////////////////////////////////////////////////

		// Add Membrane, appends MEM as jump1
//		AddMembraneMoverOP add_memb( new AddMembraneMover() );
//		add_memb->apply( pose );

		MPDockingSetupMoverOP setup( new MPDockingSetupMover() );
		setup->apply( pose );

		
		// visualize embedding
//		VisualizeEmbeddingMoverOP vis_emb( new VisualizeEmbeddingMover() );
//		vis_emb->apply( pose );
		pose.dump_pdb("after1.pdb");


////////////////////////////////////////////////////////////////////////////////

// TRIALS FOR FINDING INTERFACE USING EXISTING DOCKING MOVERS

//		DockingInitialPerturbationOP initial( new DockingInitialPerturbation( 1, true ) );
//		initial->set_randomize1( true );
//		initial->set_randomize2( false );
		
// slide together
//		Vector axis = membrane_axis( pose, 1 );
//		DockingSlideIntoContactOP slide( new DockingSlideIntoContact( 1, axis ) );
//		slide->apply( pose );

//		ScoreFunctionOP lowres_scorefxn = ScoreFunctionFactory::create_score_function( "mpframework_docking_cen_2015.wts" );
//		
//		DockingLowResOP lowresdocking( new DockingLowRes( lowres_scorefxn, 1 ) );
//		lowresdocking->set_trans_magnitude( 100 );
//		lowresdocking->set_rot_magnitude( 360 );
//		lowresdocking->apply( pose );

////////////////////////////////////////////////////////////////////////////////

		// VECTORS
		// set new center and normal
//		Vector new_center( 0, 5, 10 );
//		Vector new_normal( 0, 15, 0 );
//		Vector new_center (0, 0, 0);
//		Vector new_normal (0, 0, 15);

		// SPANNING TOPOLOGY
//		std::string spanfile("protocols/membrane/1AFO_AB.span");
		// get topology
//		SpanningTopologyOP topo( pose.conformation().membrane_info()->spanning_topology() );
//		pose.conformation().membrane_info()->show();

		// EMBEDDING
//		// get EmbeddingDef
//		EmbeddingOP embedding( new Embedding( topo, pose ) );
//		embedding->show();
//		embedding->invert();
//		embedding->show();

		// get embedding object from pose and topology
////		PoseOP pose1( new Pose( pose ) );
//		EmbeddingOP embedding( new Embedding( topo, pose ) );
//		embedding->show();

//		EmbeddingDefOP embedding( compute_structure_based_embedding(pose) );
//		embedding->show();

		// FOLDTREE
		// reorder foldtree
//		pose.fold_tree().show(std::cout);
//		core::kinematics::FoldTree foldtree = pose.fold_tree();
//		foldtree.reorder( pose.conformation().membrane_info()->membrane_rsd_num() );
//		pose.fold_tree( foldtree );
//		TR << "foldtree reordered" << std::endl;
//		pose.fold_tree().show(std::cout);
		
		// VISUALIZE EMBEDDING
//		VisualizeEmbeddingMoverOP vis_emb( new VisualizeEmbeddingMover( embedding ) );
//		VisualizeEmbeddingMoverOP vis_emb( new VisualizeEmbeddingMover() );
//		vis_emb->apply( pose );
		
		// MOVERS
//		SetMembranePositionMoverOP rt( new SetMembranePositionMover( new_center, new_normal ) );
//		rt->apply( pose );

//		TransformIntoMembraneMoverOP rt( new TransformIntoMembraneMover( new_center, new_normal, spanfile ) );
//		TransformIntoMembraneMoverOP rt( new TransformIntoMembraneMover() );
//		rt->apply( pose );
		
//		TranslationMoverOP trans = new TranslationMover( translation, jumpnum );
//		trans->apply( pose );

//		RotationMoverOP rot = new RotationMover( old_normal, new_normal, old_center, jumpnum );
//		rot->apply( pose );

//		TranslationRotationMoverOP rt = new TranslationRotationMover( old_center, old_normal, new_center, new_normal, jumpnum );
//		rt->apply( pose );

//		FlipMoverOP flip( new FlipMover(2, axis, 45) );
//		flip->apply( pose );

//		RigidBodyRandomizeMoverOP random( new RigidBodyRandomizeMover( pose, 1, partner_downstream, 180, 360, true ) );
//		random->apply( pose );

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

