// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		src/apps/pilot/jkleman/mpframework_test1.cc
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
#include <core/scoring/Energies.hh>

#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/util.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Span.hh>

#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <protocols/membrane/geometry/Embedding.hh>
#include <protocols/membrane/util.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/docking/DockingLowRes.hh>
#include <protocols/docking/DockingInitialPerturbation.hh>
#include <protocols/docking/DockingProtocol.hh>
#include <protocols/docking/DockingLowRes.hh>
#include <protocols/docking/membrane/MPDockingSetupMover.hh>
#include <protocols/docking/membrane/MPDockingMover.hh>
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/FlipMover.hh>
#include <protocols/membrane/MembranePositionFromTopologyMover.hh>
#include <protocols/membrane/OptimizeMembranePositionMover.hh>
#include <protocols/membrane/TransformIntoMembraneMover.hh>
#include <protocols/membrane/TranslationRotationMover.hh>
#include <protocols/membrane/SpinAroundPartnerMover.hh>
#include <protocols/membrane/SetMembranePositionMover.hh>
#include <protocols/membrane/visualize/VisualizeEmbeddingMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <protocols/scoring/Interface.hh>
#include <core/pose/util.hh>

// Utility Headers
#include <core/types.hh>
#include <core/pose/util.hh>
#include <utility/pointer/owning_ptr.hh>
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/conversions.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh> 
#include <protocols/docking/util.hh>

// C++ headers
#include <iostream>
#include <cstdlib> 


static thread_local basic::Tracer TR( "apps.pilot.jkleman.mpframework_test1" );

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
using namespace basic::options;

/// @brief Load Membrane Mover: Load and Initialize a Membrane Protein
/// using Rosetta's new Membrane Framework
class MPframeworkTest1Mover : public Mover {

public: 

	/// @brief Default Constructor
	MPframeworkTest1Mover() : Mover() {}

	/// @brief Get Mover Name
	std::string get_name() const { return "MPframeworkTest1Mover"; }

	/// @brief Apply Membrane Relax
	void apply( Pose & pose ) {

		// show foldtree
		pose.fold_tree().show(std::cout);

		// before move
		pose.dump_pdb("before_addmem.pdb");

////////////////////////////////////////////////////////////////////////////////

		// AddMembrane
		AddMembraneMoverOP addmem( new AddMembraneMover() );
		addmem->apply( pose );

		// membrane center
		core::Vector center = pose.conformation().membrane_info()->membrane_center();
		core::Vector normal = pose.conformation().membrane_info()->membrane_normal();
		TR << "center: " << center.to_string() << std::endl;
		TR << "normal: " << normal.to_string() << std::endl;

		// OptimizeMembrane
		OptimizeMembranePositionMoverOP opt_mem( new OptimizeMembranePositionMover() );
		opt_mem->apply( pose );

		// get center and normal again
		core::Vector center1 = pose.conformation().membrane_info()->membrane_center();
		core::Vector normal1 = pose.conformation().membrane_info()->membrane_normal();
		TR << "center1: " << center1.to_string() << std::endl;
		TR << "normal1: " << normal1.to_string() << std::endl;







//		// reorder foldtree with residue 1 as root of the foldtree
//		core::kinematics::FoldTree foldtree = pose.fold_tree();
//		foldtree.reorder( 1 );
//		pose.fold_tree( foldtree );
//		
//		// show foldtree
//		TR << "foldtree reordered" << std::endl;
//		pose.fold_tree().show(std::cout);
//
//		// MembranePositionFromTopology, structure-based
//		MembranePositionFromTopologyMoverOP initmem( new MembranePositionFromTopologyMover( true ) );
//		initmem->apply( pose );
//
//		// show foldtree
//		TR << "foldtree after MembranePositionFromTopologyMover" << std::endl;
//		pose.fold_tree().show(std::cout);
//		
//		// create scorefxn
//		ScoreFunctionOP sfxn_highres_ = core::scoring::ScoreFunctionFactory::create_score_function( "mpframework_smooth_fa_2012.wts" );
//
//		// MEMBRANE CENTER SEARCH ==============================================
//
//		// starting membrane center
//		TR << "Starting MemInfo center: " << pose.conformation().membrane_info()->membrane_center().to_string() << std::endl;
//		
//		// initialize center get its components
//		core::Vector center( pose.conformation().membrane_info()->membrane_center() );
//		Real xc = center.x();
//		Real yc = center.y();
//		core::Vector new_center;
//		Real new_z( -10 );
//		Real best_z( new_z );
//		
//		// set new starting center
//		new_center.assign( xc, yc, new_z );
//		TR << "starting new_center: " << new_center.to_string() << std::endl;
//		
//		// SetMembraneCenterMover
//		SetMembraneCenterMoverOP set_center( new SetMembraneCenterMover( new_center ) );
//		set_center->apply( pose );
//		TR << "MemInfo center: " << pose.conformation().membrane_info()->membrane_center().to_string() << std::endl;
//		
//		// score the pose
//		Real score_old = sfxn_highres_->score( pose );
//		Real score_new( score_old );
//		Real score_best( score_old );
//		TR << "original score: " << score_old << std::endl;
//		
//		// drag along the z axis and score
//		for ( Size i = 1; i <= 200; ++i ){ // 200
//
//			// drag z through the membrane, 0.1
//			new_z += 0.1;
//			
//			// set new center
//			new_center.assign( xc, yc, new_z );
//			TR << "new_center: " << new_center.to_string() << std::endl;
//			
//			// SetMembraneCenterMover
//			SetMembraneCenterMoverOP set_center1( new SetMembraneCenterMover( new_center ) );
//			set_center1->apply( pose );
//			
//			// rescore pose
//			score_new = sfxn_highres_->score( pose );
//			TR << "score: " << score_new << " original score: " << score_old << std::endl;
//			
//			// save best center
//			if ( score_new < score_best ){
//			
//				best_z = new_z;
//				score_best = score_new;
//				TR << "best center: " << new_center.to_string() << " with score " << score_new << std::endl;
//			}
//			
//			// set old score to this one
//			score_old = score_new;
//		}
//		
//		// Set center to best one that we found
//		new_center.assign( xc, yc, best_z );
//		TR << "Best center: " << new_center.to_string() << std::endl;
//		
//		// SetMembraneCenterMover
//		SetMembraneCenterMoverOP set_center2( new SetMembraneCenterMover( new_center ) );
//		set_center2->apply( pose );
//		
//		// final score from center search
//		score_new = sfxn_highres_->score( pose );
//		TR << "Final score from center search: " << score_new << std::endl;
//		
//		// final membrane center
//		TR << "Final MemInfo center: " << pose.conformation().membrane_info()->membrane_center().to_string() << std::endl;
//
//		// MEMBRANE NORMAL SEARCH ==============================================
//		
//		// starting membrane normal
//		TR << "Starting MemInfo normal: " << pose.conformation().membrane_info()->membrane_normal().to_string() << std::endl;
//		
//		// initialize normal get its components
////		core::Vector normal( pose.conformation().membrane_info()->membrane_normal() );
////		Real xn = normal.x();
////		Real yn = normal.y();
////		Real zn = normal.z();
//
//		// initializations
//		core::Vector new_normal;
//		Real angle( -45 );
//		Real new_x( 0 );
//		Real new_y( 15 * cos( numeric::conversions::radians( angle ) ) );
//		new_z = 15 * sin( numeric::conversions::radians( angle ) );
//		core::Vector best_normal( new_x, new_y, new_y );
//		
//		// set new starting center
//		new_normal.assign( new_x, new_y, new_z );
//		TR << "starting new_normal: " << new_normal.to_string() << std::endl;
//		
//		// SetMembraneNormalMover
//		SetMembraneNormalMoverOP set_normal( new SetMembraneNormalMover( new_normal ) );
//		set_normal->apply( pose );
//		TR << "MemInfo normal: " << pose.conformation().membrane_info()->membrane_normal().to_string() << std::endl;
//		
//		// score the pose
//		score_old = sfxn_highres_->score( pose );
//		score_new = score_old;
//		score_best = score_old;
//		TR << "original score: " << score_old << std::endl;
//		
//		// sample an arch over y-axis, then over x-axis
//		for ( Size j = 1; j <= 4; ++j ){
//		
//			angle = 45;
//			
//			// sample angles in 1 degree steps
//			for ( Size i = 1; i <= 180; ++i ){
//				
//				// increase angle by 1 degree
//				angle += 0.5;
//				
//				// for first iteration, sample arch over x-axis
//				if ( j == 1 ) {
//					new_x = 0;
//					new_y = 15 * cos( numeric::conversions::radians( angle ) );
//					new_z = 15 * sin( numeric::conversions::radians( angle ) );
//				}
//				
//				// for second iteration, sample arch over y-axis
//				else if ( j == 2 ) {
//					new_x = 15 * cos( numeric::conversions::radians( angle ) );
//					new_y = 0;
//					new_z = 15 * sin( numeric::conversions::radians( angle ) );
//				}
//
//				// for third iteration, sample arch over x=y
//				else if ( j == 3 ) {
//					if ( i == 1 ) { angle = -45; }
//					new_x = 15 * sin( numeric::conversions::radians( angle ) );
//					new_y = 15 * sin( numeric::conversions::radians( angle ) );
//					new_z = 15 * cos( numeric::conversions::radians( angle ) );
////					new_x = 15 * sin( numeric::conversions::radians( theta ) ) * cos( numeric::conversions::radians( angle ) );
////					new_y = 15 * sin( numeric::conversions::radians( theta ) ) * sin( numeric::conversions::radians( angle ) );
////					new_z = 15 * cos( numeric::conversions::radians( theta ) );
//				}
//
//				// for fourth iteration, sample arch over x=-y
//				else if ( j == 4 ) {
//					if ( i == 1 ) { angle = -45; }
//					new_x = 15 * sin( numeric::conversions::radians( angle ) );
//					new_y = -15 * sin( numeric::conversions::radians( angle ) );
//					new_z = 15 * cos( numeric::conversions::radians( angle ) );
////					new_x = 15 * sin( numeric::conversions::radians( theta ) ) * cos( numeric::conversions::radians( angle ) );
////					new_y = 15 * sin( numeric::conversions::radians( theta ) ) * sin( numeric::conversions::radians( angle ) );
////					new_z = 15 * cos( numeric::conversions::radians( theta ) );
//				}
//
//				// set new normal
//				new_normal.assign( new_x, new_y, new_z );
//				TR << j << " new_normal: " << new_normal.to_string() << std::endl;
//				
//				// SetMembraneNormalMover
//				SetMembraneNormalMoverOP set_normal1( new SetMembraneNormalMover( new_normal ) );
//				set_normal1->apply( pose );
//				
//				// rescore pose
//				score_new = sfxn_highres_->score( pose );
//				TR << "score: " << score_new << " original score: " << score_old << std::endl;
//				
//				// save best normal
//				if ( score_new < score_best ){
//					
//					best_normal = new_normal;
//					score_best = score_new;
//					TR << "best normal: " << new_normal.to_string() << " with score " << score_new << std::endl;
//				}
//				
//				// set old score to this one
//				score_old = score_new;
//			}
//		}
//			
//		// Set center to best one that we found
//		TR << "Best normal: " << best_normal.to_string() << std::endl;
//		
//		// SetMembraneNormalMover
//		SetMembraneNormalMoverOP set_normal2( new SetMembraneNormalMover( best_normal ) );
//		set_normal2->apply( pose );
//		
//		// final score from center search
//		score_new = sfxn_highres_->score( pose );
//		TR << "Final score from normal search: " << score_new << std::endl;
//		
//		// final membrane normal
//		TR << "Final MemInfo normal: " << pose.conformation().membrane_info()->membrane_normal().to_string() << std::endl;

		

		// iterate over normals and score

		// create MC object
//		protocols::moves::MonteCarloOP mc( new protocols::moves::MonteCarlo( pose, *sfxn_highres_, 1.0 ) );

		
		// score
//		mc->boltzmann( pose );
//		TR << "accepted? " << mc->mc_accepted() << " pose energy: " << pose.energies().total_energy() << std::endl;

		
		
		
		
		
		
		
		
		
		
//		pose.dump_pdb("after1.pdb");


////////////////////////////////////////////////////////////////////////////////

	}
};

typedef utility::pointer::shared_ptr< MPframeworkTest1Mover > MPframeworkTest1MoverOP;

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
		MPframeworkTest1MoverOP load_memb( new MPframeworkTest1Mover() );
		protocols::jd2::JobDistributor::get_instance()->go( load_memb );

		return 0; 

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

