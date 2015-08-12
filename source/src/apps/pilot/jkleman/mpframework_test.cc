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
#include <protocols/jd2/Job.hh>
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
#include <protocols/membrane/TiltMover.hh>
#include <protocols/membrane/SpinAroundPartnerMover.hh>
#include <protocols/membrane/OptimizeMembranePositionMover.hh>
#include <protocols/membrane/TransformIntoMembraneMover.hh>
#include <protocols/membrane/TranslationRotationMover.hh>
#include <protocols/membrane/SetMembranePositionMover.hh>
#include <protocols/membrane/visualize/VisualizeEmbeddingMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <protocols/scoring/Interface.hh>
#include <core/pose/PDBInfo.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>

// Utility Headers
#include <core/types.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
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
#include <utility/string_util.hh>
#include <protocols/toolbox/superimpose.hh>
#include <protocols/docking/metrics.hh>
#include <core/scoring/rms_util.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>

// C++ headers
#include <iostream>
#include <cstdlib> 

static thread_local basic::Tracer TR( "apps.pilot.jkleman.mpframework_test" );

using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace numeric;
using namespace protocols::moves;
using namespace protocols::docking;
using namespace protocols::docking::membrane;
using namespace core::kinematics;
using namespace core::conformation::membrane;
using namespace protocols::rigid;
using namespace protocols::simple_moves;
using namespace protocols::membrane;
using namespace protocols::membrane::geometry;
using namespace protocols::membrane::visualize;
using namespace basic::options;

/// @brief Load Membrane Mover: Load and Initialize a Membrane Protein
/// using Rosetta's new Membrane Framework
class MPframeworkTestMover : public Mover {

public: 

	/// @brief Default Constructor
	MPframeworkTestMover() : Mover() {}

	/// @brief Get Mover Name
	std::string get_name() const { return "MPframeworkTestMover"; }

	/////////////////////////////////////////
	

	/////////////////////////////////////////

	/// @brief Apply Membrane Relax
	void apply( Pose & pose ) {

		TR << "center: " << pose.conformation().membrane_info()->membrane_center().to_string() << std::endl;
		TR << "normal: " << pose.conformation().membrane_info()->membrane_normal().to_string()  << std::endl;
		TR << "thickness: " << pose.conformation().membrane_info()->membrane_thickness()  << std::endl;
		TR << "anchor: " << 1 << std::endl;

		core::import_pose::pose_from_pdb( pose, "/Users/julialeman/Documents/julia/git_final/Rosetta/main/source/test/protocols/membrane/3EFF_tr.pdb" );
		AddMembraneMoverOP addmem( new AddMembraneMover( "/Users/julialeman/Documents/julia/git_final/Rosetta/main/source/test/protocols/membrane/3EFF_tr.span" ) );
		addmem->apply( pose );
		




		////// FOR VISUALIZATION //////

//		// initializations
//		EmbeddingOP embedding( new Embedding() );
//		core::Vector normal( 0, 0, 1 );
//
//		// add to embedding for visualization
//		embedding->add_span_embedding( tmcom1, normal );
//		embedding->add_span_embedding( tmcom2, normal );
//		embedding->add_span_embedding( tmcom3, normal );
//		embedding->add_span_embedding( tmcom4, normal );
//
//		embedding->show();
//		
//		// visualize embeddings
//		VisualizeEmbeddingMoverOP visemb( new VisualizeEmbeddingMover( embedding ) );
//		visemb->apply( pose );


////////////////////////////////////////////////////////////////////////////////
//
//		// MPQuickRelaxProtocol
//
//		// get job for adding rmsd to scorefile
//		protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );
//		
//		// read native
//		core::pose::Pose native_;
//		if( option[ OptionKeys::in::file::native ].user() ) {
//			core::import_pose::pose_from_pdb( native_, option[ OptionKeys::in::file::native ]() );
//			addmem->apply( native_ );
//		}
//		
//		// create scorefunction
//		ScoreFunctionOP sfxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "mpframework_smooth_fa_2012.wts" );
//
//		// starting position: shake up the protein
//		using namespace core::kinematics;
//		MoveMapOP mm( new MoveMap() );
//		mm->set_bb( true );
//		mm->set_chi( true );
//		TR << "shaking up the protein for a random starting position" << std::endl;
//		
//		// get number of residues
//		Size nres( nres_protein( pose ) );
//
//		// set small and shearmover
//		using namespace protocols::simple_moves;
//		core::Real kT = 1.0;
//		SmallMoverOP small( new SmallMover( mm, kT, nres ) );
//		small->angle_max( 1.0 );
//		small->apply( pose );
//		
//		ShearMoverOP shear( new ShearMover( mm, kT, nres ) );
//		shear->angle_max( 1.0 );
//		shear->apply( pose );
//
//		// create MC object
//		protocols::moves::MonteCarloOP mc( new protocols::moves::MonteCarlo( pose, *sfxn_, 1.0 ) );
//		
//		// initialize AtomTreeMinimizer
//		core::optimization::MinimizerOptions min_opts( "lbfgs_armijo_atol", 0.01, true );
//		core::optimization::AtomTreeMinimizer atm;
//		
//		// do this for a certain number of iterations
//		// since both the packer and especially minimization takes a while, just do one
//		for ( Size i = 1; i <= 1; ++i ){
//		
//			// run small and shearmover again
//			TR << "SmallMover and ShearMover..." << std::endl;
//			small->apply( pose );
//			shear->apply( pose );
//		
//			// packing
//			TR << "Packing rotamers..." << std::endl;
//			using namespace core::pack::task;
//			PackerTaskOP repack = TaskFactory::create_packer_task( pose );
//			repack->restrict_to_repacking();
//			core::pack::pack_rotamers( pose, *sfxn_, repack );
//			
//			// minimize
//			TR << "Minimizing..." << std::endl;
//			atm.run( pose, *mm, *sfxn_, min_opts );
//			
//			// evaluate Boltzmann
//			mc->boltzmann( pose );
//			TR << "accepted? " << mc->mc_accepted() << " pose energy: " << pose.energies().total_energy() << std::endl;
//
//		} // number of iterations for search
//		
//		// superimpose poses with native
//		using namespace protocols::simple_moves;
//		SuperimposeMoverOP super( new SuperimposeMover( native_, 1, nres, 1, nres, true ) );
//		super->apply( pose );
//		
//		// calculate and store the rmsd in the score file
//		using namespace core::scoring;
//		job->add_string_real_pair("rms", bb_rmsd( pose, native_ ));

	} // apply


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

