// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/apps/pilot/jkleman/mpframework_test.cc
/// @brief  testing random small stuff in MPframework
/// @author  JKLeman (julia.koehler1982@gmail.com)

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

	// /// @brief compute last residue number of a chain
	// core::Size chain_end_res( Pose & pose, core::Size chain ) {
	//
	//  // check whether chain exists in pose
	//  if ( ! has_chain( chain, pose ) ) {
	//   TR << "??? Chain should be " << chain << "???" << std::endl;
	//   utility_exit_with_message("Cannot get chain end residue for chain that doesn't exist. Quitting");
	//  }
	//
	//  int int_chain = static_cast< int >( chain );
	//  int end_res(0);
	//  int nres( static_cast< int > ( pose.total_residue() ) );
	//
	//  // go through residues
	//  for ( int i = 1; i <= nres ; ++i ) {
	//
	//   if ( i == nres && pose.chain(i) == int_chain ) {
	//    end_res = nres;
	//   }
	//   else if ( pose.chain( i ) == int_chain && pose.chain( i+1 ) != int_chain ) {
	//    end_res = i;
	//   }
	//  }
	//
	//  return static_cast< core::Size >( end_res );
	//
	// } // chain end residue number

	/////////////////////////////////////////

	/// @brief Create a membrane foldtree with an interface
	/// @details Currently only works for two-body-docking. Both partners can have
	///   multiple chains in any order, the anchoring happens at the TM COM
	///   of each chain
	///
	///       __________________________________________
	///      |________________  _________________       |
	///      |________   iJ   ||________         |      |
	///      |        |       ||        |        |      |
	/// -------  -------  -------  -------  -------  M=root
	///  chain1   chain2   chain3   chain4 ...
	///
	///  iJ = interface jump, will be returned from the function
	///
	// core::Size create_membrane_docking_foldtree_from_partners( Pose const & pose, std::string const partners ) {
	//
	//  using namespace utility;
	//  // pose/util
	//
	//  // split partner string (AB, CDE)
	//  utility::vector1< std::string > partner( utility::string_split( partners, '_' ) );
	//
	//  // initialize partners with chains (will be 1,2 / 3,4,5)
	//  // initialize anchor points within these chains (will be 19,234 / 287,354,528)
	//  // (the anchor points are the chain TM COMs)
	//  utility::vector1< core::Size > chains1;
	//  utility::vector1< core::Size > chains2;
	//  utility::vector1< core::Size > anchors1;
	//  utility::vector1< core::Size > anchors2;
	//  utility::vector1< core::Size > cutpoints1;
	//  utility::vector1< core::Size > cutpoints2;
	//
	//  // go through first partner chainIDs, convert into chain numbers, add to vector
	//  // also get anchor points for these chains
	//  for ( core::Size i = 1; i <= partner[ 1 ].size(); ++i ){
	//
	//   // get chain, add to chains vector and get anchor point
	//   core::Size chain = get_chain_id_from_chain( partner[ 1 ][ i-1 ], pose );
	//   chains1.push_back( chain );
	//   anchors1.push_back( rsd_closest_to_chain_tm_com( pose, chain ) );
	//   cutpoints1.push_back( chain_end_res( pose, chain ) );
	//  }
	//
	//  // go through second partner chainIDs, convert into chain numbers, add to vector
	//  // also get anchor points for these chains
	//  for ( core::Size i = 1; i <= partner[ 2 ].size(); ++i ){
	//
	//   // get chain, add to chains vector and get anchor point
	//   core::Size chain = get_chain_id_from_chain( partner[ 2 ][ i-1 ], pose );
	//   chains2.push_back( chain );
	//   anchors2.push_back( rsd_closest_to_chain_tm_com( pose, chain ) );
	//   cutpoints2.push_back( chain_end_res( pose, chain ) );
	//  }
	//
	//  // create simple foldtree
	//  FoldTree ft = FoldTree();
	//  ft.simple_tree( pose.total_residue() );
	//
	//  // get membrane residue
	//  core::Size memrsd = pose.conformation().membrane_info()->membrane_rsd_num();
	//
	//  // anchor MEM on the first chain of the first partner with the cutpoint
	//  // right before the MEM residues
	//  ft.new_jump( memrsd, anchors1[ 1 ], memrsd - 1 );
	//
	//  // create jumps between the chains in partner1
	//  for ( core::Size i = 2; i <= anchors1.size(); ++i ) {
	//   ft.new_jump( anchors1[ 1 ], anchors1[ i ], cutpoints1[ i-1 ] );
	//  }
	//
	//  // create jumps between the chains in partner1
	//  for ( core::Size i = 2; i <= anchors2.size(); ++i ) {
	//   ft.new_jump( anchors2[ 1 ], anchors2[ i ], cutpoints2[ i-1 ] );
	//  }
	//
	//  // create interface jump between the partners by connecting their 1st chains
	//  // cutpoint is cutpoint of last chain in partner 1
	//  int interface_jump = ft.new_jump( anchors1[ 1 ], anchors2[ 1 ], cutpoints1[ cutpoints1.size() ] );
	//
	//  // reorder and set the foldtree in the pose to the newly created one
	//  ft.reorder( memrsd );
	//  ft.show( TR );
	//  pose.fold_tree( ft );
	//
	//  return static_cast< core::Size >( interface_jump );
	//
	// } // create_membrane_docking_foldtree_from_partners

	/////////////////////////////////////////

	/// @brief Apply Membrane Relax
	void apply( Pose & pose ) {

		core::import_pose::pose_from_pdb( pose, "/Users/julialeman/Documents/julia/git_final/Rosetta/main/source/test/protocols/membrane/3EFF_tr.pdb" );
		AddMembraneMoverOP addmem( new AddMembraneMover( "/Users/julialeman/Documents/julia/git_final/Rosetta/main/source/test/protocols/membrane/3EFF_tr.span" ) );
		addmem->apply( pose );

		Size jumpnum = create_membrane_docking_foldtree_from_partners( pose, "MK_LN" );
		TR << "jump num: " << jumpnum << std::endl;





		//  create_membrane_foldtree_anchor_pose_tmcom( pose );
		//  pose.fold_tree().show( TR );
		//  core::Size pose_anchor( rsd_closest_to_pose_tm_com( pose ) );
		//  TR << "pose anchor " << pose_anchor << std::endl;
		//
		//  for ( int i = 1; i <= static_cast< int >( get_chains( pose ).size() - 1 ); ++i ) {
		//   core::Size chain_anchor ( rsd_closest_to_chain_tm_com( pose, i ) );
		//   TR << "chain anchor " << chain_anchor << std::endl;
		//  }


		////// FOR VISUALIZATION //////

		// initializations
		//  EmbeddingOP embedding( new Embedding() );
		//  core::Vector normal( 0, 0, 1 );
		//
		//  // add to embedding for visualization
		//  embedding->add_span_embedding( tm_com, normal );
		//  embedding->show();
		//
		//  // visualize embeddings
		//  VisualizeEmbeddingMoverOP visemb( new VisualizeEmbeddingMover( embedding ) );
		//  visemb->apply( pose );

	}
	////////////////////////////////////////////////////////////////////////////////
	//
	//  // MPQuickRelaxProtocol
	//
	//  // get job for adding rmsd to scorefile
	//  protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );
	//
	//  // read native
	//  core::pose::Pose native_;
	//  if( option[ OptionKeys::in::file::native ].user() ) {
	//   core::import_pose::pose_from_pdb( native_, option[ OptionKeys::in::file::native ]() );
	//   addmem->apply( native_ );
	//  }
	//
	//  // create scorefunction
	//  ScoreFunctionOP sfxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "mpframework_smooth_fa_2012.wts" );
	//
	//  // starting position: shake up the protein
	//  using namespace core::kinematics;
	//  MoveMapOP mm( new MoveMap() );
	//  mm->set_bb( true );
	//  mm->set_chi( true );
	//  TR << "shaking up the protein for a random starting position" << std::endl;
	//
	//  // get number of residues
	//  Size nres( nres_protein( pose ) );
	//
	//  // set small and shearmover
	//  using namespace protocols::simple_moves;
	//  core::Real kT = 1.0;
	//  SmallMoverOP small( new SmallMover( mm, kT, nres ) );
	//  small->angle_max( 1.0 );
	//  small->apply( pose );
	//
	//  ShearMoverOP shear( new ShearMover( mm, kT, nres ) );
	//  shear->angle_max( 1.0 );
	//  shear->apply( pose );
	//
	//  // create MC object
	//  protocols::moves::MonteCarloOP mc( new protocols::moves::MonteCarlo( pose, *sfxn_, 1.0 ) );
	//
	//  // initialize AtomTreeMinimizer
	//  core::optimization::MinimizerOptions min_opts( "lbfgs_armijo_atol", 0.01, true );
	//  core::optimization::AtomTreeMinimizer atm;
	//
	//  // do this for a certain number of iterations
	//  // since both the packer and especially minimization takes a while, just do one
	//  for ( Size i = 1; i <= 1; ++i ){
	//
	//   // run small and shearmover again
	//   TR << "SmallMover and ShearMover..." << std::endl;
	//   small->apply( pose );
	//   shear->apply( pose );
	//
	//   // packing
	//   TR << "Packing rotamers..." << std::endl;
	//   using namespace core::pack::task;
	//   PackerTaskOP repack = TaskFactory::create_packer_task( pose );
	//   repack->restrict_to_repacking();
	//   core::pack::pack_rotamers( pose, *sfxn_, repack );
	//
	//   // minimize
	//   TR << "Minimizing..." << std::endl;
	//   atm.run( pose, *mm, *sfxn_, min_opts );
	//
	//   // evaluate Boltzmann
	//   mc->boltzmann( pose );
	//   TR << "accepted? " << mc->mc_accepted() << " pose energy: " << pose.energies().total_energy() << std::endl;
	//
	//  } // number of iterations for search
	//
	//  // superimpose poses with native
	//  using namespace protocols::simple_moves;
	//  SuperimposeMoverOP super( new SuperimposeMover( native_, 1, nres, 1, nres, true ) );
	//  super->apply( pose );
	//
	//  // calculate and store the rmsd in the score file
	//  using namespace core::scoring;
	//  job->add_string_real_pair("rms", bb_rmsd( pose, native_ ));

	////////////////////////////////////////////////////////////////////////////////

	// VECTORS
	// set new center and normal
	//  Vector new_center( 0, 5, 10 );
	//  Vector new_normal( 0, 15, 0 );
	//  Vector new_center (0, 0, 0);
	//  Vector new_normal (0, 0, 15);

	// SPANNING TOPOLOGY
	//  std::string spanfile("protocols/membrane/1AFO_AB.span");
	// get topology
	//  SpanningTopologyOP topo( pose.conformation().membrane_info()->spanning_topology() );
	//  pose.conformation().membrane_info()->show();

	// EMBEDDING
	//  // get EmbeddingDef
	//  EmbeddingOP embedding( new Embedding( topo, pose ) );
	//  embedding->show();
	//  embedding->invert();
	//  embedding->show();

	// get embedding object from pose and topology
	////  PoseOP pose1( new Pose( pose ) );
	//  EmbeddingOP embedding( new Embedding( topo, pose ) );
	//  embedding->show();

	//  EmbeddingDefOP embedding( compute_structure_based_embedding(pose) );
	//  embedding->show();

	// FOLDTREE
	// reorder foldtree
	//  pose.fold_tree().show(std::cout);
	//  core::kinematics::FoldTree foldtree = pose.fold_tree();
	//  foldtree.reorder( pose.conformation().membrane_info()->membrane_rsd_num() );
	//  pose.fold_tree( foldtree );
	//  TR << "foldtree reordered" << std::endl;
	//  pose.fold_tree().show(std::cout);

	// VISUALIZE EMBEDDING
	//  VisualizeEmbeddingMoverOP vis_emb( new VisualizeEmbeddingMover( embedding ) );
	//  VisualizeEmbeddingMoverOP vis_emb( new VisualizeEmbeddingMover() );
	//  vis_emb->apply( pose );

	// MOVERS
	//  SetMembranePositionMoverOP rt( new SetMembranePositionMover( new_center, new_normal ) );
	//  rt->apply( pose );

	//  TransformIntoMembraneMoverOP rt( new TransformIntoMembraneMover( new_center, new_normal, spanfile ) );
	//  TransformIntoMembraneMoverOP rt( new TransformIntoMembraneMover() );
	//  rt->apply( pose );

	//  TranslationMoverOP trans = new TranslationMover( translation, jumpnum );
	//  trans->apply( pose );

	//  RotationMoverOP rot = new RotationMover( old_normal, new_normal, old_center, jumpnum );
	//  rot->apply( pose );

	//  TranslationRotationMoverOP rt = new TranslationRotationMover( old_center, old_normal, new_center, new_normal, jumpnum );
	//  rt->apply( pose );

	//  FlipMoverOP flip( new FlipMover(2, axis, 45) );
	//  flip->apply( pose );

	//  RigidBodyRandomizeMoverOP random( new RigidBodyRandomizeMover( pose, 1, partner_downstream, 180, 360, true ) );
	//  random->apply( pose );

	////////////////////////////////////////////////////////////////////////////////

	// VECTORS
	// set new center and normal
	//  Vector new_center( 0, 5, 10 );
	//  Vector new_normal( 0, 15, 0 );
	//  Vector new_center (0, 0, 0);
	//  Vector new_normal (0, 0, 15);

	// SPANNING TOPOLOGY
	//  std::string spanfile("protocols/membrane/1AFO_AB.span");
	// get topology
	//  SpanningTopologyOP topo( pose.conformation().membrane_info()->spanning_topology() );
	//  pose.conformation().membrane_info()->show();

	// EMBEDDING
	//  // get EmbeddingDef
	//  EmbeddingOP embedding( new Embedding( topo, pose ) );
	//  embedding->show();
	//  embedding->invert();
	//  embedding->show();

	// get embedding object from pose and topology
	////  PoseOP pose1( new Pose( pose ) );
	//  EmbeddingOP embedding( new Embedding( topo, pose ) );
	//  embedding->show();

	//  EmbeddingDefOP embedding( compute_structure_based_embedding(pose) );
	//  embedding->show();

	// FOLDTREE
	// reorder foldtree
	//  pose.fold_tree().show(std::cout);
	//  core::kinematics::FoldTree foldtree = pose.fold_tree();
	//  foldtree.reorder( pose.conformation().membrane_info()->membrane_rsd_num() );
	//  pose.fold_tree( foldtree );
	//  TR << "foldtree reordered" << std::endl;
	//  pose.fold_tree().show(std::cout);

	// VISUALIZE EMBEDDING
	//  VisualizeEmbeddingMoverOP vis_emb( new VisualizeEmbeddingMover( embedding ) );
	//  VisualizeEmbeddingMoverOP vis_emb( new VisualizeEmbeddingMover() );
	//  vis_emb->apply( pose );

	// MOVERS
	//  SetMembranePositionMoverOP rt( new SetMembranePositionMover( new_center, new_normal ) );
	//  rt->apply( pose );

	//  TransformIntoMembraneMoverOP rt( new TransformIntoMembraneMover( new_center, new_normal, spanfile ) );
	//  TransformIntoMembraneMoverOP rt( new TransformIntoMembraneMover() );
	//  rt->apply( pose );

	//  TranslationMoverOP trans = new TranslationMover( translation, jumpnum );
	//  trans->apply( pose );

	//  RotationMoverOP rot = new RotationMover( old_normal, new_normal, old_center, jumpnum );
	//  rot->apply( pose );

	//  TranslationRotationMoverOP rt = new TranslationRotationMover( old_center, old_normal, new_center, new_normal, jumpnum );
	//  rt->apply( pose );

	//  FlipMoverOP flip( new FlipMover(2, axis, 45) );
	//  flip->apply( pose );

	//  RigidBodyRandomizeMoverOP random( new RigidBodyRandomizeMover( pose, 1, partner_downstream, 180, 360, true ) );
	//  random->apply( pose );

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

