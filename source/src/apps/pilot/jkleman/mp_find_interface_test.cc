// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/apps/pilot/jkleman/mp_find_interface_test.cc
/// @author  JKLeman (julia.koehler1982@gmail.com)

// Unit Headers
#include <devel/init.hh>

// Package Headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/internal_util.hh>

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
#include <protocols/membrane/geometry/util.hh>

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
#include <core/pose/util.hh>
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

// C++ headers
#include <iostream>
#include <cstdlib>

static basic::Tracer TR( "apps.pilot.jkleman.mp_find_interface_test" );

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
class MPFindInterfaceTestMover : public Mover {

public:

	/// @brief Default Constructor
	MPFindInterfaceTestMover() : Mover() {}

	/// @brief Get Mover Name
	std::string get_name() const { return "MPFindInterfaceTestMover"; }

	/// @brief Apply Membrane Relax
	void apply( Pose & pose ) {

		TR << "center: " << pose.conformation().membrane_info()->membrane_center(pose.conformation()).to_string() << std::endl;
		TR << "normal: " << pose.conformation().membrane_info()->membrane_normal(pose.conformation()).to_string()  << std::endl;
		TR << "thickness: " << pose.conformation().membrane_info()->membrane_thickness()  << std::endl;
		TR << "anchor: " << 1 << std::endl;

		// show foldtree
		pose.fold_tree().show(std::cout);

		// before move
		pose.dump_pdb("before_addmem.pdb");

		////////////////////////////////////////////////////////////////////////////////

		// MPFindInterfaceMover

		// get job
		protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );

		// initialize
		SpanningTopologyOP topo_;
		SpanningTopologyOP topo_up_( new SpanningTopology() );
		SpanningTopologyOP topo_down_( new SpanningTopology() );
		utility::vector1< int > jumps_;
		int jump_;
		ScoreFunctionOP sfxn_lowres_;
		std::string partners_;
		core::pose::Pose native_;

		// get partners
		if ( option[ OptionKeys::docking::partners ].user() ) {
			partners_ = option[ OptionKeys::docking::partners ]();
		}

		// read in native and superimpose the first partner of the pose with the native
		if ( option[ OptionKeys::in::file::native ].user() ) {
			core::import_pose::pose_from_file( native_, option[ OptionKeys::in::file::native ]() , core::import_pose::PDB_file);

			// call AddMembraneMover on native for RMSD calculation
			AddMembraneMoverOP addmem( new AddMembraneMover() );
			addmem->apply( native_ );

			// get partner 1: ALL CHAINS MUST BE IN PDB CONSECUTIVELY!!!
			utility::vector1< std::string > partners( utility::string_split( partners_, '_' ) );
			utility::vector1< std::string > partner1( utility::split( partners[1] ) );

			// get residue range for superposition: get start
			Size start(0);
			for ( Size i = 1; i <= pose.size(); ++i ) {
				if ( start == 0 &&
						partner1[1] == utility::to_string( pose.pdb_info()->chain( i ) ) ) {
					start = i;
				}
			}

			// get end
			Size end(0);
			for ( Size j = pose.size(); j >= 1; --j ) {
				if ( end == 0 &&
						partner1[partner1.size()] == utility::to_string( pose.pdb_info()->chain( j ) ) ) {
					end = j;
				}
			}
			TR << "range start: " << start << ", end: " << end << std::endl;

			// superimpose partner 1 with starting pose for easier rmsd calculation
			using namespace protocols::simple_moves;
			SuperimposeMoverOP super( new SuperimposeMover( native_, start, end, start, end, true ) );
			super->apply( pose );

		} // if native given
		pose.dump_pdb("pose_superimposed_to_native.pdb");
		TR << "dumped superimposed pose" << std::endl;

		// get foldtree from partners (setup_foldtree) and movable jump
		setup_foldtree( pose, partners_, jumps_);
		jump_ = jumps_[1];
		TR << "jump_ from foldtree: " << jump_ << std::endl;

		// get anchor point for membrane from jump; I think this residue closest to
		// COM of the upstream partner
		int anchor = pose.fold_tree().upstream_jump_residue( jump_ );
		TR << "anchor point: " << std::endl;

		// Add Membrane, appends MEM as jump1
		AddMembraneMoverOP add_memb( new AddMembraneMover( anchor, 0 ) );
		add_memb->apply( pose );

		// reorder foldtree to have membrane at root
		core::kinematics::FoldTree ft = pose.fold_tree();
		ft.reorder( pose.conformation().membrane_info()->membrane_rsd_num() );
		pose.fold_tree( ft );
		TR << "reordered foltree: " << std::endl;
		pose.fold_tree().show( std::cout );
		TR << "jump: " << jump_ << std::endl;

		// get topology
		topo_ = pose.conformation().membrane_info()->spanning_topology();

		// split_topology_by_jump_noshift
		split_topology_by_jump_noshift( pose, jump_, topo_, topo_up_, topo_down_ );

		// setup lowres scorefunction
		if ( option[ OptionKeys::mp::dock::lowres ].user() ) {
			sfxn_lowres_ = core::scoring::ScoreFunctionFactory::create_score_function( "mpframework_docking_cen_2015.wts" );
		} else if ( option[ OptionKeys::mp::dock::highres ].user() ) {
			sfxn_lowres_ = core::scoring::ScoreFunctionFactory::create_score_function( "mpframework_docking_fa_2015.wts" );
		} else {
			utility_exit_with_message( "Have to specify score function: either use -mp::dock:lowres or -mp::dock:highres flag!" );
		}

		TR << "Sampling protein-protein interface in the membrane..." << std::endl;

		// compute embedding for partners (compute structure-based embedding with split topologies)
		EmbeddingDefOP emb_up( compute_structure_based_embedding( pose, *topo_up_ ) );
		EmbeddingDefOP emb_down( compute_structure_based_embedding( pose, *topo_down_ ) );

		// create random starting position by moving partners apart in the membrane
		TR << "creating random starting position..." << std::endl;
		SpinAroundPartnerMoverOP move_apart( new SpinAroundPartnerMover( jump_, 200 ) );
		move_apart->apply( pose );

		// create MC object
		protocols::moves::MonteCarloOP mc( new protocols::moves::MonteCarlo( pose, *sfxn_lowres_, 1.0 ) );

		// do this for a certain number of iterations
		for ( Size i = 1; i <= 10; ++i ) {

			// SPIN MOVER
			// get a random spin angle between 0 and 360 degrees
			int spin_angle( numeric::random::random_range( 0, 360 ) );
			TR << "=========================SPIN MOVER==========================" << std::endl;
			TR << "random spin angle: " << spin_angle << std::endl;

			// spin downstream partner around spin angle
			emb_up = compute_structure_based_embedding( pose, *topo_up_ );
			emb_down = compute_structure_based_embedding( pose, *topo_down_ );
			RigidBodyDeterministicSpinMoverOP spin( new RigidBodyDeterministicSpinMover(
				jump_, emb_down->normal(), emb_down->center(), spin_angle ) );
			spin->apply( pose );

			// slide into contact
			TR << "=========================SLIDE-INTO-CONTACT==========================" << std::endl;
			DockingSlideIntoContactOP slide( new DockingSlideIntoContact( jump_ ) );
			slide->apply( pose );
			mc->boltzmann( pose );
			TR << "accepted? " << mc->mc_accepted() << " pose energy: " << pose.energies().total_energy() << std::endl;

			// SPIN-AROUND-PARTNER-MOVER!!!
			TR << "=========================SPIN-AROUND-PARTNER MOVER==========================" << std::endl;
			SpinAroundPartnerMoverOP around( new SpinAroundPartnerMover( jump_, 100 ) );
			around->apply( pose );

			// slide into contact
			TR << "=========================SLIDE-INTO-CONTACT==========================" << std::endl;
			slide->apply( pose );
			mc->boltzmann( pose );
			TR << "accepted? " << mc->mc_accepted() << " pose energy: " << pose.energies().total_energy() << std::endl;

			// TILT MOVER
			// get axis perpendicular to downstream partner embedding normal and
			// perpendicular to vector between embedding centers of two partners
			TR << "=========================TILT MOVER==========================" << std::endl;
			TiltMoverOP tilt( new TiltMover( jump_ ) );
			tilt->apply( pose );

			// slide into contact
			TR << "=========================SLIDE-INTO-CONTACT==========================" << std::endl;
			slide->apply( pose );
			mc->boltzmann( pose );
			TR << "accepted? " << mc->mc_accepted() << " pose energy: " << pose.energies().total_energy() << std::endl;

			// FLIP MOVER
			// FlipMover calculates embedding automatically
			TR << "=========================FLIP MOVER==========================" << std::endl;
			protocols::membrane::FlipMoverOP flip( new FlipMover( jump_ ) );
			flip->set_random_membrane_flip_angle();
			flip->apply( pose );

			// slide into contact
			TR << "=========================SLIDE-INTO-CONTACT==========================" << std::endl;
			slide->apply( pose );
			mc->boltzmann( pose );
			TR << "accepted? " << mc->mc_accepted() << " pose energy: " << pose.energies().total_energy() << std::endl;

			// SPIN-AROUND-PARTNER-MOVER!!!
			TR << "=========================SPIN-AROUND-PARTNER MOVER==========================" << std::endl;
			SpinAroundPartnerMoverOP around1( new SpinAroundPartnerMover( jump_, 100 ) );
			around1->apply( pose );

			// slide into contact
			TR << "=========================SLIDE-INTO-CONTACT==========================" << std::endl;
			slide->apply( pose );
			mc->boltzmann( pose );
			TR << "accepted? " << mc->mc_accepted() << " pose energy: " << pose.energies().total_energy() << std::endl;

		} // number of iterations for search

		// calculate and store the rmsd in the score file
		using namespace protocols::docking;
		job->add_string_real_pair("Lrms", calc_Lrmsd( pose, native_, jumps_ ));

	} // apply


	//////////////////////////////////////////////////////////////////////////////

	//   pose.dump_pdb( out + ".pdb" );
	//   pose.dump_scored_pdb( out + ".pdb", *sfxn_lowres_ );

	// FOURTH MOVER: MPDOCKINGMOVER IN LOWRES ONLY
	//   TR << "=========================DOCK MOVER==========================" << std::endl;
	//   protocols::docking::membrane::MPDockingMoverOP mpdock( new protocols::docking::membrane::MPDockingMover( true, false ) );
	//   mpdock->set_cycles_inner( 1 );
	//   mpdock->set_cycles_outer( 1 );
	//   mpdock->apply( pose );
	//   out += "d";
	////   pose.dump_pdb( out + ".pdb" );
	//   pose.dump_scored_pdb( out + ".pdb", *sfxn_lowres_ );

	// add all movers into random mover = Mover container
	//   randmover->add_mover( spin, 0.4 );
	// //  randmover->add_mover( tilt, 0.4 );
	//   randmover->add_mover( flip, 0.2 );

	// slide into contact
	//  DockingSlideIntoContactOP slide( new DockingSlideIntoContact( jump_ ) );
	//  slide->apply( pose );

	//   // evaluate boltzmann criterion and print scores
	//   mc->boltzmann( pose );
	//   TR << "accepted? " << mc->mc_accepted() << " pose energy: " << pose.energies().total_energy() << std::endl;
	//  }

	// visualize embedding
	//  VisualizeEmbeddingMoverOP vis_emb( new VisualizeEmbeddingMover() );
	//  vis_emb->apply( pose );
	//  pose.dump_pdb("after1.pdb");


	////////////////////////////////////////////////////////////////////////////////

	//  MPDockingSetupMoverOP setup( new MPDockingSetupMover() );
	//  setup->apply( pose );

	// TRIALS FOR FINDING INTERFACE USING EXISTING DOCKING MOVERS

	//  DockingInitialPerturbationOP initial( new DockingInitialPerturbation( 1, true ) );
	//  initial->set_randomize1( true );
	//  initial->set_randomize2( false );

	// slide together
	//  Vector axis = membrane_axis( pose, 1 );
	//  DockingSlideIntoContactOP slide( new DockingSlideIntoContact( 1, axis ) );
	//  slide->apply( pose );

	//  ScoreFunctionOP lowres_scorefxn = ScoreFunctionFactory::create_score_function( "mpframework_docking_cen_2015.wts" );
	//
	//  DockingLowResOP lowresdocking( new DockingLowRes( lowres_scorefxn, 1 ) );
	//  lowresdocking->set_trans_magnitude( 100 );
	//  lowresdocking->set_rot_magnitude( 360 );
	//  lowresdocking->apply( pose );

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

typedef utility::pointer::shared_ptr< MPFindInterfaceTestMover > MPFindInterfaceTestMoverOP;

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
		MPFindInterfaceTestMoverOP load_memb( new MPFindInterfaceTestMover() );
		protocols::jd2::JobDistributor::get_instance()->go( load_memb );

		return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
