// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/apps/pilot/jkleman/mpframework_test1.cc
/// @brief  testing random small stuff in MPframework
/// @author  JKLeman (julia.koehler1982@gmail.com)

// Unit Headers
#include <devel/init.hh>

// Package Headers
#include <protocols/jd2/JobDistributor.hh>
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
#include <protocols/membrane/AddMPLigandMover.hh>
#include <protocols/membrane/FlipMover.hh>
#include <protocols/membrane/MembranePositionFromTopologyMover.hh>
#include <protocols/membrane/OptimizeProteinEmbeddingMover.hh>
#include <protocols/membrane/OptimizeMembranePositionMover.hh>
#include <protocols/membrane/TransformIntoMembraneMover.hh>
#include <protocols/membrane/TranslationRotationMover.hh>
#include <protocols/membrane/SpinAroundPartnerMover.hh>
#include <protocols/membrane/SetMembranePositionMover.hh>
#include <protocols/membrane/visualize/VisualizeEmbeddingMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <protocols/scoring/Interface.hh>
#include <core/scoring/methods/RG_Energy_Fast.hh>
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


static THREAD_LOCAL basic::Tracer TR( "apps.pilot.jkleman.mpframework_test1" );

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

	////////////////////////////////////////////////////////////////////////////////

	void move_apart( Pose & pose, int jump, core::Vector axis ) {

		using namespace protocols::rigid;

		RigidBodyTransMoverOP mover( new RigidBodyTransMover( axis, jump ) );
		mover->step_size( 500 );
		mover->apply( pose );
	}

	void move_together( Pose & pose, int jump, core::Vector axis ) {

		using namespace protocols::rigid;
		using namespace core::scoring::methods;
		using namespace protocols::scoring;

		RigidBodyTransMoverOP mover( new RigidBodyTransMover( axis, jump ) );

		// split pose by jump
		Pose pose1, pose2;
		partition_pose_by_jump( pose, jump, pose1, pose2 );

		// compute radii of upstream and downstream partners
		RG_Energy_Fast rg = RG_Energy_Fast();
		core::Real rgyr1 = rg.calculate_rg_score( pose1 );
		core::Real rgyr2 = rg.calculate_rg_score( pose2 );

		// compute first step size
		core::Real step1 = axis.length() - rgyr1 - rgyr2 - 10;
		mover->step_size( step1 );
		mover->apply( pose );

		// small steps to move closer, compute nres in interface
		mover->step_size( 1.0 );
		Interface interface = Interface( jump );
		interface.calculate( pose );

		while ( interface.interface_nres() <= 2 ) {
			mover->apply( pose );
			interface.calculate( pose );
		}

	}


	////////////////////////////////////////////////////////////////////////////////


	/// @brief Apply Membrane Relax
	void apply( Pose & pose ) {

		using namespace protocols::membrane;

		// show foldtree
		pose.fold_tree().show(std::cout);

		//  // before move
		//  pose.dump_pdb("before_addmem.pdb");

		////////////////////////////////////////////////////////////////////////////////

		// AddMembrane
		AddMembraneMoverOP addmem( new AddMembraneMover() );
		addmem->apply( pose );

		AddMPLigandMoverOP add_ligand( new AddMPLigandMover( 118, pose.size()-1 ) );
		add_ligand->apply( pose );

		pose.fold_tree().show( TR );


		//  // TransformIntoMembrane
		//  TransformIntoMembraneMoverOP transform( new TransformIntoMembraneMover() );
		//  transform->apply( pose );

		//  // OptimizeMembranePositionMover
		//  OptimizeProteinEmbeddingMoverOP opt( new OptimizeProteinEmbeddingMover() );
		//  opt->apply( pose );
		//
		//  // compute tilt angle and distance from membrane center
		//  utility::vector1< core::Real > tilt_and_depth( pose_tilt_angle_and_center_distance( pose ) );







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

