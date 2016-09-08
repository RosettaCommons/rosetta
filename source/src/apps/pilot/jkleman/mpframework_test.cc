// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <basic/options/keys/in.OptionKeys.gen.hh>
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
#include <numeric/numeric.functions.hh>

// C++ headers
#include <iostream>
#include <cstdlib>

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.jkleman.mpframework_test" );

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

		// MPDomainAssembly protocol
		// read in N-terminal PDB

		// read in center domain (can also be refined TMdomain model)

		// read in C-terminal domain model

		// read in missing residues/make sure there aren't any inside the domains
		// missing residues are allowed between the domains up to 15 residues,
		// longer domains should be modeled de novo
		// spit out a warning for this but don't crash!

		// transform TMdomain into membrane and optimize embedding

		// combine poses to single pose

		// set foldtree root in center of TMdomain

		// make sure loops file makes sense: all PDBs and loops files should be
		// renumbered to match the single consecutive pose numbering scheme

		// do loop modeling between domains 1 and 2
		// use KIC without fragments as of yet, it's easier for now
		// protocols/loops/loop_closure/kinematic_closure/KinematicMover

		// do loop modeling between domains 2 and 3

		// relax entire structure

		// optionally: sample flexibility at regions X


		TR << pose.size() << std::endl;




	}

};

////////////////////////////////////////////////////////////////////////////////

/// @brief Remove membrane residue(s)
void remove_membrane_residue( Pose & pose ) {

	// initialize vector of membrane residues
	utility::vector1< core::Size > mem_rsds;

	// make a copy of the pose
	Pose pose_cp( pose );

	// go through pose and find membrane residues
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( pose.residue( i ).name3() == "MEM" ) {
			mem_rsds.push_back( i );
		}
	}

	// remove residues from the end to the start, because otherwise the
	// residue numbers will be incorrect because we changed them by removing
	// a residue
	for ( core::Size i = mem_rsds.size(); i >= 1; --i ) {
		pose_cp.delete_polymer_residue( i );
	}


	// TODO: need to remove the jump from the foldtree, possibly with  /// @brief  Useful for removing a loop modeling jump+cut
	//  void
	//  delete_jump_and_intervening_cutpoint(
	//            int jump_begin,
	//            int jump_end
	//            );
	//
	//   /// @brief  Useful for removing a loop modeling jump+cut
	//  void
	//  delete_jump_and_intervening_cutpoint(
	//            int const jump_number
	//            );

	pose = pose_cp;

} // remove membrane residue

/////////////////////////////////////////

utility::vector1< bool > interface_between_chains( Pose & pose ) {

	// initialize bool of false for each residue
	utility::vector1< bool > intf( nres_protein( pose ), false );

	// iterate over residues
	for ( core::Size r1 = 1; r1 <= nres_protein( pose ); ++r1 ) {

		// get chain of residue 1
		core::Size chain1 = pose.residue( r1 ).chain();

		// get CA coordinates of that residue
		core::Vector coord1 = pose.residue( r1 ).xyz( "CA" );

		// iterate over residues again
		for ( core::Size r2 = r1; r2 <= nres_protein( pose ); ++r2 ) {

			// get chain of residue 2
			core::Size chain2 = pose.residue( r2 ).chain();

			if ( chain1 != chain2 ) {

				// get CA coordinates of that residue
				core::Vector coord2 = pose.residue( r2 ).xyz( "CA" );

				// if distance between two residues is <8A, set interface
				// of that residue to true
				core::Real distance = ( coord1 - coord2 ).length();
				if ( distance < 8.0 ) {
					intf[ r1 ] = true;
					intf[ r2 ] = true;
				}
			} // different chains
		} // residues
	} // residues

	return intf;

} // interface between chains


/////////////////////////////////////////

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

