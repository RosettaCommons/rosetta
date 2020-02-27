// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/thoevenaars/test_jacobi_loop_closure.cc
/// @brief An app for testing the Jacobi loop closure protocol
/// @author teunhoevenaars (teunhoevenaars@gmail.com)

// devel headers
#include <devel/init.hh>

// protocol headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/loops/loop_closure/jacobi/JacobiLoopClosureMover.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/tree/BondedAtom.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/import_pose/import_pose.hh>
#include <core/id/TorsionID.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <core/conformation/Conformation.hh>

// temp
#include <core/conformation/Residue.hh>
#include <core/kinematics/RT.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <core/chemical/ResidueProperties.fwd.hh>
#include <core/chemical/ResidueConnection.hh>

// utility headers
#include <utility/excn/Exceptions.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <utility/options/OptionCollection.hh>
#include <basic/options/option_macros.hh>
#include <src/core/kinematics/tree/BondedAtom.hh>

static basic::Tracer TR("test_jacobi_loop_closure");


/// @brief Indicate which commandline flags are relevant to this application.
void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;



	option.add_relevant( in::file::s );
	option.add_relevant( in::file::l );

}


/// @brief Program entry point.
int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;


		//this won't compile until you fill in brief and default yourself

		devel::init( argc, argv );
		register_options();

		// later in case user can insert both file and loop
		//if ( ( ! option [ in::file::l ].user() ) && ( ! option [ in::file::s ].user() ) ) {
		//    utility_exit_with_message("Please specify either -s or -l to specify the input PDB.");
		//}

		//// INITIALIZATION
		// Read sample pose
		utility::vector1< std::string > filename = {"/Users/teunhoevenaars/Desktop/2sil.clean.pdb"};
		core::pose::PoseOP mypose = (core::import_pose::pose_from_file(filename[ 1 ]));

		// create sample loop
		protocols::loops::Loop loop( 254, 265, 265 );

		// set all angles as in Matlab
		/*for ( core::Size i=loop.start(); i != loop.stop()+2; ++i ) {
		for ( core::Size j=1; j != 4; ++j ) {
		core::id::AtomID const atomnum ( j, i );
		// set torsion angles, indirectly because seems that set_torsion is only processed under conditions
		core::id::DOF_ID const dofphi(atomnum, core::id::PHI);
		if (j != 2) { //phi angle related to first atom is the omega angle of the previous residue
		core::id::TorsionID const torsion_id(i, core::id::BB, j);
		mypose->set_dof(dofphi, 117 * numeric::constants::f::degrees_to_radians);
		} else{
		core::id::TorsionID const torsion_id(i, core::id::BB, j );
		mypose->set_dof(dofphi, 170 * numeric::constants::f::degrees_to_radians);
		}
		// set theta (apparantly in radians)
		core::id::DOF_ID const doftheta(atomnum, core::id::THETA);
		mypose->set_dof(doftheta, 60 * numeric::constants::f::degrees_to_radians);
		// set D
		core::id::DOF_ID const dofd(atomnum, core::id::D);
		mypose->set_dof(dofd, 1.5);
		}
		}*/

		protocols::loops::set_single_loop_fold_tree( *mypose, loop );
		protocols::loops::add_single_cutpoint_variant( *mypose, loop );

		/* for ( core::Size i = loop.start(); i <= loop.stop(); ++i ) {
		mypose->set_phi( i, 50 );
		mypose->set_psi( i, 100 );
		} */
		//mypose->set_phi( loop.stop(), mypose->phi(loop.stop()) + 0 );

		// start PyMol link
		//protocols::moves::PyMOLObserverOP the_observer = protocols::moves::AddPyMOLObserver( *mypose, true, 0 );
		//the_observer->pymol().apply( *mypose);

		//// INITIALIZATION
		// Initialize mover
		protocols::loops::loop_closure::jacobi::JacobiLoopClosureMoverOP jacobi_mover (new protocols::loops::loop_closure::jacobi::JacobiLoopClosureMover( loop ) );

		//// APPLY
		// apply mover
		jacobi_mover->apply(*mypose);
		//protocols::jd2::JobDistributor::get_instance()->go( jacobi_mover );

		TR << "end of loop closure " << std::endl;
		//


	} catch ( utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;
}
