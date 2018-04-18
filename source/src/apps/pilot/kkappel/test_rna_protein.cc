// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <iostream>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <devel/init.hh>
#include <utility/pointer/owning_ptr.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <numeric/random/random.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMolMover.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/pose/util.hh>
#include <protocols/jd2/JobDistributor.hh>

void test_protein_rna() {
	utility::vector1< std::string > input_files = basic::options::option[ basic::options::OptionKeys::in::file::s ]();
	core::pose::PoseOP pose;
	if ( input_files.size() > 0 ) {
		pose = core::import_pose::pose_from_pdb( input_files[ 1 ] );
	} else {
		std::cout << "No input structure provided!" << std::endl;
		return;
	}

	// Set up the pymol viewer
	protocols::moves::PyMolObserverOP observer = protocols::moves::AddPyMolObserver( *pose, true, 0 );

	core::optimization::MinimizerOptions min_opts( "lbfgs_armijo_atol", 0.01, true );
	core::optimization::AtomTreeMinimizer atm;
	core::scoring::ScoreFunctionOP sfxn;
	sfxn = core::scoring::get_score_function();
	core::kinematics::MoveMap mm;
	mm.set_bb( true );
	mm.set_chi( true );
	pose->fold_tree().show(std::cout);
	atm.run( *pose, mm, *sfxn, min_opts );
}

int main( int argc, char ** argv ) {

	try {
		devel::init( argc, argv );
		test_protein_rna();
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;



	//// utility::vector1< std::string > filenames = basic::options::option[ basic::options::OptionKeys::in::file::s ]();
	//// if ( filenames.size() > 0 ) {
	//// std::cout << "You entered: " << filenames[ 1 ] << " as the PDB file to be read" << std::endl;
	//// } else {
	////  std::cout << "You didn’t provide a PDB file with the -in::file::s option" << std::endl;
	////  return 1;
	//// }
	//// core::pose::PoseOP mypose = core::import_pose::pose_from_pdb( filenames[ 1 ] );
	//// mypose->fold_tree().show(std::cout);
	//        core::scoring::ScoreFunctionOP sfxn;
	//        sfxn = core::scoring::get_score_function();
	// //protocols::bootcamp::BootCampMoverOP mover( new protocols::bootcamp::BootCampMover );
	// protocols::bootcamp::BootCampMoverOP mover( new protocols::bootcamp::BootCampMover(sfxn) );
	//// mover->set_score_function( sfxn );
	// protocols::jd2::JobDistributor::get_instance()->go(mover);
	//
	//
	//// core::scoring::ScoreFunctionOP sfxn;
	//// sfxn = core::scoring::get_score_function();
	//// core::Real score = ( *sfxn )( *mypose );
	//// //Print out the score for the pose
	//// std::cout << score << std::endl;
	////
	////// this stuff should go in the bootcamp mover
	//// //Set up the monte carlo mover
	//// protocols::moves::MonteCarloOP mc = protocols::moves::MonteCarloOP( new protocols::moves::MonteCarlo( *mypose, *sfxn, 1 ) );
	////
	//// // Set up minimization
	//// core::kinematics::MoveMap mm;
	//// mm.set_bb( true );
	//// mm.set_chi( true );
	//// core::optimization::MinimizerOptions min_opts( "lbfgs_armijo_atol", 0.01, true );
	//// core::optimization::AtomTreeMinimizer atm;
	//// // Make a copy of the pose to speed things up for viewing with pymol
	//// core::pose::Pose copy_pose;
	////
	//// // Setup pymol viewing
	//// protocols::moves::PyMolObserverOP observer = protocols::moves::AddPyMolObserver( *mypose, true, 0 );
	////
	//// core::Real avg_energy( 0 );
	////
	//// // Make a fold tree from the secondary structure of the pose
	//// core::scoring::dssp::Dssp dssp( *mypose );
	//// dssp.insert_ss_into_pose( *mypose );
	//// core::kinematics::FoldTree foldtree_ss = protocols::bootcamp::fold_tree_from_ss( *mypose );
	//// mypose->fold_tree( foldtree_ss );
	//// std::cout << "Successfully reset the fold tree from the ss!" << std::endl;
	//// core::pose::correctly_add_cutpoint_variants( *mypose );
	//// sfxn->set_weight( core::scoring::linear_chainbreak, 1);
	////
	//// // The main loop
	//// for ( int i = 1; i<=300; ++i ) {
	////  //Copy the pose for viewing
	////  core::Size total_residues = mypose->total_residue();
	////  // Pick a random residue in the pose
	////  core::Size randres = static_cast< core::Size > ( numeric::random::uniform() * total_residues + 1);
	////  core::Real pert1 = numeric::random::gaussian();
	////  core::Real pert2 = numeric::random::gaussian();
	////  // Make perturbations to the pose
	////  core::Real orig_phi = mypose->phi( randres );
	////  core::Real orig_psi = mypose->psi( randres );
	////  mypose->set_phi( randres, orig_phi + pert1 );
	////  mypose->set_psi( randres, orig_psi + pert2 );
	////  // Run the packer
	////  core::pack::task::PackerTaskOP repack_task = core::pack::task::TaskFactory::create_packer_task( *mypose );
	////  repack_task->restrict_to_repacking();
	////  core::pack::pack_rotamers( *mypose, *sfxn, repack_task );
	////  copy_pose = *mypose;
	////  // Do the minimization
	////  atm.run( copy_pose, mm, *sfxn, min_opts );
	////  *mypose = copy_pose;
	////  // Check the metropolis criterion
	////  mc->boltzmann( *mypose ); // This will modify the pose if rejected
	////  avg_energy += mypose->energies().total_energy();
	////  //std::cout << "Average Energy: " << avg_energy / i << std::endl;
	////  //std::cout << "Score: " << (*sfxn)(*mypose) << std::endl;
	////  if ( i%100 == 0 ) {
	////   mc->show_counters();
	////   std::cout << "Average Energy: " << avg_energy / i << std::endl;
	////  }
	//// }
	//// return 0;
}
