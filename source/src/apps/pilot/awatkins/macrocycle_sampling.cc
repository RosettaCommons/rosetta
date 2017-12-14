// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   oop_conformation.cc
/// @brief  Creates an OOP dimer and toys around with its dihedrals.
/// @author Watkins


// Package headers
#include <devel/init.hh>

// Project headers
#include <core/types.hh>
#include <core/io/carbohydrates/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>
//#include <core/pose/PDBInfo.hh>
//#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
//#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
//#include <core/id/TorsionID.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
//#include <core/pack/task/PackerTask.hh>
//#include <core/pack/task/TaskFactory.hh>

#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/RandomTorsionMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/docking/DockingInitialPerturbation.hh>
#include <protocols/docking/util.hh>
#include <protocols/relax/AtomCoordinateCstMover.hh>

#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopBuilder.hh>
#include <protocols/loop_modeling/LoopProtocol.hh>
#include <protocols/loop_modeling/samplers/LegacyKicSampler.hh>
#include <protocols/loop_modeling/refiners/MinimizationRefiner.hh>
#include <protocols/loop_modeling/refiners/RotamerTrialsRefiner.hh>
#include <protocols/loop_modeling/refiners/RepackingRefiner.hh>
#include <protocols/kinematic_closure/KicMover.hh>
#include <protocols/kinematic_closure/perturbers/RamaPerturber.hh>
#include <protocols/kinematic_closure/perturbers/FragmentPerturber.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/util.hh>



// Utility headers
#include <utility/vector1.hh>
#include <utility/string_util.hh>

// Numeric headers
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>

// C++ headers
#include <iostream>
//#include <algorithm>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <protocols/ncbb/oop/OopCreatorMover.hh>
#include <protocols/simple_moves/oop/OopRandomSmallMover.hh>

using namespace basic::options;
using namespace utility;
using namespace core;
using namespace core::chemical;
using namespace kinematics;
using namespace scoring;
using namespace import_pose;
using namespace pose;
using namespace protocols;
using namespace simple_moves;
using namespace numeric::random;


utility::vector1< Real >
load_bb_torsions_into_vector(
	core::pose::Pose & pose
) {
	utility::vector1< Real > bb_vec;
	for ( core::Size ii = 1; ii <= pose.size(); ++ii ) {
		for ( core::Size jj = 1; jj <= pose.residue( ii ).mainchain_torsions().size(); ++jj ) {
			bb_vec.push_back( pose.residue( ii ).mainchain_torsions()[ jj ] );
		}
	}
	return bb_vec;
}

Real
compute_deficit(
	utility::vector1< Real > v1,
	utility::vector1< Real > v2
) {
	runtime_assert( v1.size() == v2.size() );

	Real deficit = 0;
	for ( core::Size ii = 1; ii <= v1.size(); ++ii ) {
		deficit += v2[ ii ] - v1[ ii ];
	}
	return deficit;
}

int
main( int argc, char *argv[] )
{
	try {
		// initialize core
		devel::init( argc, argv );


		scoring::ScoreFunctionOP score_fxn = scoring::get_score_function();

		score_fxn->set_weight( core::scoring::atom_pair_constraint, 1 );
		score_fxn->set_weight( core::scoring::angle_constraint, 1.0 );
		score_fxn->set_weight( core::scoring::dihedral_constraint, 1 );//10.0 );

		Pose pose;
		std::string filename = option[OptionKeys::in::file::s].value()[1];
		import_pose::pose_from_file( pose, filename , core::import_pose::PDB_file);
		filename.erase( filename.find( ".pdb" ) );


		// okay, we want to keep the macrocycle closed
		// this sampling strategy works on anything but is essential for macrocycles
		// 1. pick random residue
		// 2. perturb dihedral
		// 3. that perturbation value is "debt"
		// 4. perturb rest of dihedrals in other direction to resolve debt

		//kinematics::MoveMapOP mm( new kinematics::MoveMap );
		// Because of gripes from the AtomTree, let's avoid the question
		// of the un-settable BB torsions on resi 1 and npos
		//mm->set_bb( true );
		/*for ( core::Size  ii = 2; ii <= pose.size() - 1; ++ii ) {
		mm->set_bb( ii, true );
		}
		RandomTorsionMoverOP tor_mover( new RandomTorsionMover( mm, 90, 1 ) );
		RandomTorsionMoverOP small_tor_mover( new RandomTorsionMover( mm, 5, 10 ) );

		utility::vector1< Real > init_torsions;
		utility::vector1< Real > final_torsions;

		moves::MonteCarloOP pert_mc( new moves::MonteCarlo( pose, *score_fxn, 3 ) );

		Real best_score = 100000;
		for ( core::Size ii = 1; ii <= 1; ++ii ) {
		init_torsions = load_bb_torsions_into_vector( pose );
		tor_mover->apply( pose );
		//final_torsions = load_bb_torsions_into_vector( pose );
		Real deficit;
		core::Size jj = 0;
		do {
		small_tor_mover->apply( pose );
		final_torsions = load_bb_torsions_into_vector( pose );
		deficit = compute_deficit( init_torsions, final_torsions );
		pose.dump_pdb( filename + "_" + utility::to_string( ++jj ) + ".pdb" );
		std::cout << "deficit: " <<  deficit << std::endl;
		} while ( deficit > 5 || deficit < -5 );

		//std::cout << "pre mc->boltzmann" << std::endl;
		//pert_mc->show_state();
		if ( pert_mc->boltzmann( pose ) ) {
		Real test_score( pose.energies().total_energies().dot( score_fxn->weights() ) );
		if( test_score <= best_score ){
		best_score = test_score;
		pert_mc->reset( pose );
		}
		}
		//std::cout << "post mc->boltzmann" << std::endl;
		//pert_mc->show_state();

		//std::cout << "Score " << ( *score_fxn )( pose ) << std::endl;
		if ( ! ( ii % 100 ) ) pose.dump_pdb( filename + "_" + utility::to_string( ii ) + ".pdb" );
		}
		pert_mc->recover_low( pose );

		*/

		auto cutpoint = static_cast<core::Size>( ( 2 + pose.size()-1 ) / 2 );
		//protocols::loops::Loop loop( loop_start, loop_end, cutpoint );
		//TR << "i.e. " << loop << std::endl;
		loops::LoopsOP loops( new protocols::loops::Loops );
		loops->push_back( 2, pose.size()-1, cutpoint );//loop );
		std::cout << (*loops) << std::endl;

		// Add constraints to native except loop

		pose::PoseOP native_pose( new core::pose::Pose( pose ) );
		native_pose->conformation().delete_residue_range_slow(2, pose.size()-1);
		native_pose->conformation().detect_disulfides();

		relax::AtomCoordinateCstMover coord_cst;
		coord_cst.set_refstruct( native_pose );
		coord_cst.apply( pose );

		// fragment initialization (not necessary in this case)
		utility::vector1< core::fragment::FragSetOP > frag_libs;

		//setup of looprelax_mover
		using namespace std;
		using protocols::loop_modeling::LoopMoverOP;
		using protocols::loop_modeling::LoopProtocol;
		using protocols::loop_modeling::LoopProtocolOP;
		using protocols::loop_modeling::refiners::RepackingRefiner;
		using protocols::loop_modeling::refiners::RotamerTrialsRefiner;
		using protocols::loop_modeling::refiners::MinimizationRefiner;
		using protocols::kinematic_closure::KicMover;
		using protocols::kinematic_closure::KicMoverOP;
		using protocols::kinematic_closure::perturbers::RamaPerturber;
		using protocols::kinematic_closure::perturbers::FragmentPerturber;

		std::cout << "Beginning full-atom KIC sampling..." << endl;

		core::Size sfxn_cycles = option[ OptionKeys::loops::refine_outer_cycles ]();
		core::Size temp_cycles = 10 * loops->loop_size();
		core::Size repack_period = 20;

		if ( option[OptionKeys::loops::max_inner_cycles].user() ) {
			core::Size max_cycles = option[OptionKeys::loops::max_inner_cycles]();
			temp_cycles = std::max(temp_cycles, max_cycles);
		}
		if ( option[OptionKeys::loops::fast] ) {
			sfxn_cycles = 3;
			temp_cycles = loops->loop_size();
		}
		if ( option[OptionKeys::run::test_cycles] ) {
			sfxn_cycles = 3;
			temp_cycles = 3;
		}
		if ( option[OptionKeys::loops::repack_period].user() ) {
			repack_period = option[OptionKeys::loops::repack_period]();
		}

		LoopProtocolOP protocol( new LoopProtocol );
		KicMoverOP kic_mover( new KicMover );
		kic_mover->clear_perturbers();
		kic_mover->add_perturber(kinematic_closure::perturbers::PerturberOP( new RamaPerturber ));//to emulate legacy KIC behavior

		if ( option[OptionKeys::loops::ramp_rama].user() ) {
			protocol->set_rama_term_ramping(true);
		}
		if ( option[OptionKeys::loops::ramp_fa_rep].user() ) {
			protocol->set_repulsive_term_ramping(true);
		}

		protocol->set_loops(*loops);
		protocol->set_score_function(score_fxn);
		protocol->set_sfxn_cycles(sfxn_cycles);
		protocol->set_temp_cycles(temp_cycles);
		protocol->set_mover_cycles(2);
		protocol->add_mover(kic_mover);
		protocol->add_mover(LoopMoverOP( new RepackingRefiner(repack_period) ));
		protocol->add_mover(LoopMoverOP( new RotamerTrialsRefiner ));
		protocol->add_mover(LoopMoverOP( new MinimizationRefiner ));
		protocol->apply(pose);


		pose.dump_pdb( filename+"_final.pdb" );

	} catch (utility::excn::Exception const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
