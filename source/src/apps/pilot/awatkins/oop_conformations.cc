// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/docking/DockingInitialPerturbation.hh>
#include <protocols/docking/util.hh>

// Utility headers
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/random/random.hh>

// C++ headers
#include <iostream>
//#include <algorithm>

#include <protocols/ncbb/oop/OopCreatorMover.hh>
#include <protocols/simple_moves/oop/OopRandomSmallMover.hh>

using namespace basic::options;
using namespace basic::options::OptionKeys;

namespace oop_conformations {
	// pert options
	RealOptionKey const step ( "oop_conformations::step" );
	RealOptionKey const dihedral_start ( "oop_conformations::dihedral_start" );
	RealOptionKey const dihedral_end ( "oop_conformations::dihedral_end" );
	BooleanOptionKey const oop_optimize( "oop_conformations::oop_optimize" );
	// BooleanOptionKey const correct_hbs_dihedrals ( "hbs_creator::correct_hbs_dihedrals" ); to be implemented if possible
	
}

void
idealize( core::pose::Pose & pose, core::scoring::ScoreFunctionOP scorefxn )
{
	using namespace core;
	using namespace protocols;
	
	moves::SequenceMoverOP pert_sequence( new moves::SequenceMover() );
	moves::MonteCarloOP pert_mc( new moves::MonteCarlo( pose, *scorefxn, 0.2 ) );
	
	kinematics::MoveMapOP pert_pep_mm( new kinematics::MoveMap() );
	
	utility::vector1< Size > oop_pre_positions;

	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( pose.residue( i ).has_variant_type( chemical::OOP_PRE ) ) {
			oop_pre_positions.push_back( i );
		} else {
			pert_pep_mm->set_bb( i, true );
		}
	}
	
	simple_moves::SmallMoverOP pert_pep_small( new simple_moves::SmallMover( pert_pep_mm, 0.2, 1 ) );
	pert_pep_small->angle_max( 'H', 2.0 );
	pert_pep_small->angle_max( 'L', 2.0 );
	pert_pep_small->angle_max( 'E', 2.0 );
	
	pert_sequence->add_mover( pert_pep_small );
	
	if ( oop_pre_positions.size() > 0 ) {
		simple_moves::oop::OopRandomSmallMoverOP opm( new simple_moves::oop::OopRandomSmallMover ( oop_pre_positions, 2.0 ) );
		moves::RepeatMoverOP pert_pep_repeat( new moves::RepeatMover( opm, oop_pre_positions.size() * 50 ) );
		pert_sequence->add_mover( pert_pep_repeat );
	}
	moves::TrialMoverOP pert_trial( new moves::TrialMover( pert_sequence, pert_mc ) );
	
	pert_trial->apply( pose );
   	pert_mc->recover_low( pose );
}

int
main( int argc, char *argv[] )
{
	try {
		using namespace std;
		using namespace utility;
		using namespace core;
		using namespace kinematics;
		using namespace scoring;
		using namespace import_pose;
		using namespace pose;
		using namespace protocols;
		using namespace simple_moves;
		using namespace oop;
		using namespace basic;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		
		option.add( oop_conformations::step, "Degree increment" ).def(120);
		option.add( oop_conformations::dihedral_start, "Start dihedral" ).def(30);
		option.add( oop_conformations::dihedral_end, "End dihedral." ).def(360);
		option.add( oop_conformations::oop_optimize, "Do a final conformational sampling. Default true" ).def(true);

		// initialize core
		devel::init( argc, argv );

		Pose pose;

		// Make an oop.
		core::pose::make_pose_from_sequence( pose, "AAAA", "fa_standard" );
		
		utility::vector1< Size > plus_pos;
		utility::vector1< Size > empty;
		plus_pos.push_back( 1 );
		plus_pos.push_back( 3 );
		protocols::ncbb::oop::OopCreatorMoverOP OC_mover( new protocols::ncbb::oop::OopCreatorMover(plus_pos, empty, empty, empty, empty, 0, 0, false, true, false, true ));
		OC_mover->apply( pose );

		kinematics::MoveMapOP mvmp( new kinematics::MoveMap );
		ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( "mm_std_cart" );
		if ( scorefxn->has_zero_weight( core::scoring::atom_pair_constraint ) )
			scorefxn->set_weight( core::scoring::atom_pair_constraint, 1.0 );
		
		// activate if doing non AAAA
		//mvmp->set_chi( 1, true );
		
		mvmp->set_bb( true );
		protocols::simple_moves::MinMover mnmvr( mvmp, scorefxn, "lbfgs_armijo_nonmonotone", 0.0001, true );
		mnmvr.cartesian( true );
		
		// MUST minimize initially or pose is ridiculous
		mnmvr.apply( pose );

		Real step = option[oop_conformations::step].value();
		Real start = option[oop_conformations::dihedral_start].value();
		Real end = option[oop_conformations::dihedral_end].value();
		
		utility::vector1< Real > angles( 9, start );
		utility::vector1< Real > best_angles( 8, 0.0 );
		
		Real lowest_energy = 10000;
		Size p = 1;
		while ( angles[ 9 ] == start ) {
			pose.set_phi( 1, angles[1] );
			pose.set_phi( 2, angles[3] );
			pose.set_phi( 3, angles[5] );
			pose.set_phi( 4, angles[7] );
			pose.set_psi( 1, angles[2] );
			pose.set_psi( 2, angles[4] );
			pose.set_psi( 3, angles[6] );
			pose.set_psi( 4, angles[8] );
			
			std::cout << "Best energy so far is " << lowest_energy << ". Evaluating (" << angles[1] << ", " << angles[2] << "), (" << angles[3] << ", " << angles[4] << "), (" << angles[5] << ", " << angles[6] << "), (" << angles[7] << ", " << angles[8] << ")." << std::endl;
			// score the pose
			Real orig_ener = (*scorefxn)( pose );

			std::cout << "original gives " << orig_ener;
			mnmvr.apply( pose );
			std::cout << " minimizes to  " << ( (*scorefxn)( pose ) );
			
			if ( option[oop_conformations::oop_optimize].value() && ( (*scorefxn)( pose ) ) < 100 ) {
				idealize( pose, scorefxn );
				std::cout << " oop moves to  " << ( (*scorefxn)( pose ) );
				mnmvr.apply( pose );
				std::cout << " and more minimization to " << ( (*scorefxn)( pose ) );
			}
			std::cout << std::endl;
			
			Real energy = (*scorefxn)( pose );
			if ( energy < lowest_energy ) {
				lowest_energy = energy;
				Real real_phi1 = pose.phi( 1 );
				Real real_psi1 = pose.psi( 1 );
				Real real_phi2 = pose.phi( 2 );
				Real real_psi2 = pose.psi( 2 );
				Real real_phi3 = pose.phi( 3 );
				Real real_psi3 = pose.psi( 3 );
				Real real_phi4 = pose.phi( 4 );
				Real real_psi4 = pose.psi( 4 );
				std::cout << "new min " << lowest_energy << " found at (" << real_phi1 << ", " << real_psi1 << "), (" << real_phi2 << ", " << real_psi2 << "), (" << real_phi3 << ", " << real_psi3 << "), (" << real_phi4 << ", " << real_psi4 << ")!" << std::endl;
				best_angles[1] = real_phi1;
				best_angles[2] = real_psi1;
				best_angles[3] = real_phi2;
				best_angles[4] = real_psi2;
				best_angles[5] = real_phi3;
				best_angles[6] = real_psi3;
				best_angles[7] = real_phi4;
				best_angles[8] = real_psi4;
			}
			
			angles[ 1 ] += step;
			while ( angles[ p ] >= end ) {
				angles[ p ] = start;
				angles[ ++p ] += step;
				if ( angles[ p ] < end ) p = 1;
			}
			
		}
		
		std::cout << "Final lowest energy is " << lowest_energy << " found at (" << best_angles[1] << ", " << best_angles[2] << "), (" << best_angles[3] << ", " << best_angles[4] << "), (" << best_angles[5] << ", " << best_angles[6] << "), (" << best_angles[7] << ", " << best_angles[8] << ")!" << std::endl;
		
		pose.set_phi( 1, best_angles[1] );
		pose.set_phi( 2, best_angles[3] );
		pose.set_phi( 3, best_angles[5] );
		pose.set_phi( 4, best_angles[7] );
		pose.set_psi( 1, best_angles[2] );
		pose.set_psi( 2, best_angles[4] );
		pose.set_psi( 3, best_angles[6] );
		pose.set_psi( 4, best_angles[8] );
		
		pose.dump_pdb( "final_no_final_opt.pdb" );
		idealize( pose, scorefxn );
		mnmvr.apply( pose );
		pose.dump_pdb( "final_final_opt.pdb" );
		
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
