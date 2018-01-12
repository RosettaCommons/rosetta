// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   cov_hbs.cc
/// @brief  Sidechain conjugation to acryl amides
/// @author Andy Watkins (amw579@nyu.edu)

// includes
#include <iostream>
#include <fstream>
#include <string>

#include <devel/init.hh>

#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/ncbb/util.hh>

#include <core/import_pose/import_pose.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/Residue.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/id/TorsionID.hh>
#include <core/id/types.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/func/FadeFunc.hh>
#include <core/scoring/func/SumFunc.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Numeric Headers
#include <numeric/conversions.hh>

// Mover headers
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/RotamerTrialsMover.hh>
#include <protocols/minimization_packing/TaskAwareMinMover.hh>
#include <protocols/simple_moves/BackboneMover.fwd.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/RandomTorsionMover.hh>
#include <protocols/ncbb/a3b_hbs/A3BHbsPatcher.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/excn/Exceptions.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <numeric/random/random.hh>

using namespace protocols;
using namespace protocols::moves;
using namespace protocols::simple_moves;
using namespace protocols::minimization_packing;
using namespace protocols::simple_moves::a3b_hbs;

using namespace basic;
using namespace basic::options;
using namespace basic::options::OptionKeys;


void
add_hbond_constraints( core::pose::Pose & pose ) {

	using namespace core;
	using namespace id;
	using namespace scoring;
	using namespace constraints;
	using namespace func;

	HarmonicFuncOP harm( new core::scoring::func::HarmonicFunc( 2, .5 ) );
	CircularHarmonicFuncOP pi( new core::scoring::func::CircularHarmonicFunc( 3.14159, .05 ) );
	CircularHarmonicFuncOP zero( new core::scoring::func::CircularHarmonicFunc( 0, .05 ) );

	pose.add_constraint( AtomPairConstraintOP( new AtomPairConstraint(
		*new AtomID( pose.residue( 1 ).atom_index( "NL1" ), 1 ),
		*new AtomID( pose.residue( 4 ).atom_index( "1HI2" ), 4 ), harm ) ) );
	pose.add_constraint( AtomPairConstraintOP( new AtomPairConstraint(
		*new AtomID( pose.residue( 1 ).atom_index( "1HL2" ), 1 ),
		*new AtomID( pose.residue( 4 ).atom_index( "OL" ), 4 ) , harm ) ) );

	pose.add_constraint( DihedralConstraintOP( new DihedralConstraint(
		*new AtomID( pose.residue( 4 ).atom_index( "OL" ), 4 ),
		*new AtomID( pose.residue( 1 ).atom_index( "1HL2" ), 1 ),
		*new AtomID( pose.residue( 1 ).atom_index( "2HL2" ), 1 ),
		*new AtomID( pose.residue( 1 ).atom_index( "CK2" ), 1 ), pi ) ) );

	pose.add_constraint( DihedralConstraintOP( new DihedralConstraint(
		*new AtomID( pose.residue( 4 ).atom_index( "1HI2" ), 4 ),
		*new AtomID( pose.residue( 1 ).atom_index( "NL1" ), 1 ),
		*new AtomID( pose.residue( 1 ).atom_index( "CK1" ), 1 ),
		*new AtomID( pose.residue( 1 ).atom_index( "CK2" ), 1 ), pi ) ) );
	/*
	pose.add_constraint( AtomPairConstraintOP( new AtomPairConstraint(
	*new AtomID( pose.residue( 2 ).atom_index( "OI1" ), 2 ),
	*new AtomID( pose.residue( 9 ).atom_index( "1HL3" ), 9 ) , harm ) ) );
	pose.add_constraint( AtomPairConstraintOP( new AtomPairConstraint(
	*new AtomID( pose.residue( 2 ).atom_index( "NI2" ), 2 ),
	*new AtomID( pose.residue( 9 ).atom_index( "1HL2" ), 9 ) , harm ) ) );
	pose.add_constraint( AtomPairConstraintOP( new AtomPairConstraint(
	*new AtomID( pose.residue( 2 ).atom_index( "1HL" ), 2 ),
	*new AtomID( pose.residue( 9 ).atom_index( "OL1" ), 9 ) , harm ) ) );

	pose.add_constraint( AtomPairConstraintOP( new AtomPairConstraint(
	*new AtomID( pose.residue( 3 ).atom_index( "1HL3" ), 3 ),
	*new AtomID( pose.residue( 8 ).atom_index( "OI1" ), 8 ) , harm ) ) );
	pose.add_constraint( AtomPairConstraintOP( new AtomPairConstraint(
	*new AtomID( pose.residue( 3 ).atom_index( "1HL2" ), 3 ),
	*new AtomID( pose.residue( 8 ).atom_index( "NI2" ), 8 ) , harm ) ) );
	pose.add_constraint( AtomPairConstraintOP( new AtomPairConstraint(
	*new AtomID( pose.residue( 3 ).atom_index( "OL1" ), 3 ),
	*new AtomID( pose.residue( 8 ).atom_index( "1HL" ), 8 ) , harm ) ) );

	pose.add_constraint( AtomPairConstraintOP( new AtomPairConstraint(
	*new AtomID( pose.residue( 4 ).atom_index( "OL" ), 4 ),
	*new AtomID( pose.residue( 7 ).atom_index( "2HL2" ), 7 ) , harm ) ) );
	pose.add_constraint( AtomPairConstraintOP( new AtomPairConstraint(
	*new AtomID( pose.residue( 4 ).atom_index( "1HI2" ), 4 ),
	*new AtomID( pose.residue( 7 ).atom_index( "NL1" ), 7 ) , harm ) ) );

	pose.add_constraint( AtomPairConstraintOP( new AtomPairConstraint(
	*new AtomID( pose.residue( 5 ).atom_index( "1HI2" ), 5 ),
	*new AtomID( pose.residue( 6 ).atom_index( "NL1" ), 6 ) , harm ) ) );
	pose.add_constraint( AtomPairConstraintOP( new AtomPairConstraint(
	*new AtomID( pose.residue( 5 ).atom_index( "OI1" ), 5 ),
	*new AtomID( pose.residue( 6 ).atom_index( "2HL2" ), 6 ) , harm ) ) );
	*/
}

int main ( int argc, char* argv[] )
{
	try {
		//option[ chemical::patch_selectors ].push_back( "CTERM_AMIDATION" );

		devel::init(argc, argv);

		using namespace core;
		using namespace utility;
		using namespace scoring;
		using namespace pose;
		using namespace core::chemical;
		using namespace conformation;
		using namespace func;
		using namespace constraints;

		using namespace core::id;
		using namespace core::pack;
		using namespace core::pack::task;

		Pose pose;

		ScoreFunctionOP scorefxn = get_score_function();
		ScoreFunctionOP guaranteed_cart_scorefxn = ScoreFunctionFactory::create_score_function( "talaris2013_cart" );
		//scorefxn->set_weight( hbond_sc, 5 );

		//Get the residue set we are drawing from.
		core::chemical::ResidueTypeSetCOP residue_set_cap = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

		ResidueType const & pna = residue_set_cap->name_map( "APN" );
		ResidueType const & pnc = residue_set_cap->name_map( "CPN" );
		ResidueType const & png = residue_set_cap->name_map( "GPN" );
		ResidueType const & pnt = residue_set_cap->name_map( "TPN" );
		ResidueType const & pnu = residue_set_cap->name_map( "UPN" );

		Residue res_pna_nterm( residue_set_cap->name_map( "APN:NtermProteinFull" ), true );
		Residue res_pna_cterm( residue_set_cap->name_map( "APN:CtermProteinFull" ), true );
		Residue res_pna( pna, true );
		Residue res_pnc( pnc, true );
		Residue res_png( png, true );
		Residue res_pnt( pnt, true );
		Residue res_pnu( pnu, true );
		Residue res_pnu_nterm( residue_set_cap->name_map( "UPN:NtermProteinFull" ), true );
		Residue res_pnu_cterm( residue_set_cap->name_map( "UPN:CtermProteinFull" ), true );

		pose.append_residue_by_jump( res_pna, 1 );
		pose.append_residue_by_bond( res_pnc, true );
		pose.append_residue_by_jump( res_png, 2 );
		pose.append_residue_by_bond( res_pnu, true );

		//std::cout << pose.fold_tree();

		//pose.pdb_info()->attach_to( pose.conformation() );
		//pose.pdb_info()->show(std::cout);
		/*pose.append_residue_by_jump( res_pna_nterm, 1 );
		pose.append_residue_by_bond( res_pnc, true );
		pose.append_residue_by_bond( res_png, true );
		pose.append_residue_by_bond( res_pnt, true );
		pose.append_residue_by_bond( res_pnu, true );
		pose.append_residue_by_jump( res_pna, 2 );
		pose.append_residue_by_bond( res_pna, true );
		pose.append_residue_by_bond( res_pnc, true );
		pose.append_residue_by_bond( res_png, true );
		pose.append_residue_by_bond( res_pnu, true );
		*/
		protocols::rigid::RigidBodyTransMover trans_mover( pose, 1 );
		trans_mover.step_size( 4 );
		trans_mover.apply( pose );
		rigid::RigidBodyPerturbMoverOP big_pert_dock_rbpm( new rigid::RigidBodyPerturbMover(1, 1.0, 1.5 ) );
		rigid::RigidBodyPerturbMoverOP pert_dock_rbpm( new rigid::RigidBodyPerturbMover(1, 1.0, 0.5 ) );

		add_hbond_constraints( pose );

		kinematics::MoveMapOP mm( new kinematics::MoveMap );
		mm->set_bb( true );

		kinematics::MoveMapOP min_mm( new kinematics::MoveMap );
		//min_mm->set_chi( true );
		min_mm->set_jump( true );

		// SET INITIAL CONFORMATION
		for ( Size resi = 1; resi <= pose.size(); ++resi ) {
			mm->set( id::TorsionID( resi, id::BB, 5 ), false );
			min_mm->set( id::TorsionID( resi, id::BB, 5 ), false );

			pose.set_torsion( id::TorsionID( resi, id::BB, 1 ), 70 );
			pose.set_torsion( id::TorsionID( resi, id::BB, 2 ), 70 );
			pose.set_torsion( id::TorsionID( resi, id::BB, 3 ), 95 );
			pose.set_torsion( id::TorsionID( resi, id::BB, 4 ), -5 );
			pose.set_torsion( id::TorsionID( resi, id::BB, 5 ), 180 );

			pose.set_torsion( id::TorsionID( resi, id::CHI, 1 ), 0 );
			pose.set_torsion( id::TorsionID( resi, id::CHI, 2 ), 180 );
			pose.set_torsion( id::TorsionID( resi, id::CHI, 3 ), -90 );

		}

		RandomTorsionMoverOP small_tor_mover( new RandomTorsionMover( mm, 1, 10 ) );
		MinMoverOP min( new MinMover( min_mm, guaranteed_cart_scorefxn, "lbfgs_armijo_nonmonotone", 0.001, true ) );
		min->cartesian( true );
		moves::MonteCarloOP pert_mc( new moves::MonteCarlo( pose, *scorefxn, 3 ) );

		min->apply( pose );
		pose.dump_pdb( "min.pdb" );

		Real apc = 0.1;
		scorefxn->set_weight( atom_pair_constraint, apc );
		scorefxn->set_weight( angle_constraint, apc );
		scorefxn->set_weight( dihedral_constraint, apc );

		Real best_score = 100000;
		for ( Size ii = 1; ii <= 1000; ++ii ) {
			std::cout << "PNA Creator: Round " << ii << " / 1000" << std::endl;
			std::cout << "PNA Creator: score " << ( *scorefxn )( pose ) << " best " << best_score << std::endl;

			for ( Size jj = 1; jj <= 10; ++jj ) {
				small_tor_mover->apply( pose );
				//big_pert_dock_rbpm->apply( pose );
			}

			if ( pert_mc->boltzmann( pose ) ) {
				min->apply( pose );
				Real test_score( pose.energies().total_energies().dot( scorefxn->weights() ) );
				if ( test_score <= best_score ) {
					best_score = test_score;
					pert_mc->reset( pose );
				}
			}

			if ( ! ( ii % 100 ) ) { pose.dump_pdb( "out_" + utility::to_string( ii ) + ".pdb" ); apc += 0.1; }
		}
		//pert_mc->recover_low( pose );


		std::cout << "about to dump" << std::endl;
		pose.dump_pdb("out.pdb" );

		std::cout << "PNA Creator: Removing constraints and remodeling" << std::endl;
		core::scoring::constraints::ConstraintCOPs cs = pose.constraint_set()->get_all_constraints();
		for ( Size ii = 1; ii <= cs.size(); ii++ ) {
			pose.remove_constraint( cs[ii], true );
		}

		best_score = 100000;
		for ( Size ii = 1; ii <= 100; ++ii ) {
			std::cout << "PNA Creator: Round " << ii << " / 100" << std::endl;
			std::cout << "PNA Creator: score " << ( *scorefxn )( pose ) << " best " << best_score << std::endl;

			for ( Size jj = 1; jj <= 100; ++jj ) {
				small_tor_mover->apply( pose );
				//pert_dock_rbpm->apply( pose );
			}

			if ( pert_mc->boltzmann( pose ) ) {
				Real test_score( pose.energies().total_energies().dot( scorefxn->weights() ) );
				if ( test_score <= best_score ) {
					min->apply( pose );
					best_score = test_score;
					pert_mc->reset( pose );
				}
			}

			if ( ! ( ii % 10 ) ) pose.dump_pdb( "new_" + utility::to_string( ii ) + ".pdb" );
		}
		pert_mc->recover_low( pose );
		pose.dump_pdb( "new.pdb" );

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
