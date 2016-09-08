// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   a3b_testt.cc
/// @brief  Miscellany with beta aas
/// @author Andy Watkins (amw579@nyu.edu)

// includes
#include <iostream>
#include <fstream>
#include <string>

#include <devel/init.hh>

#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/ncbb/util.hh>

#include <core/import_pose/import_pose.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/Residue.hh>

#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerMap.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/id/TorsionID.hh>
#include <core/id/types.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>

#include <protocols/simple_moves/MinMover.hh>
#include <protocols/ncbb/SecStructFinder.hh>
#include <protocols/ncbb/SecStructMinimizeMultiFunc.hh>

#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/excn/Exceptions.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <numeric/random/random.hh>

using namespace core;
using namespace utility;
using namespace scoring;
using namespace pose;
using namespace chemical;
using namespace conformation;
using namespace protocols;
using namespace basic;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::pack;
using namespace core::pack::task;
using namespace core::id;
using namespace core::scoring;
using namespace constraints;
using namespace func;

int main ( int argc, char* argv[] )
{
	try {

		devel::init(argc, argv);

		//now do initialization stuff.
		TaskFactoryOP task_factory( new TaskFactory );
		task_factory->push_back( operation::TaskOperationCOP( new operation::InitializeFromCommandline ) );
		//need these to keep pack_rotamers from redesigning the residue.
		operation::RestrictToRepackingOP rtrop = operation::RestrictToRepackingOP( new operation::RestrictToRepacking );
		task_factory->push_back( rtrop );

		//ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( "mm_std_fa_elec_dslf_fa13_split_unfolded" );
		ScoreFunctionOP scorefxn = get_score_function();

		scorefxn->set_weight( dihedral_constraint, 1.0 );
		scorefxn->set_weight( atom_pair_constraint, 1.0 );
		//Get the residue set we are drawing from.
		core::chemical::ResidueTypeSetCOP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );

		Pose pose;

		pose.append_residue_by_jump( Residue( rts->name_map( "GLN" ), 1 ), 1 );
		pose.append_residue_by_bond( Residue( rts->name_map( "B3E" ), 1 ), true );
		pose.append_residue_by_bond( Residue( rts->name_map( "GLY" ), 1 ), true );
		pose.append_residue_by_bond( Residue( rts->name_map( "PHE" ), 1 ), true );
		pose.append_residue_by_bond( Residue( rts->name_map( "SER" ), 1 ), true );
		pose.append_residue_by_bond( Residue( rts->name_map( "B3D" ), 1 ), true );
		pose.append_residue_by_bond( Residue( rts->name_map( "LEU" ), 1 ), true );
		pose.append_residue_by_bond( Residue( rts->name_map( "TRP" ), 1 ), true );
		pose.append_residue_by_bond( Residue( rts->name_map( "LYS" ), 1 ), true );
		pose.append_residue_by_bond( Residue( rts->name_map( "B3L" ), 1 ), true );
		pose.append_residue_by_bond( Residue( rts->name_map( "LEU" ), 1 ), true );
		pose.append_residue_by_bond( Residue( rts->name_map( "SER:MethylatedCtermProteinFull" ), 1 ), true );

		for ( Size ii = 1; ii <= pose.size(); ++ii ) {
			if ( pose.residue_type( ii ).is_beta_aa() ) {
				id::TorsionID bb1( ii, id::BB, 1 ); //phi
				id::TorsionID bb2( ii, id::BB, 2 ); //theta
				id::TorsionID bb3( ii, id::BB, 3 ); //psi
				id::TorsionID bb4( ii, id::BB, 4 ); //omega
				// -115, 60, -105
				// -111.207, 84.5012, -98.5471
				pose.set_torsion( bb1, -111.207  );
				pose.set_torsion( bb2,   84.5012 );
				pose.set_torsion( bb3,  -98.5471 );
				pose.set_torsion( bb4,  180 );
			} else {
				pose.set_phi(   ii, -57 );
				pose.set_psi(   ii, -48 );
				pose.set_omega( ii, 180 );
			}
			for ( Size jj = 1; jj <= pose.residue_type( ii ).nchi(); ++jj ) {
				pose.set_torsion( id::TorsionID( ii, id::CHI, jj ),  180 );
			}
		}

		kinematics::MoveMapOP movemap( new kinematics::MoveMap );
		movemap->set_bb( false );
		movemap->set_chi( true );
		protocols::simple_moves::MinMover minmover( movemap, scorefxn, "lbfgs_armijo_nonmonotone", 0.0001, true );

		pose.dump_pdb ( "B3A_initial.pdb");
		Real score = ( *scorefxn ) ( pose );
		std::cout << "Initial score is " << score << std::endl;

		minmover.apply( pose );
		Real curr_energy = ( *scorefxn ) ( pose );
		std::cout << "Minned chi:" << curr_energy << std::endl;

		pose.dump_pdb ( "B3A_chiminned.pdb");

		// Constrain h bonds
		for ( Size ii = 1; ii <= pose.size()-4; ++ii ) {

			std::string ca = "CA";
			if ( pose.residue_type( ii ).is_beta_aa() ) {
				ca = "CM";
			}
			pose.add_constraint( DihedralConstraintOP( new DihedralConstraint(
				*new AtomID( pose.residue(  ii  ).atom_index(  ca  ), ii ),
				*new AtomID( pose.residue(  ii  ).atom_index( "C" ), ii ),
				*new AtomID( pose.residue(  ii+1  ).atom_index( "N" ), ii+1 ),
				*new AtomID( pose.residue(  ii+1  ).atom_index( "CA" ), ii+1 ),
				CircularHarmonicFuncOP( new CircularHarmonicFunc( 3.14159, 0.01 ) )
				) ) );


			pose.add_constraint( AtomPairConstraintOP( new AtomPairConstraint(
				*new AtomID( pose.residue(  ii  ).atom_index( "O" ), ii ),
				*new AtomID( pose.residue( ii+4 ).atom_index( "H" ), ii ),
				HarmonicFuncOP( new HarmonicFunc( 1.8, 0.2 ) ) ) ) );
		}

		movemap->set_bb( true );
		minmover.apply( pose );

		std::cout << "Pose with torsions:" << std::endl;
		for ( Size ii = 1; ii <= pose.size(); ++ii ) {
			for ( Size jj = 1; jj <= pose.residue( ii ).mainchain_torsions().size(); ++jj ) {
				std::cout << pose.residue( ii ).mainchain_torsions()[ jj ] << ", ";
			}
			std::cout << std::endl;
		}
		std::cout << "has energy " << curr_energy << std::endl;
		pose.dump_pdb ( "B3A_bbminned.pdb");

		// Replace first and third residue.
		ResidueType const & pre_type = rts->get_residue_type_with_variant_added( pose.residue_type( 1 ), chemical::HBS_PRE );
		ResidueType const & post_type = rts->get_residue_type_with_variant_added( pose.residue_type( 3 ), chemical::HBS_POST );

		// It's losing the chis of the original residue. Stupid stupid stupid.
		utility::vector1< Real > old_chis;
		for ( Size ii = 1; ii <= pose.residue_type( ii ).nchi(); ++ii ) {
			old_chis.push_back( pose.residue( ii ).chi( ii ) );
		}
		pose.replace_residue( 1, Residue( pre_type, 1 ), true );
		for ( Size ii = 1; ii <= pose.residue_type( 1 ).nchi(); ++ii ) {
			pose.set_torsion( id::TorsionID( 1, id::CHI, ii ),  old_chis[ ii ] );
		}
		pose.replace_residue( 3, Residue( post_type, 1 ), true );
		core::pose::ncbb::initialize_ncbbs( pose );
		pose.conformation().declare_chemical_bond( 1, "CYH", 3, "CZH" );

		pose.set_torsion( id::TorsionID( 1, id::BB, 5 ), -60 );


		pose.dump_pdb ( "B3A_hbsed.pdb");

		movemap->set_bb( false );
		movemap->set_chi( false );

		movemap->set_bb( 1, true );
		movemap->set_bb( 2, true );
		movemap->set_bb( 3, true );
		minmover.apply( pose );

		std::cout << pose.residue_type( 1 );

		std::cout << pose.residue( 1 );
		pose.dump_pdb ( "B3A_final.pdb");


		// New way.
		// Cribbed from the SecStructFinder, but I don't wanna loop.
		Pose minpose( pose );
		core::optimization::MinimizerMap min_map;
		min_map.setup( minpose, *movemap );

		( *scorefxn ) ( minpose );
		protocols::ncbb::SecStructMinimizeMultiFunc ssmmf( minpose, *scorefxn, min_map, "AAABAAA", "AAABAAA" );

		core::optimization::MinimizerOptions minoptions( "lbfgs_armijo_nonmonotone", 0.0001, true, false, false ); // investigate the final bools?
		minpose.energies().set_use_nblist( minpose, min_map.domain_map(), false );
		core::optimization::Minimizer minimizer( ssmmf, minoptions );

		//minmover.apply ( minpose );
		utility::vector1< Real > dihedrals_for_minimization;
		dihedrals_for_minimization.push_back( -57 );
		dihedrals_for_minimization.push_back( -48 );
		dihedrals_for_minimization.push_back(-111.207 );
		dihedrals_for_minimization.push_back(  84.5012 );
		dihedrals_for_minimization.push_back( -98.5471 );

		minimizer.run( dihedrals_for_minimization );

		minpose.energies().reset_nblist();

		scorefxn->set_weight( core::scoring::dihedral_constraint, 0.0 );
		Real minscore = ( ( *scorefxn ) ( minpose ) );

		std::cout << " and minned is " << minscore << std::endl;
		minpose.dump_pdb ( "B3A_real.pdb");



	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}
