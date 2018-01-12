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

#include <protocols/minimization_packing/MinMover.hh>

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

// init options
namespace a3b_test
{
RealOptionKey const bin_size("a3b_test::bin_size"); //the number of degrees to step when enumerating each bond angle, e.g. 10 means 10 degree bins.
IntegerOptionKey const n_trials("a3b_test::n_trials"); //the number of degrees to step when enumerating each bond angle, e.g. 10 means 10 degree bins.
}

int main ( int argc, char* argv[] )
{
	try {
		option.add(a3b_test::bin_size, "The number of degrees to step when enumerating each bond angle, e.g. 5 means 5 degree bins, or 72 per bond. Defaults to 10.").def(10);

		devel::init(argc, argv);

		//first, load the file of residue types to get min energies for.

		Real bin_size( option[ a3b_test::bin_size ]() );
		Size n_trials( option[ a3b_test::n_trials ]() );
		std::cout << "Using " << bin_size << " degree bins on angles " << std::endl;
		//  bool save_best_for_last(option[a3b_test::save_best_for_last]());

		//now do initialization stuff.
		TaskFactoryOP task_factory( new TaskFactory );
		task_factory->push_back( operation::TaskOperationCOP( new operation::InitializeFromCommandline ) );
		//need these to keep pack_rotamers from redesigning the residue.
		operation::RestrictToRepackingOP rtrop = operation::RestrictToRepackingOP( new operation::RestrictToRepacking );
		task_factory->push_back( rtrop );

		ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( "mm_std" );

		//Get the residue set we are drawing from.
		core::chemical::ResidueTypeSetCOP residue_set_cap = core::chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );

		Pose pose;  //master pose of whatever residue we are working on now.
		pose.clear();

		ResidueType const & restype_first = residue_set_cap->name_map( "ALA:AcetylatedNtermProteinFull" );
		ResidueType const & internal_B3A = residue_set_cap->name_map( "B3A" );
		ResidueType const & internal_ALA = residue_set_cap->name_map( "ALA" );
		ResidueType const & restype_last = residue_set_cap->name_map( "ALA:MethylatedCtermProteinFull" );
		Residue res_first( restype_first, true );
		Residue res_int_B3A( internal_B3A, true );
		Residue res_int_ALA( internal_ALA, true );
		Residue res_last( restype_last, true );
		pose.append_residue_by_jump( res_first, 1 );
		pose.append_residue_by_bond( res_int_ALA, true );
		pose.append_residue_by_bond( res_int_ALA, true );
		pose.append_residue_by_bond( res_int_B3A, true );
		pose.append_residue_by_bond( res_int_ALA, true );
		pose.append_residue_by_bond( res_int_ALA, true );
		pose.append_residue_by_bond( res_last, true );

		pose.set_phi( 2, -57);
		pose.set_phi( 3, -57);
		pose.set_phi( 5, -57);
		pose.set_phi( 6, -57);
		pose.set_phi( 7, -57);
		pose.set_psi( 1, -48);
		pose.set_psi( 2, -48);
		pose.set_psi( 3, -48);
		pose.set_psi( 5, -48);
		pose.set_psi( 6, -48);
		pose.set_omega( 1, 180);
		pose.set_omega( 2, 180);
		pose.set_omega( 3, 180);
		pose.set_omega( 4, 180);
		pose.set_omega( 5, 180);
		pose.set_omega( 6, 180);
		pose.set_omega( 7, 180);

		id::TorsionID bb1( 4, id::BB, 1 ); //phi
		id::TorsionID bb2( 4, id::BB, 2 ); //theta
		id::TorsionID bb3( 4, id::BB, 3 ); //psi

		//make a minmover, let it min everything
		kinematics::MoveMapOP movemap( new kinematics::MoveMap );
		movemap->set_bb( true);
		protocols::minimization_packing::MinMover minmover( movemap, scorefxn, "lbfgs_armijo_nonmonotone", 0.0001, true );//"dfpmin_strong_wolfe", 0.0001, true );

		// iterate over possible sets of bond angles, test them all! Record the best.

		pose.dump_pdb ( "B3A_initial.pdb");
		Real score = ( *scorefxn ) ( pose );
		std::cout << "Initial score is " << score << std::endl;
		Real best_score = 1000000;
		Pose best_pose = pose;
		utility::vector1< Pose > best_poses;
		for ( Size i = 1; i <= n_trials; ++i ) {
			Pose currpose = pose;
			currpose.set_torsion( bb1, numeric::random::uniform() * 360 - 180 );
			currpose.set_torsion( bb2, numeric::random::uniform() * 360 - 180 );
			currpose.set_torsion( bb3, numeric::random::uniform() * 360 - 180 );
			minmover.apply ( currpose );
			Real curr_energy = ( *scorefxn ) ( currpose );
			std::cout << "Pose with torsions ( " << currpose.torsion( bb1 ) << ", " << currpose.torsion( bb2 ) << ", " << currpose.torsion( bb3 ) << " ) has energy " << curr_energy << std::endl;
			if ( curr_energy < best_score ) {
				best_score = curr_energy;
				best_pose = currpose;
				if ( best_poses.size() < n_trials/10 ) {
					best_poses.push_back( currpose );
				} else {
					for ( Size j = 1; j <= n_trials/10; ++j ) {
						if ( ( *scorefxn ) ( best_poses[ j ] ) > ( curr_energy ) ) {
							best_poses[ j ] = currpose;
						}
					}
				}
			}
		}

		best_pose.dump_pdb ( "B3A_best.pdb");
		for ( Size j = 1; j <= n_trials/10; ++j ) {
			std::stringstream fn;
			fn << "B3A_best_" << j << ".pdb";
			best_pose.dump_pdb( fn.str() );
		}
	} catch (utility::excn::Exception const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}
