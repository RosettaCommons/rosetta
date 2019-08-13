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

		Size const pose_len = 30;

		option.add(a3b_test::bin_size, "The number of degrees to step when enumerating each bond angle, e.g. 5 means 5 degree bins, or 72 per bond. Defaults to 10.").def(10);

		devel::init(argc, argv);

		//first, load the file of residue types to get min energies for.

		//Real bin_size( option[ a3b_test::bin_size ]() );
		//Size n_trials( option[ a3b_test::n_trials ]() );
		//std::cout << "Using " << bin_size << " degree bins on angles " << std::endl;
		//  bool save_best_for_last(option[a3b_test::save_best_for_last]());

		//now do initialization stuff.
		TaskFactoryOP task_factory( new TaskFactory );
		task_factory->push_back( utility::pointer::make_shared< operation::InitializeFromCommandline >() );
		//need these to keep pack_rotamers from redesigning the residue.
		operation::RestrictToRepackingOP rtrop = utility::pointer::make_shared< operation::RestrictToRepacking >();
		task_factory->push_back( rtrop );

		//ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( "r" );
		ScoreFunctionOP scorefxn = get_score_function();

		//Get the residue set we are drawing from.
		core::chemical::ResidueTypeSetCOP residue_set_cap = core::chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );

		Pose pose;  //master pose of whatever residue we are working on now.
		pose.clear();

		//ResidueType const & restype_first = residue_set_cap->name_map( "META_POLYARAMID_ALA" );
		ResidueType const & internal_META_POLYARAMID_ALA = residue_set_cap->name_map( "META_POLYARAMID_ALA" );
		ResidueType const & internal_META_POLYARAMID_CYS = residue_set_cap->name_map( "META_POLYARAMID_CYS" );
		ResidueType const & internal_META_POLYARAMID_GLY = residue_set_cap->name_map( "META_POLYARAMID_GLY" );
		ResidueType const & internal_META_POLYARAMID_GLU = residue_set_cap->name_map( "META_POLYARAMID_GLU" );
		ResidueType const & internal_META_POLYARAMID_MET = residue_set_cap->name_map( "META_POLYARAMID_MET" );
		//ResidueType const & restype_last = residue_set_cap->name_map( "ALA:MethylatedCtermProteinFull" );
		//Residue res_first( restype_first, true );
		Residue res_int_META_POLYARAMID_ALA( internal_META_POLYARAMID_ALA, true );
		Residue res_int_META_POLYARAMID_CYS( internal_META_POLYARAMID_CYS, true );
		Residue res_int_META_POLYARAMID_GLY( internal_META_POLYARAMID_GLY, true );
		Residue res_int_META_POLYARAMID_GLU( internal_META_POLYARAMID_GLU, true );
		Residue res_int_META_POLYARAMID_MET( internal_META_POLYARAMID_MET, true );
		//Residue res_int_ALA( internal_ALA, true );
		//Residue res_last( restype_last, true );
		pose.append_residue_by_jump( res_int_META_POLYARAMID_ALA, 1 );
		for ( Size ii = 1; ii <= pose_len - 1; ++ii ) {
			//pose.append_residue_by_bond( res_int_META_POLYARAMID_ALA, true );
			pose.append_residue_by_bond( res_int_META_POLYARAMID_CYS, true );
			pose.append_residue_by_bond( res_int_META_POLYARAMID_GLY, true );
			pose.append_residue_by_bond( res_int_META_POLYARAMID_GLU, true );
			pose.append_residue_by_bond( res_int_META_POLYARAMID_MET, true );
		}

		//pose.append_residue_by_bond( res_last, true );
		using namespace core::id;

		for ( Size ii = 1; ii <= pose.size() - 1; ++ii ) {
			pose.set_torsion( TorsionID( ii, BB, 5 ),  180 );
		}

		for ( Real phi = -120; phi <= 180; phi += 60 ) {

			for ( Size ii = 2; ii <= pose.size(); ++ii ) {
				pose.set_torsion( TorsionID( ii, BB, 1 ),  phi );
			}

			for ( Real psi = -120; psi <= 180; psi += 60 ) {
				for ( Size ii = 1; ii <= pose.size() - 1; ++ii ) {
					pose.set_torsion( TorsionID( ii, BB, 4 ),  psi );
				}

				auto this_pose = Pose(pose);
				std::stringstream namestr;
				namestr << "ARA_" << phi << "_" << psi << "_initial.pdb";
				this_pose.dump_pdb( namestr.str() );


				//make a minmover, let it min everything
				kinematics::MoveMapOP movemap( new kinematics::MoveMap );
				movemap->set_bb( true );
				for ( Size ii = 1; ii <= pose.size() - 1; ++ii ) {
					movemap->set_bb( ii, 5, false );
				}
				for ( Size ii = 1; ii <= pose.size(); ++ii ) {
					movemap->set_bb( ii, 2, false );
					movemap->set_bb( ii, 3, false );
				}

				movemap->set_chi( true );

				scorefxn->show(std::cout, this_pose);
				protocols::minimization_packing::MinMover minmover( movemap, scorefxn, "lbfgs_armijo_nonmonotone", 0.0001, true );//"dfpmin_strong_wolfe", 0.0001, true );
				//minmover.cartesian(true);
				minmover.apply(this_pose);
				scorefxn->show(std::cout, this_pose);

				std::stringstream name2str;
				name2str << "ARA_" << phi << "_" << psi << "_min.pdb";
				this_pose.dump_pdb( name2str.str() );
			}
		}

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;

}
