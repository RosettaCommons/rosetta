// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license.
// (c) The Rosetta software is developed by the contributing members of the
// (c) Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org.
// (c) Questions about this can be addressed to University of Washington UW
// (c) TechTransfer, email: license@u.washington.edu.

/// @file   build_a3b.cc
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

#include <core/kinematics/MoveMap.hh>

#include <core/id/TorsionID.hh>
#include <core/id/types.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>

#include <protocols/simple_moves/MinMover.hh>

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

static numeric::random::RandomGenerator RG(6117547);

int main ( int argc, char* argv[] )
{
    devel::init(argc, argv);

	//first, load the file of residue types to get min energies for.

	//now do initialization stuff.
	TaskFactoryOP task_factory = new TaskFactory;
	task_factory->push_back( new operation::InitializeFromCommandline );
	//need these to keep pack_rotamers from redesigning the residue.
	operation::RestrictToRepackingOP rtrop = new operation::RestrictToRepacking;
	task_factory->push_back( rtrop );

	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( "mm_std" );

	//Get the residue set we are drawing from.
	core::chemical::ResidueTypeSetCAP residue_set_cap = core::chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );

	Pose pose;  //master pose of whatever residue we are working on now.
	pose.clear();

	ResidueType const & restype_first = residue_set_cap->name_map( "PHE:NtermProteinFull:a3b_hbs_pre" );
	ResidueType const & internal_B3A = residue_set_cap->name_map( "B3A" );
	ResidueType const & internal_ALA = residue_set_cap->name_map( "ALA" );
	ResidueType const & internal_GLY = residue_set_cap->name_map( "GLY:a3b_hbs_post" );
	ResidueType const & restype_last = residue_set_cap->name_map( "ASN:MethylatedCtermProteinFull" );
	Residue res_first( restype_first, true );
	Residue res_int_B3A( internal_B3A, true );
	Residue res_int_ALA( internal_ALA, true );
	Residue res_int_GLY( internal_GLY, true );
	Residue res_last( restype_last, true );
	pose.append_residue_by_jump( res_first, 1 );
	pose.append_residue_by_bond( res_int_ALA, true );
	pose.append_residue_by_bond( res_int_GLY, true );
	pose.append_residue_by_bond( res_int_B3A, true );
	pose.append_residue_by_bond( res_int_ALA, true );
	pose.append_residue_by_bond( res_int_ALA, true );
	pose.append_residue_by_bond( res_int_ALA, true );
	pose.append_residue_by_bond( res_int_B3A, true );
	pose.append_residue_by_bond( res_int_ALA, true );
	pose.append_residue_by_bond( res_int_ALA, true );
	pose.append_residue_by_bond( res_int_ALA, true );
	pose.append_residue_by_bond( res_int_B3A, true );
	pose.append_residue_by_bond( res_int_ALA, true );
	pose.append_residue_by_bond( res_int_ALA, true );
	pose.append_residue_by_bond( res_int_ALA, true );
	pose.append_residue_by_bond( res_last, true );

    for ( Size i = 1; i <= pose.n_residue(); ++i ) {
        if ( pose.residue(i).type().is_beta_aa() ) { //( i == 4 || i == 8 || i == 12 ) {
            id::TorsionID bb1( i, id::BB, 1 ); //phi
            id::TorsionID bb2( i, id::BB, 2 ); //theta
            id::TorsionID bb3( i, id::BB, 3 ); //psi
            id::TorsionID bb4( i, id::BB, 4 ); //omg
            
            pose.set_torsion( bb1, -104 );
            pose.set_torsion( bb2, 64 );
            pose.set_torsion( bb3, -116 );
            pose.set_torsion( bb4, 180 );
        } else {
            pose.set_phi( i, -57);
            pose.set_omega( i, 180);
            pose.set_psi( i, -48);
        }
    }
    
    core::pose::ncbb::add_a3b_hbs_constraint( pose, 1 );
	
	//make a minmover, let it min everything
    std::cout << "Movemap creating" << std::endl;
	kinematics::MoveMapOP movemap( new kinematics::MoveMap );
    std::cout << "BB true" << std::endl;
	movemap->set_bb( true);
    std::cout << "Minmover creating" << std::endl;
	protocols::simple_moves::MinMover minmover( movemap, scorefxn, "lbfgs_armijo_nonmonotone", 0.0001, true );//"dfpmin_strong_wolfe", 0.0001, true );

	// iterate over possible sets of bond angles, test them all! Record the best.
    std::cout << "Dump initial" << std::endl;
	pose.dump_pdb ( "B3A_initial.pdb");
    std::cout << "Score" << std::endl;
	Real score = ( *scorefxn ) ( pose );
    std::cout << "Minimize" << std::endl;
	minmover.apply ( pose );
    std::cout << "Dump final" << std::endl;
	pose.dump_pdb ( "B3A_min.pdb");
}
