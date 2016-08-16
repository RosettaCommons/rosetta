// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief


// libRosetta headers
//#include <basic/options/option.hh>

#include <core/types.hh>


#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>

#include <devel/init.hh>

//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/kinematics/FoldTree.hh>

#include <core/io/pdb/pdb_writer.hh>

#include <core/pose/Pose.hh>


#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <utility/vector1.hh>

#include <numeric/random/random.hh>


// C++ headers
//#include <cstdlib>
#include <iostream>
#include <string>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <utility/options/keys/BooleanOptionKey.hh>


using namespace core;

using utility::vector1;


///////////////////////////////////////////////////////////////////////////////
void
scoring_test( pose::Pose & pose )
{
	using namespace core::scoring;
	ScoreFunctionOP scorefxn;
	scorefxn = new ScoreFunction();
	scorefxn->set_weight( fa_atr, 1.0 );
	scorefxn->set_weight( fa_pair, 1.0 );
	scorefxn->set_weight( fa_rep, 1.0 );

	(*scorefxn)(pose);
	/// Now handled automatically.  scorefxn->accumulate_residue_total_energies( pose );

	conformation::Residue const & rsd1 = pose.residue(1);
	conformation::Residue const & rsd2 = pose.residue(2);

	EnergyMap emap;

	scorefxn->eval_ci_2b( rsd1, rsd2, pose, emap );

	std::cout << "fa_atr between 1 and 2: " << emap[fa_atr] << std::endl;
	std::cout << "fa_pair between 1 and 2: " << emap[fa_pair] << std::endl;
	std::cout << "fa_rep between 1 and 2: " << emap[fa_rep] << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
void
small_moves_test( pose::Pose & pose )
{
	using namespace protocols::moves;
	using namespace io::pdb;
	// setup the move objects
	protocols::simple_moves::SmallMover small_mover;
	small_mover.apply( pose );
	core::Real phi = small_mover.new_phi();
	core::Real psi = small_mover.new_psi();
	std::cout << "phi: " << phi << " psi: " << psi << std::endl;

	dump_pdb( pose, "small_move.pdb" );
}

///////////////////////////////////////////////////////////////////////////////
void
shear_moves_test( pose::Pose & pose )
{
	using namespace protocols::moves;
	using namespace io::pdb;
	// setup the move objects
	protocols::simple_moves::ShearMover shear_mover;
	shear_mover.apply( pose );
	core::Real phi = shear_mover.new_phi();
	core::Real psi = shear_mover.new_psi();
	std::cout << "phi: " << phi << " psi: " << psi << std::endl;

	dump_pdb( pose, "shear_move.pdb" );
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
	using namespace chemical;
	using namespace pose;
	using namespace kinematics;
	using namespace scoring;
	using namespace io::pdb;
	// options, random initialization
	devel::init( argc, argv );

	// after devel::init to regenerate as constant seed after initial initialization
	numeric::random::RandomGenerator::initializeRandomGenerators(
		1000, numeric::random::_RND_TestRun_, "ran3" );

	Pose pose, start_pose;
	ResidueTypeSet const & residue_set( *(ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD ) ) );
	core::import_pose::pose_from_file( pose, residue_set, "input/test_in.pdb" , core::import_pose::PDB_file);

	scoring_test( pose );

	FoldTree f( pose.total_residue() );
	pose.fold_tree( f );
	std::cout << pose.fold_tree() << std::endl;

	std::cout << pose.sequence() << std::endl;

	start_pose = pose;

	small_moves_test( pose );
	pose = start_pose;
	shear_moves_test( pose );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
