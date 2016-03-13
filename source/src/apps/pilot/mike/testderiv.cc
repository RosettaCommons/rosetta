// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


// libRosetta headers

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/moves/rigid_body_moves.hh>

#include <protocols/loops/ccd_closure.hh>
#include <protocols/relax_protocols.hh>

#include <core/types.hh>

#include <core/scoring/sasa.hh>

#include <basic/prof.hh> // profiling

#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSelector.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>

//#include <core/chemical/residue_io.hh>

#include <core/scoring/etable/Etable.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>

#include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/AtomID_Map.hh>

#include <core/mm/MMTorsionLibrary.hh>
#include <core/mm/MMTorsionLibrary.fwd.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>

#include <basic/basic.hh>

#include <basic/database/open.hh>

#include <devel/init.hh>

#include <core/io/pdb/pdb_writer.hh>

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//silly using/typedef


//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <basic/Tracer.hh>
using basic::T;
using basic::Error;
using basic::Warning;


using namespace core;
using namespace protocols;

using utility::vector1;




///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
    try {
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// setup
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	devel::init(argc, argv);

	// is this also done inside devel::init?
	numeric::random::RandomGenerator::initializeRandomGenerators(
		 111111, numeric::random::_RND_TestRun_, "ran3");

	using namespace protocols::moves;
	using namespace scoring;

	pose::PoseOP pose ( new pose::Pose );
	std::cerr << "READING start.pdb" << std::endl;
	core::import_pose::pose_from_file( *pose, "start.pdb" , core::import_pose::PDB_file); // default is standard fullatom residue_set

	kinematics::MoveMap mm;
	// setup moving dofs
	//for ( int i=30; i<= 35; ++i ) {
		mm.set_bb (  true );
		mm.set_chi(  true );
	//}

	using namespace optimization;
	using pose::Pose;
	using id::AtomID;
	using id::DOF_ID;
	using id::PHI;
	using id::THETA;
	using id::D;

	scoring::ScoreFunction scorefxn;
	scorefxn.set_weight( scoring::envsmooth, 1.0 );
	//scorefxn.set_weight( scoring::fa_elec , 1.0 );
	//scorefxn.set_weight( scoring::fa_rep , 1.0 );

	(scorefxn)(*pose);
	scorefxn.show(std::cout, *pose);

	// setup the options
	MinimizerOptions options( "linmin", 10.0, true ,
													true , false );
	AtomTreeMinimizer minimizer;
	std::cout << "MINTEST: p_aa_pp" << "\n";
	std::cout << "start score: " << scorefxn( *pose ) << "\n";
	minimizer.run( *pose, mm, scorefxn, options );
	pose->dump_pdb( "min_intrares.pdb" );
	std::cout << "end score: " << scorefxn( *pose ) << "\n";

//	core::scoring::ScoreFunctionOP scorefxn( new ScoreFunction() );
//
////	// aiming for standard packer weights
//	scorefxn->set_weight( envsmooth, 0.80 );
//
//	(*scorefxn)(*pose);
//	scorefxn->show(std::cout, *pose);
//
////	using namespace optimization;
////	core::scoring::ScoreFunctionOP scorefxn( new ScoreFunction() );
////	scorefxn->set_weight( fa_atr, 0.80 );
////	scorefxn->set_weight( fa_rep, 0.44 );
////	scorefxn->set_weight( fa_sol, 0.65 );
////	scorefxn->set_weight( fa_pair, 0.49 );
////	scorefxn->set_weight( fa_dun, 0.56 );
////	scorefxn->set_weight( rama, 0.2 );
////	scorefxn->set_weight( omega, 0.5 );
////	scorefxn->set_weight( hbond_lr_bb, 1.17 );
////	scorefxn->set_weight( hbond_sr_bb, 1.17 );
////	scorefxn->set_weight( hbond_bb_sc, 1.17 );
////	scorefxn->set_weight( hbond_sc   , 1.10 );
//
//	// setup the options
//	MinimizerOptions options( "lbfgs_armijo_nonmonotone", 0.000001, true ,
//													true , false );
//
//
//	kinematics::MoveMap mm;
//	mm.set_bb (  true );
//	mm.set_chi(  true );
//	AtomTreeMinimizer minimizer;
//	std::cout << "MINTEST: p_aa_pp" << "\n";
//	std::cout << "start score: " << (*scorefxn)( *pose ) << "\n";
//	int starttime = time(NULL);
//	int startclock = clock();
//	minimizer.run( *pose, mm, *scorefxn, options );
//	int endtime = time(NULL);
//	int endclock = clock();
//	std::cout << "Time: " << (endtime - starttime) << "   " << (endclock - startclock) << std::endl;
//
//	(*scorefxn)(*pose);
//	scorefxn->show(std::cout, *pose);
    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
    return 0;
}
