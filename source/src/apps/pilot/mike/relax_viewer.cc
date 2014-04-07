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
/// @author Mike Tyka
/// @author James Thompson

// libRosetta headers

#include <protocols/viewer/viewers.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
// AUTO-REMOVED #include <protocols/simple_moves/PackRotamersMover.hh>
// AUTO-REMOVED #include <protocols/moves/rigid_body_moves.hh>

// AUTO-REMOVED #include <protocols/loops/ccd_closure.hh>
// AUTO-REMOVED #include <protocols/relax_protocols.hh>


#include <core/types.hh>

// AUTO-REMOVED #include <core/scoring/sasa.hh>

// AUTO-REMOVED #include <basic/prof.hh> // profiling

// AUTO-REMOVED #include <core/chemical/AtomTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/MMAtomTypeSet.hh>

#include <core/chemical/AA.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/ResidueSelector.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>
// AUTO-REMOVED #include <core/chemical/VariantType.hh>
// AUTO-REMOVED

// AUTO-REMOVED #include <core/scoring/etable/Etable.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/Ramachandran.hh>
// AUTO-REMOVED #include <core/pack/dunbrack/RotamerLibrary.hh>
// AUTO-REMOVED #include <core/scoring/hbonds/HBondSet.hh>
// AUTO-REMOVED #include <core/scoring/hbonds/hbonds.hh>
// AUTO-REMOVED #include <core/scoring/etable/count_pair/CountPairFunction.hh>

// AUTO-REMOVED #include <core/pack/rotamer_trials.hh>
// AUTO-REMOVED #include <core/pack/pack_rotamers.hh>
// AUTO-REMOVED #include <core/pack/task/PackerTask.hh>
// AUTO-REMOVED #include <core/pack/task/TaskFactory.hh>

// AUTO-REMOVED #include <core/kinematics/FoldTree.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
// AUTO-REMOVED #include <core/id/AtomID_Map.hh>
// AUTO-REMOVED #include <core/id/AtomID_Map.Pose.hh>

#include <core/pose/Pose.hh>

#include <basic/options/util.hh> //option.hh>
//#include <basic/options/after_opts.hh>

// AUTO-REMOVED #include <basic/basic.hh>
#include <basic/Tracer.hh>
// AUTO-REMOVED #include <basic/database/open.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

// AUTO-REMOVED #include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>


// C++ headers
// AUTO-REMOVED #include <fstream>
#include <iostream>
#include <string>


// option key includes

// AUTO-REMOVED #include <basic/options/keys/phil.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>
#include <protocols/relax/ClassicRelax.hh>



using basic::T;
using basic::Error;
using basic::Warning;

using namespace core;
using namespace protocols;

using utility::vector1;
using io::pdb::dump_pdb;


///////////////////////////////////////////////////////////////////////////////
void
relax_test()
{
	using namespace protocols::moves;
	using namespace scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// std::string pdbfile = "1rdgA.pdb";
	std::string pdbfile = basic::options::start_file();//option[ phil::s ]();
	std::cout << "relax_viewer!" << std::endl;

	pose::Pose pose;
	std::cerr << "READING " << pdbfile << std::endl;
	core::import_pose::pose_from_pdb( pose, pdbfile ); // default is standard fullatom residue_set

	std::cerr << "SETUP SCORE FUNCTION" << std::endl;
	core::scoring::ScoreFunctionOP scorefxn( new ScoreFunction() );

	protocols::viewer::add_conformation_viewer( pose.conformation(), "start_pose" );

	// aiming for standard packer weights
	scorefxn->set_weight( fa_atr, 0.80 );
	scorefxn->set_weight( fa_rep, 0.44 );
	scorefxn->set_weight( fa_sol, 0.65 );
	scorefxn->set_weight( fa_pair, 0.49 );
	scorefxn->set_weight( fa_dun, 0.56 );
	scorefxn->set_weight( rama, 0.2 );
	scorefxn->set_weight( hbond_lr_bb, 1.17 );
	scorefxn->set_weight( hbond_sr_bb, 1.17 );
	scorefxn->set_weight( hbond_bb_sc, 1.17 );
	scorefxn->set_weight( hbond_sc   , 1.10 );

	protocols::relax::ClassicRelax myrelax( scorefxn);
	moves::MonteCarloOP mc = myrelax.get_mc(pose );
	protocols::viewer::add_monte_carlo_viewer( *mc );

	std::cerr << "RUNNING RELAX" << std::endl;
	myrelax.apply(pose );
	std::cerr << "DONE" << std::endl;
}



///////////////////////////////////////////////////////////////////////////////

void*
my_main( void* )
{
	numeric::random::RandomGenerator::initializeRandomGenerators(
		 111111, numeric::random::_RND_TestRun_, "ran3");

	relax_test();
	return 0;
}

int
main( int argc, char * argv [] )
{
    try {
    	// options, random initialization
    	devel::init( argc, argv );
    	protocols::viewer::viewer_main( my_main );
    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
    return 0;
}
