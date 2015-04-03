// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief
/// @author jk

#include <iostream>
#include <iomanip>

// Protocol Headers
#include <protocols/rigid/RigidBodyMover.hh>

// Core Headers
#include <core/conformation/Residue.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/pose/PDBInfo.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <basic/options/util.hh>
#include <core/id/AtomID_Map.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace core;
using namespace core::scoring;
using namespace core::optimization;


/// General testing code
int
main( int argc, char * argv [] )
{

	try {


	devel::init(argc, argv);

	std::cout << "Starting to build fasta sequence" << std::endl;

	// create pose for wild type native
	pose::Pose native_pose;

	//read in pdb file from command line
	std::string const input_pdb_name ( basic::options::start_file() );
	core::import_pose::pose_from_pdb( native_pose, input_pdb_name );

 	scoring::ScoreFunctionOP scorefxn( get_score_function() );
	scorefxn->set_weight( core::scoring::fa_dun, 0.1 );


	//output file
	utility::io::ozstream seq_outstream;
	seq_outstream.open( "sequence.out", std::ios::out );

	for( core::Size i = 1, nres = native_pose.total_residue(); i <= nres; i++ ) {

		std::ostringstream data_stream;
		chemical::AA const wt_aa( native_pose.residue(i).aa() );
		data_stream << oneletter_code_from_aa( wt_aa );

		seq_outstream << data_stream.str();

	}

	seq_outstream.close();
	seq_outstream.clear();

	std::cout << "Finished printing sequence" << std::endl;

	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


