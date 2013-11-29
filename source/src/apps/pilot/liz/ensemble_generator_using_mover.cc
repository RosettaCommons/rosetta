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
/// @author Liz Kellogg ekellogg@u.washington.edu

// libRosetta headers




#include <core/types.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/scoring/sasa.hh>
#include <core/scoring/rms_util.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
//#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/Func.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/PDBInfo.hh>
#include <core/io/pdb/pose_io.hh>

#include <devel/init.hh>

#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/keys/OptionKeys.hh>

#include <basic/basic.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/ProteinSilentStruct.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

//protocols
//#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/ShakeStructureMover.hh>

#include <utility/file/FileName.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>
#include <time.h>
using basic::T;
using basic::Warning;
using basic::Error;

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

//Auto Headers
#include <core/import_pose/import_pose.hh>

numeric::random::RandomGenerator RG(15434); // <- Magic number, do not change it!!!
int
main( int argc, char* argv [] )
{
    try {
	using namespace core;
	using namespace core::pose;
	using namespace utility;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace protocols::moves;
	// options, random initialization. MAKE SURE THIS COMES BEFORE OPTIONS
	devel::init( argc, argv );
	pose::Pose pose;


	vector1<file::FileName> files;
	if(basic::options::option[in::file::s].user()){
		std::cout << "using -s option" << std::endl;
		files=option[in::file::s]();
	}else if( option[in::file::l].user()){
		std::cout << "using -l option " << std::endl;

		utility::vector1<file::FileName> list = basic::options::option[ in::file::l ]();
		for(unsigned int h=1;h<=list.size();h++){
			utility::io::izstream pdbs(list[h]);
			std::string fname;
			while(pdbs >> fname){
				files.push_back(fname);
			}
		}
	}
	for(unsigned int f=1; f<=files.size();f++){
		core::import_pose::pose_from_pdb(pose, files[f]);
		int num_struct = 10;
		ShakeStructureMover ssm; //don't initialize with a scorefunction,
		//scorefunction with rama,omega,fa_dun,p_aa_pp terms will be created for you
		ssm.set_sc_min(true); //sc min after perturbing backbone
		//ssm.set_mc_temperature(5);
		ssm.set_ensemble_diversity(1.0);
		ssm.set_rmsd_target_tolerance(0.1);
		std::ostringstream curr;

		while(num_struct > 0){
		  std::cout << "structure num. " << num_struct << std::endl;
		  pose::Pose init(pose);
		  ssm.apply(init);
		  //dump init
		  num_struct--;
			curr << num_struct;
			init.dump_pdb("test_"+curr.str()+".pdb");
		}

	}
    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
    }
    return 0;
}
