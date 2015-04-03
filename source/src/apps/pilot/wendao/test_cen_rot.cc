// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file test_bbmc.cc
/// @brief test centroid rot model
/// @author Yuan Liu


// Core Headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>

#include <core/chemical/AA.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>

// 
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/dunbrack/cenrot/SingleResidueCenrotLibrary.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
//#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/Minimizer.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <numeric/constants.hh>
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>


#include <protocols/simple_moves/BBGaussianMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/PyMolMover.hh>

#include <sstream>
#include <iostream>
#include <string>

using namespace core;
using namespace core::import_pose;
using namespace core::io::pdb;
using namespace basic::options;
using namespace basic::options::OptionKeys;

using namespace core::pose;
using namespace core::chemical;
using namespace core::scoring;
using namespace core::kinematics;
using namespace protocols::moves;

// create a TaskFactory with the resfile
using namespace core::pack::task;
using namespace core::pack::dunbrack;
using namespace core::pack::dunbrack::cenrot;

static thread_local basic::Tracer TR( "pilot.wendao.cenrot" );

OPT_KEY(Boolean, debug_cenrot_min)

int main( int argc, char * argv [] ) {
	NEW_OPT(debug_cenrot_min, "debug cenrot min", false);

	devel::init(argc, argv);
	
	ResidueTypeSetCAP rsd_set;
	rsd_set=ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	protocols::simple_moves::SwitchResidueTypeSetMover to_fullatom("fa_standard");
	protocols::simple_moves::SwitchResidueTypeSetMover to_centroid("centroid");
	protocols::simple_moves::SwitchResidueTypeSetMover to_cenrot("centroid_rot");

	Pose p;
	//read in a fullatom pose
	core::import_pose::pose_from_pdb( p, *rsd_set, option[ in::file::native ]() );
	//convert it to cenrot
	to_cenrot.apply(p);
	core::scoring::ScoreFunctionOP score_fxn = core::scoring::ScoreFunctionFactory::create_score_function( "test" );

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if (!option[debug_cenrot_min]()) {
		// obsolete
		std::cout << "min_init" << std::endl;
		protocols::simple_moves::MinMover minmover;
		minmover.score_function(*score_fxn);
		minmover.min_type("lbfgs_armijo_nonmonotone");
		minmover.tolerance(1e-4);
		core::kinematics::MoveMapOP final_mm = new core::kinematics::MoveMap();
		//final_mm->set_chi( true );
		final_mm->set_bb( true );
		minmover.movemap(final_mm);

		std::cout << "min_apply" << std::endl;
		minmover.apply(p);
		
	}
	else {
		core::kinematics::MoveMap mm;
		mm.set_bb( true );
		//only for sidechain
		//mm.set_chi ( true );
		//mm.set( core::id::THETA, true );
		//mm.set( core::id::D, true );

		score_fxn->show(TR,p);
		TR.flush();

		core::optimization::MinimizerOptions options( "lbfgs_armijo_nonmonotone", 0.0001, true, true, true );
		core::optimization::CartesianMinimizer minimizer;
		std::cout << "CART MINTEST: " << "\n";
		long t1=clock();
		minimizer.run( p, mm, *score_fxn, options );
		long t2=clock();
		double time = ((double)t2 - t1) / CLOCKS_PER_SEC;
		std::cout << "end score: " << (*score_fxn)(p) << "\n";
		std::cout << "MIN TIME: " << time << " sec \n";
	}

	(*score_fxn)(p);
	score_fxn->show(TR,p);
	TR.flush();
	p.dump_pdb("minimized.pdb");
	return 0;
}
