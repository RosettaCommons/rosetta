// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file fa_to_cenrot_score.cc
/// @brief test centroid rot model
/// @author Yuan Liu


// Core Headers
#include <core/chemical/AA.hh>
#include <core/chemical/Atom.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>

#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/PDBInfo.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/TenANeighborGraph.hh>

//
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/dunbrack/cenrot/SingleResidueCenrotLibrary.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <numeric/constants.hh>
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
//#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/Minimizer.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <core/scoring/dssp/Dssp.hh>

#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/moves/Mover.fwd.hh>

//sampling
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/BBGaussianMover.hh>

#include <sstream>
#include <iostream>
#include <string>
#include <fstream>
#include <map>

/////////////////////////////////////////////////////////////////
using namespace core;
using namespace core::import_pose;
using namespace core::io::pdb;
using namespace basic::options;
using namespace basic::options::OptionKeys;

using namespace core::pose;
using namespace core::chemical;
using namespace core::conformation;
using namespace core::scoring;
using namespace core::kinematics;
using namespace protocols::moves;

using namespace core::pack::task;
using namespace core::pack::dunbrack;
using namespace core::pack::dunbrack::cenrot;

//////////////////////////////////////////////////////////////////
static numeric::random::RandomGenerator RG(62331911);
basic::Tracer TR("pilot.wendao.cenrot");

void relax_cenrot_pose(core::pose::PoseOP &native_pose,
	core::pose::Pose & p, std::string const &tag);

OPT_KEY(String, cenrot_score)
OPT_KEY(Real, relax_temp)
OPT_KEY(Integer, relax_step_per_cycle)
OPT_KEY(Integer, relax_cycle_number)

int main( int argc, char * argv [] ) {
	NEW_OPT(cenrot_score, "cenrot score weight file", "test.wts");

	NEW_OPT(relax_temp, "temp", 1.0);
	NEW_OPT(relax_step_per_cycle, "step", 100);
	NEW_OPT(relax_cycle_number, "cycle", 1);

	devel::init(argc, argv);

	ResidueTypeSetCAP rsd_set;
	rsd_set=ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	//SwitchResidueTypeSetCenrot to_cenrot;
	protocols::simple_moves::SwitchResidueTypeSetMover to_cenrot("centroid_rot");

	//load native pdb
	PoseOP native_pose;
	if (option[in::file::native].user()) {
		native_pose = new Pose();
		core::import_pose::pose_from_pdb( *native_pose, *rsd_set, option[ in::file::native ]() );
	}

	if (option[ in::file::s ].user()) {
		Size npdbs = option[ in::file::s ]().size();
		for (Size npdb=1; npdb<=npdbs; npdb++)
		{
			PoseOP pose = new Pose();
			Pose &p(*pose);
			pose_from_pdb( p, *rsd_set, option[ in::file::s ]()[npdb] );

			to_cenrot.apply(p); //switch_to_residue_type_set_cenrot(p);

			relax_cenrot_pose(native_pose, p, option[ in::file::s ]()[npdb]);
		}
	}
	return 0;
}

void relax_cenrot_pose(core::pose::PoseOP &native_pose, core::pose::Pose & p, std::string const &tag)
{
	using namespace core::pack::task;
	using namespace protocols;

	pose::Pose native_p(p);

	//score function
	core::scoring::ScoreFunctionOP score_bb_fxn;
	core::scoring::ScoreFunctionOP score_sc_fxn;
	score_bb_fxn = core::scoring::ScoreFunctionFactory::create_score_function(option[cenrot_score]);
	score_sc_fxn = core::scoring::ScoreFunctionFactory::create_score_function(option[cenrot_score]);

	//repack mover
	TaskFactoryOP main_task_factory = new TaskFactory;
	operation::RestrictToRepackingOP rtrop = new operation::RestrictToRepacking;
	main_task_factory->push_back( rtrop );
	protocols::simple_moves::PackRotamersMoverOP packrotamersmover(new protocols::simple_moves::PackRotamersMover());
	packrotamersmover->task_factory(main_task_factory);
	packrotamersmover->score_function(score_sc_fxn);
	//protocols::simple_moves::RotamerTrialsMoverOP packrotamersmover(new protocols::simple_moves::RotamerTrialsMover());
	//packrotamersmover->task_factory(main_task_factory);
	//packrotamersmover->score_function(score_sc_fxn);

	//gaussian mover
	simple_moves::BBG8T3AMoverOP bbgmover(new protocols::simple_moves::BBG8T3AMover());

	//monte carlo
	MonteCarloOP mc = new MonteCarlo(p, *score_bb_fxn, option[relax_temp]);

	//combo
	moves::SequenceMoverOP combo( new moves::SequenceMover() );
	combo->add_mover(bbgmover);
	combo->add_mover(packrotamersmover);
	moves::TrialMoverOP trial ( new moves::TrialMover(combo, mc) );
	moves::RepeatMoverOP run( new moves::RepeatMover(trial, option[relax_step_per_cycle]) );

	//this is not real relax
	for (Size i=1; i<=option[relax_cycle_number]; i++) {
		run->apply(p);

		score_bb_fxn->show(TR,p);
		TR << "score: " << (*score_bb_fxn)(p) << " rmsd: " << core::scoring::CA_rmsd(p, native_p) << std::endl;
		TR.flush();

		std::ostringstream outputfn;
		mc->show_counters();
		mc->reset_counters();

		outputfn << "traj_" << i << ".pdb";
		p.dump_pdb(outputfn.str());
	}
}

