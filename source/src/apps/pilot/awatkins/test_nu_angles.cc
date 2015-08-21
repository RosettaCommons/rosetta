// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/ScoringManager.hh>

#include <core/chemical/VariantType.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <utility/pointer/owning_ptr.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Mover headers
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/ncbb/util.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/ResidueFactory.hh>

// Filter headers
#include <basic/MetricValue.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
//#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>

#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/PackstatCalculator.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>

// C++ headers
#include <string>
#include <sstream>

//The original author used a lot of using declarations here.  This is a stylistic choice.
// Namespaces
using namespace core;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace pose;
using namespace protocols;
using namespace protocols::moves;
using namespace protocols::simple_moves;
using namespace protocols::rigid;
using namespace protocols::toolbox;
using namespace protocols::toolbox::pose_metric_calculators;
using namespace core::pack::task;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::id;
using basic::T;
using basic::Error;
using basic::Warning;
using utility::file::FileName;

// tracer - used to replace cout
static basic::Tracer TR("test_nu_angles" );

int
main( int argc, char* argv[] )
{
	try {

		//utility::vector1< core::Size > empty_vector(0);

		// init command line options
		devel::init(argc, argv);

		Pose pose;
		ScoreFunctionOP score_fxn = get_score_function();
		score_fxn->set_weight( ring_close, 1 );
		ResidueTypeSetCOP rts( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD ) );

		pose.append_residue_by_jump( *new Residue( rts->name_map("A04"), true ), 1 );
		pose.append_residue_by_bond( *new Residue( rts->name_map("A98"), true ), true );
		pose.append_residue_by_bond( *new Residue( rts->name_map("B02"), true ), true );
		pose.append_residue_by_bond( *new Residue( rts->name_map("B06"), true ), true );
		pose.append_residue_by_bond( *new Residue( rts->name_map("C01"), true ), true );
		pose.append_residue_by_bond( *new Residue( rts->name_map("B19"), true ), true );
		pose.append_residue_by_bond( *new Residue( rts->name_map("B95"), true ), true );
		pose.append_residue_by_bond( *new Residue( rts->name_map("C00"), true ), true );
		pose.append_residue_by_bond( *new Residue( rts->name_map("C12"), true ), true );

		for ( Size i = 1; i <= pose.total_residue(); ++i ) {
			pose.set_phi( i, -150 );
			pose.set_psi( i,  150 );
			pose.set_omega( i,  180 );
		}
		pose.dump_scored_pdb( "out.pdb", *score_fxn );

		kinematics::MoveMapOP mm( new kinematics::MoveMap );
		mm->set_bb( true );
		mm->set_chi( true );
		mm->set_nu( true );

		MinMoverOP min( new MinMover( mm, score_fxn, "lbfgs_armijo_nonmonotone", 0.0001, true ) );
		min->apply( pose );

		pose.dump_scored_pdb( "min.pdb", *score_fxn );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}//main

