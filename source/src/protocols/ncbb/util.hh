// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


// Utilities used by various NCBB based design and dock-design applications.


// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Conformation.hh>

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
#include <protocols/simple_moves/BackboneMover.fwd.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/hbs/HbsPatcher.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>

// Filter headers
#include <basic/MetricValue.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
//#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>

#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/PackstatCalculator.hh>

// Utility Headers
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
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
using namespace protocols::simple_moves::hbs;
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

namespace protocols {
namespace ncbb {

void
count_uniq_char( std::string pattern, Size & num, utility::vector1<char> & uniqs );

Size
give_dihedral_index(
	Size n,
	utility::vector1< char > uniqs,
	std::string dihedral_pattern,
	std::string alpha_beta_pattern
);

void
ncbb_design_main_loop( Size loop_num, Size pert_num, Pose pose, TrialMoverOP pert_trial, utility::vector1<Size> designable_positions, Size pep_start, Size pep_end, TaskAwareMinMoverOP desn_ta_min, ScoreFunctionOP score_fxn, MonteCarloOP mc );

void
calculate_statistics( protocols::jd2::JobOP curr_job, core::pose::Pose pose, core::scoring::ScoreFunctionOP score_fxn  );

void
setup_pert_foldtree( core::pose::Pose & pose );

void
setup_filter_stats();

void
init_common_options( utility::tag::TagCOP tag, basic::datacache::DataMap &data, ScoreFunctionOP score_fxn_, Real & mc_temp_, Real & pert_mc_temp_, Real & pert_dock_rot_mag_, Real & pert_dock_trans_mag_, Real & pert_pep_small_temp_, Real & pert_pep_small_H_, Real & pert_pep_small_L_, Real & pert_pep_small_E_, Real & pert_pep_shear_temp_, Real & pert_pep_shear_H_, Real & pert_pep_shear_L_, Real & pert_pep_shear_E_, Size & pert_pep_num_rep_, Size & pert_num_, Size & dock_design_loop_num_, bool & no_design_, bool & final_design_min_, bool & use_soft_rep_, bool & mc_initial_pose_, bool & pymol_, bool & keep_history_ );

void
final_design_min( core::pose::Pose & pose, ScoreFunctionOP score_fxn_, core::pack::task::TaskFactoryOP desn_tf );

}
}
