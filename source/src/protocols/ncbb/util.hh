// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


// Utilities used by various NCBB based design and dock-design applications.


// Project Headers
#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

#include <protocols/jd2/Job.fwd.hh>

// Mover headers
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/moves/TrialMover.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/simple_moves/TaskAwareMinMover.fwd.hh>

// Utility Headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/Tracer.hh>
// C++ headers
#include <string>
#include <sstream>

namespace protocols {
namespace ncbb {

void
count_uniq_char( std::string pattern, core::Size & num, utility::vector1<char> & uniqs );

core::Size
give_dihedral_index(
	core::Size n,
	utility::vector1< char > uniqs,
	std::string dihedral_pattern,
	std::string alpha_beta_pattern
);

core::Size
get_number_dihedrals(
	utility::vector1< char > uniqs,
	std::string const & dihedral_pattern,
	std::string const & alpha_beta_pattern
);

void
ncbb_design_main_loop(
	core::Size loop_num,
	core::Size pert_num,
	core::pose::Pose pose,
	protocols::moves::TrialMoverOP pert_trial,
	utility::vector1<core::Size> designable_positions,
	core::Size pep_start,
	core::Size pep_end,
	protocols::simple_moves::TaskAwareMinMoverOP desn_ta_min,
	core::scoring::ScoreFunctionOP score_fxn,
	protocols::moves::MonteCarloOP mc
);

void
calculate_statistics( protocols::jd2::JobOP curr_job, core::pose::Pose pose, core::scoring::ScoreFunctionOP score_fxn  );

void
setup_pert_foldtree( core::pose::Pose & pose );

void
setup_filter_stats();

void
init_common_options(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &data,
	core::scoring::ScoreFunctionOP score_fxn_,
	core::Real & mc_temp_,
	core::Real & pert_mc_temp_,
	core::Real & pert_dock_rot_mag_,
	core::Real & pert_dock_trans_mag_,
	core::Real & pert_pep_small_temp_,
	core::Real & pert_pep_small_H_,
	core::Real & pert_pep_small_L_,
	core::Real & pert_pep_small_E_,
	core::Real & pert_pep_shear_temp_,
	core::Real & pert_pep_shear_H_,
	core::Real & pert_pep_shear_L_,
	core::Real & pert_pep_shear_E_,
	core::Size & pert_pep_num_rep_,
	core::Size & pert_num_,
	core::Size & dock_design_loop_num_,
	bool & no_design_,
	bool & final_design_min_,
	bool & use_soft_rep_,
	bool & mc_initial_pose_,
	bool & pymol_,
	bool & keep_history_
);

void
final_design_min( core::pose::Pose & pose, core::scoring::ScoreFunctionOP score_fxn_, core::pack::task::TaskFactoryOP desn_tf );

}
}
