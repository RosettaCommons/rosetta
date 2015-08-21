/*
* FoldFromLoops_functions.hh
*
*  Created on: May 1, 2009
*      Author: bcorreia
*/

#ifndef FOLDFROMLOOPS_FUNCTIONS_HH_
#define FOLDFROMLOOPS_FUNCTIONS_HH_

#include <core/types.hh>


#include <core/fragment/FragSet.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <utility/vector1.hh>


namespace devel {
namespace fold_from_loops {


bool is_loop (
	protocols::loops::Loops & loops,
	core::Size & residue
);


std::vector<core::Size> define_cut_points (

	protocols::loops::Loops & loops,
	core::pose::Pose & nat_pose
);


void fold_tree_cutpoints_generator(
	protocols::loops::Loops & loops ,
	std::vector<core::Size> & cutpoints,
	core::pose::Pose & pose,
	core::kinematics::FoldTree & f
);


void CA_cst_generator(core::pose::Pose & pose,
	core::scoring::constraints::ConstraintSetOP & cst,
	protocols::loops::Loops & loops,
	std::vector<core::Size> & cut_points
);

void CA_cst_generator(core::pose::Pose & pose,
	core::scoring::constraints::ConstraintSetOP & cst,
	protocols::loops::Loops & loops
);


bool is_loop_neighbor( protocols::loops::Loops & loops,
	core::Size & residue,
	core::Size & range
);


void define_movemap_extending_chain(
	core::kinematics::MoveMapOP & movemap,
	core::pose::Pose & pose,
	protocols::loops::Loops & loops
);


void define_movemap(
	core::kinematics::MoveMapOP & movemap,
	core::pose::Pose & pose,
	protocols::loops::Loops & loops
);

void extending_chain(
	core::kinematics::MoveMapOP & movemap,
	core::pose::Pose & pose
);

void get_fragments(
	core::fragment::FragSetOP & fragset_large_,
	core::fragment::FragSetOP & fragset_small_
);


void new_pose_generator(
	core::pose::Pose & target_loops,
	core::pose::Pose & nat_prot,
	protocols::loops::Loops & loops,
	std::vector<core::Size> & cut_points
);


void refresh_cutpoints(
	core::pose::Pose & pose,
	std::vector<core::Size> & cut_points
);


void copying_side_chains_swap_loop (
	core::pose::Pose & swap_loops,
	core::pose::Pose & fold_pose,
	protocols::loops::Loops & loops,
	core::kinematics::MoveMapOP & movemap
);


void design_excluding_swap_loops (
	core::pose::Pose & fold_pose,
	protocols::loops::Loops & loops_in,
	core::scoring::ScoreFunctionOP & scorefxn_fa

);


void exclude_loop_residues( core::pose::Pose & pose,
	utility::vector1< bool > & residues_to_mutate,
	utility::vector1< bool > & allowed_aas ,
	core::pack::task::PackerTaskOP & task,
	protocols::loops::Loops & loops
);

bool is_cut( std::vector<core::Size> & cut_points,
	core::Size & residue
);


void copying_side_chains(
	core::pose::Pose & nat_pose,
	core::pose::Pose & fold_pose,
	protocols::loops::Loops & loops,
	core::kinematics::MoveMapOP & movemap

);


}
}

#endif /* FOLDFROMLOOPS_FUNCTIONS_HH_ */
