#ifndef INCLUDED_protocols_glycopeptide_docking_utils_hh
#define INCLUDED_protocols_glycopeptide_docking_utils_hh

/// @file protocols/glycopeptide_docking/utils.cc
/// @brief util functions for peptide/glycopeptide sampling and refinement
/// @author Sai Pooja Mahajan (saipooja@gmail.com)
//
#include <protocols/glycopeptide_docking/GlycopeptideDockingFlags.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <map>

namespace protocols {
namespace glycopeptide_docking {
/// @brief Calculate distance between donor and substrate.
core::Real
calculate_sampled_distance(core::pose::Pose const &pose, core::Size const glycosylation_residue, core::Size const donor_residue);

/// @brief Calculate additional distance metrics for glycosylation

/// TODO: generalize
std::map<std::string,core::Real>
calculate_additional_glycosylation_metrics(core::pose::Pose const &pose, core::Size glycosylation_residue, core::Size const donor_residue );

/// @brief Write pdbs with output string with specified name and prefix, suffix and decoy number to match the fianl decoy.
void
write_debug_pdb(core::pose::Pose const &pose, core::Size const nstruct_max, core::Size const nstruct_index, std::string name);

/// @brief Glycosylate residues (sugar_residues) with
/// sugars specified in sugar_names. This feature is experimental.
void
glycosylate_residues(core::pose::Pose &pose,utility::vector1<core::Size> const &sugar_residues,utility::vector1<std::string> &sugar_names );

void
/// @details The foldtree for glycosylation can be set up in two ways.
/// The "outward" foldtree or the "docking" foldtree.\n The docking foldtree
/// is the standard docking foldtree and is setup across the enzyme-peptide or enzyme
/// -glycopeptide interface.\n
/// The outward foldtree is anchored at an "anchor residue"
/// specified with the anchor_residue flag. In the absence of this specification the
/// residue_to_glycosylate is also the anchor residue. This foldtree goes outwards in all
/// directions from the anchored position. The outward foldtree is
/// especially useful for glycopeptide glycosylation i.e. glycosylating a pre-glycosylated
/// peptide at a new site. In such a case, if the binding site of one sugar on the peptide
/// is known, the preglycosylated site on the peptide substrate can be specified as the
/// anchor residue. This allows for a floppy tail type of sampling.
/// Another useful choice for anchor residue is the site of glycosylation, especially when
/// the user is sampling with the site of glycosylation already glycosylated.
/// In the current implementation, a specialized foldtree is setup for a pose. User-specified foldtree is
/// not supported.
setup_glycosylation_foldtree( core::pose::Pose & pose,
	protocols::glycopeptide_docking::GlycopeptideDockingFlagsOP flags,
	core::kinematics::FoldTreeOP ft_docking);

/// @brief Record metrics such as important distances,
/// rmsds and interaction energies (to match publication)
void
record_pose_metrics( core::pose::Pose & pose,
	protocols::glycopeptide_docking::GlycopeptideDockingFlagsOP flags,
	utility::vector1< int > const jumps,
	core::pose::PoseOP ref_pose );

}
}
#endif //protocols_glycopeptide_docking_utils_HH
