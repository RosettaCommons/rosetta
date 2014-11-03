// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Align a random jump to template
/// @detailed
/// @author Yifan Song

#ifndef INCLUDED_protocols_hybridization_util_hh
#define INCLUDED_protocols_hybridization_util_hh

#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/conformation/Residue.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/Edge.hh>

#include <core/types.hh>

#include <list>

#include <protocols/loops/Loops.hh>
#include <protocols/loops/Loop.hh>
#include <core/sequence/SequenceAlignment.hh>

#include <set>

namespace protocols {
//namespace comparative_modeling {
namespace hybridization {

core::Size get_num_residues_nonvirt( core::pose::Pose const & pose );
core::Size get_num_residues_prot( core::pose::Pose const & pose );

// constraint loading and generation
void setup_centroid_constraints(
	core::pose::Pose &pose,
	utility::vector1 < core::pose::PoseCOP > templates,
	utility::vector1 < core::Real > template_weights,
	std::string cen_cst_file,
	std::set< core::Size > ignore_res_for_AUTO = std::set<core::Size>());

void setup_fullatom_constraints(
	core::pose::Pose &pose,
	utility::vector1 < core::pose::PoseCOP > templates,
	utility::vector1 < core::Real > template_weights,
	std::string cen_cst_file,
	std::string fa_cst_file  );

void generate_centroid_constraints(
	core::pose::Pose &pose,
	utility::vector1 < core::pose::PoseCOP > templates,
	utility::vector1 < core::Real > template_weights,
	std::set< core::Size > ignore_res = std::set<core::Size>());

void generate_fullatom_constraints(
	core::pose::Pose &pose,
	utility::vector1 < core::pose::PoseCOP > templates,
	utility::vector1 < core::Real > template_weights );

// user-defined coord csts
void setup_user_coordinate_constraints( core::pose::Pose &pose, utility::vector1<core::Size> reses_to_cst );

// restrict all but interface
void setup_interface_coordinate_constraints( core::pose::Pose &pose, utility::vector1<bool> ignore_res );

// restrict all but interface
void setup_interface_atompair_constraints( core::pose::Pose &pose, utility::vector1<bool> ignore_res );

// input strand pairings
void add_strand_pairs_cst(core::pose::Pose & pose, utility::vector1< std::pair< core::Size, core::Size > > const strand_pairs);

// ligand/DNA
void add_non_protein_cst(core::pose::Pose & pose, core::pose::Pose & tmpl, core::Real const self_cst_weight, core::Real const het_prot_cst_weight);

bool discontinued_upper(core::pose::Pose const & pose, core::Size const seqpos);

bool discontinued_lower(core::pose::Pose const & pose, core::Size const seqpos);

std::list < core::Size > downstream_residues_from_jump(core::pose::Pose const & pose, core::Size const jump_number);

// atom_map: from mod_pose to ref_pose
void
get_superposition_transformation(
								 core::pose::Pose const & mod_pose,
								 core::pose::Pose const & ref_pose,
								 core::id::AtomID_Map< core::id::AtomID > const & atom_map,
								 numeric::xyzMatrix< core::Real > &R, numeric::xyzVector< core::Real > &preT, numeric::xyzVector< core::Real > &postT );

core::Size atom_map_valid_size(
							   core::pose::Pose const & pose,
							   core::id::AtomID_Map< core::id::AtomID > const & atom_map
							   );

void
partial_align(
		core::pose::Pose & pose,
		core::pose::Pose const & ref_pose,
		core::id::AtomID_Map< core::id::AtomID > const & atom_map,
		bool iterate_convergence,
		utility::vector1<core::Real> distance_thresholds,
		core::Real min_coverage );

void
partial_align(
			  core::pose::Pose & pose,
			  core::pose::Pose const & ref_pose,
			  core::id::AtomID_Map< core::id::AtomID > const & atom_map,
			  std::list <core::Size> const & residue_list,
			  bool iterate_convergence = false,
				utility::vector1<core::Real> distance_thresholds=utility::vector1<core::Real>(0),
				core::Real min_coverage = 0.2);

core::id::AtomID_Map< core::id::AtomID >
update_atom_map(
				core::pose::Pose & pose,
				core::pose::Pose const & ref_pose,
				core::id::AtomID_Map< core::id::AtomID > const & atom_map,
				core::Real distance_squared_threshold
				);

core::Size
natom_aligned(
			  core::pose::Pose & pose,
			  core::pose::Pose const & ref_pose,
			  core::id::AtomID_Map< core::id::AtomID > const & atom_map,
			  core::Real distance_squared_threshold = 4.0
			  );

void
get_superposition_transformation(
								 core::pose::Pose const & mod_pose,
								 core::pose::Pose const & ref_pose,
								 core::id::AtomID_Map< core::id::AtomID > const & atom_map,
								 numeric::xyzMatrix< core::Real > &R, numeric::xyzVector< core::Real > &preT, numeric::xyzVector< core::Real > &postT );

void
apply_transformation(
					 core::pose::Pose & mod_pose,
					 std::list <core::Size> const & residue_list,
					 numeric::xyzMatrix< core::Real > const & R, numeric::xyzVector< core::Real > const & preT, numeric::xyzVector< core::Real > const & postT
					 );

core::fragment::FragSetOP
create_fragment_set( core::pose::Pose const & pose, core::Size len, core::Size nfrag );

core::fragment::FragSetOP
create_fragment_set_no_ssbias( core::pose::Pose const & pose, core::Size len, core::Size nfrag, char force_ss='D' );

core::fragment::FragSetOP
create_fragment_set_no_ssbias( core::pose::Pose const & pose, std::set<core::Size> user_pos, core::Size len, core::Size nfrag, char force_ss='D' );



protocols::loops::Loops
renumber_with_pdb_info( protocols::loops::Loops & template_chunk, core::pose::PoseCOP template_pose );

core::Real get_gdtmm( core::pose::Pose & native, core::pose::Pose & pose, core::sequence::SequenceAlignmentOP & aln );

} // hybridize
//} // comparative_modeling
} // protocols

#endif

