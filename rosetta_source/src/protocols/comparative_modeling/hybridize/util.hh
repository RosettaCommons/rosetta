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

#ifndef INCLUDED_protocols_comparative_modeling_hybridize_util_hh
#define INCLUDED_protocols_comparative_modeling_hybridize_util_hh

#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/Edge.hh>

#include <core/types.hh>

#include <list>

namespace protocols {
namespace comparative_modeling {
namespace hybridize {
			
using namespace core;
using namespace kinematics;

// constraint loading and generation
void setup_centroid_constraints( 
	core::pose::Pose &pose,
	utility::vector1 < core::pose::PoseCOP > templates,
	utility::vector1 < core::Real > template_weights,
	std::string cen_cst_file );

void setup_fullatom_constraints(
	core::pose::Pose &pose,
	utility::vector1 < core::pose::PoseCOP > templates,
	utility::vector1 < core::Real > template_weights,
	std::string cen_cst_file,
	std::string fa_cst_file  );

void generate_centroid_constraints( 
	core::pose::Pose &pose,
	utility::vector1 < core::pose::PoseCOP > templates,
	utility::vector1 < core::Real > template_weights );

void generate_fullatom_constraints(
	core::pose::Pose &pose,
	utility::vector1 < core::pose::PoseCOP > templates,
	utility::vector1 < core::Real > template_weights );


bool discontinued_upper(core::pose::Pose const & pose, Size const seqpos);

bool discontinued_lower(core::pose::Pose const & pose, Size const seqpos);

std::list < Size > downstream_residues_from_jump(core::pose::Pose const & pose, Size const jump_number);
	
// atom_map: from mod_pose to ref_pose
void
get_superposition_transformation(
								 pose::Pose const & mod_pose,
								 pose::Pose const & ref_pose,
								 id::AtomID_Map< id::AtomID > const & atom_map,
								 numeric::xyzMatrix< core::Real > &R, numeric::xyzVector< core::Real > &preT, numeric::xyzVector< core::Real > &postT );

void
partial_align(
			  core::pose::Pose & pose,
			  core::pose::Pose const & ref_pose,
			  id::AtomID_Map< id::AtomID > const & atom_map,
			  std::list <Size> const & residue_list,
			  bool iterate_convergence = false,
			  core::Real distance_squared_threshold = 4.0
			  );

core::id::AtomID_Map< core::id::AtomID >
update_atom_map(
				core::pose::Pose & pose,
				core::pose::Pose const & ref_pose,
				id::AtomID_Map< id::AtomID > const & atom_map,
				core::Real distance_squared_threshold
				);

Size
natom_aligned(
			  core::pose::Pose & pose,
			  core::pose::Pose const & ref_pose,
			  id::AtomID_Map< id::AtomID > const & atom_map,
			  core::Real distance_squared_threshold = 4.0
			  );

void
get_superposition_transformation(
								 pose::Pose const & mod_pose,
								 pose::Pose const & ref_pose,
								 core::id::AtomID_Map< core::id::AtomID > const & atom_map,
								 numeric::xyzMatrix< core::Real > &R, numeric::xyzVector< core::Real > &preT, numeric::xyzVector< core::Real > &postT );

void
apply_transformation(
				pose::Pose & mod_pose,
				std::list <Size> const & residue_list,
				numeric::xyzMatrix< core::Real > const & R, numeric::xyzVector< core::Real > const & preT, numeric::xyzVector< core::Real > const & postT
				);

} // hybridize 
} // comparative_modeling 
} // protocols

#endif

