// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/rms/rms_util.hh
/// @brief  RMS stuff from rosetta++
/// @author James Thompson
/// @author Ian Davis
/// @date   Wed Aug 22 12:10:37 2007
///

#ifndef core_scoring_rms_util_HH
#define core_scoring_rms_util_HH

#include <core/id/AtomID.fwd.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

#include <list>
#include <map>

//Auto Headers
#include <core/id/AtomID_Map.fwd.hh>
#include <utility/vector1_bool.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <iostream>

namespace core {
namespace scoring {

using namespace ObjexxFCL;

typedef std::list< core::Size > ResidueSelection;
typedef utility::vector1< core::Size > ResidueSelectionVector;

void invert_exclude_residues( core::Size nres, utility::vector1<int> const& exclude_list, ResidueSelection& );

ResidueSelection invert_exclude_residues( core::Size nres, utility::vector1<int> const& exclude_list );

extern core::Real native_CA_rmsd( const core::pose::Pose &native_pose ,  const core::pose::Pose &pose );

extern core::Real native_CA_gdtmm( const core::pose::Pose &native_pose ,  const core::pose::Pose &pose );




/// @brief RMSD between residues, accounting for automorphisms
/// (symmetries).  Does NOT include H atoms -- they add lots of extra
/// symmetries.
core::Real
automorphic_rmsd(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	bool superimpose
);


//////////////////////////////////////////////////////////////////////////////
// Predicate functions to use with rmsd_no_super() and rmsd_with_super()

bool
is_protein_CA(
	core::pose::Pose const & pose1,
	core::pose::Pose const & ,//pose2,
	core::Size resno,
	core::Size atomno
);

bool
is_protein_CA_or_CB(
	core::pose::Pose const & pose1,
	core::pose::Pose const & ,//pose2,
	core::Size resno,
	core::Size atomno
);

bool
is_protein_backbone(
	core::pose::Pose const& pose1,
	core::pose::Pose const&,
	core::Size resno,
	core::Size atomno
);

bool
is_protein_backbone_including_O(
	core::pose::Pose const& pose1,
	core::pose::Pose const&,
	core::Size resno,
	core::Size atomno
);

bool
is_protein_sidechain_heavyatom(
	core::pose::Pose const & pose1,
	core::pose::Pose const & ,//pose2,
	core::Size resno,
	core::Size atomno
);

bool
is_polymer_heavyatom(
	core::pose::Pose const & pose1,
	core::pose::Pose const & ,//pose2,
	core::Size resno,
	core::Size atomno
);

bool
is_ligand_heavyatom(
	core::pose::Pose const & pose1,
	core::pose::Pose const & ,//pose2,
	core::Size resno,
	core::Size atomno
);

bool
is_ligand_heavyatom_residues(
		core::conformation::Residue const & residue1,
		core::conformation::Residue const &, // residue2
		core::Size atomno
);

bool
is_heavyatom(
	core::pose::Pose const & pose1,
	core::pose::Pose const & ,//pose2,
	core::Size resno,
	core::Size atomno
);

bool
is_nbr_atom(
	core::pose::Pose const & pose1,
	core::pose::Pose const & ,//pose2,
	core::Size resno,
	core::Size atomno
);

//////////////////////////////////////////////////////////////////////////////

core::Real
CA_rmsd(
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2
);

// compute rmsd for residues between start and end
core::Real
CA_rmsd(
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2,
	Size start,
	Size end
);

core::Real
bb_rmsd(
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2
);

// compute rmsd for only backbone residues (excluding carboxyl oxygen) between start and end
core::Real
bb_rmsd(
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2,
	Size start,
	Size end,
	utility::vector1< Size > const& exclude
);

core::Real
bb_rmsd_including_O(
        const core::pose::Pose & pose1,
        const core::pose::Pose & pose2
);

// compute rmsd for only backbone residues (including carboxyl oxygen) between start and end
core::Real
bb_rmsd_including_O(
    const core::pose::Pose & pose1,
    const core::pose::Pose & pose2,
    Size start,
    Size end,
    utility::vector1< Size > const& exclude
);

// compute rmsd for residues between start and end
core::Real
CA_rmsd(
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2,
	Size start,
	Size end,
	utility::vector1< Size > const& exclude
);

// compute rmsd for residues for residues in list
core::Real
CA_rmsd(
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2,
	std::list<Size> residue_selection
);

core::Real
all_atom_rmsd(
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2
);

core::Real
all_atom_rmsd(
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2,
	std::list< Size > residue_selection
);

core::Real
nbr_atom_rmsd(
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2
);

void
fill_rmsd_coordinates(
	int & natoms,
	ObjexxFCL::FArray2D< core::Real > & p1a,
	ObjexxFCL::FArray2D< core::Real > & p2a,
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2,
	std::string atom_name
);

// other model-quality related functions

/// @brief Calculates a C-alpha maxsub-based superposition between pose1 and pose2, returns
/// the number of residues superimposed past a certain threshold. See maxsub.hh and maxsub.cc
/// for more information.
int
CA_maxsub(
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2,
	Real rms = 4.0
);

int
CA_maxsub(
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2,
	std::list<Size> residue_selection,
	Real rms = 4.0
);


int xyz_maxsub(
	FArray2D< core::Real > p1a,
	FArray2D< core::Real > p2a,
	int natoms
);

/// @brief Calculate gdtmm based on the given sets of xyz coordinates in p1a and p2a.
core::Real
xyz_gdtmm(
	FArray2D< core::Real > p1a,
	FArray2D< core::Real > p2a,
	core::Real& m_1_1,
	core::Real& m_2_2,
	core::Real& m_3_3,
	core::Real& m_4_3,
	core::Real& m_7_4
);

/// @brief Calculate gdtmm based on the given sets of xyz coordinates in p1a and p2a.
core::Real
xyz_gdtmm(
	FArray2D< core::Real > p1a,
	FArray2D< core::Real > p2a
);


/// @brief Calculate gdtmm score based on the C-alpha positions in pose1 and pose2.
core::Real
CA_gdtmm(
	core::pose::Pose const& pose1,
	core::pose::Pose const& pose2
);

/// @brief Calculate gdtmm score based on the C-alpha positions in pose1 and pose2. Also returns the
/// five components of the gdtmm score.
core::Real
CA_gdtmm(
	core::pose::Pose const& pose1,
	core::pose::Pose const& pose2,
	core::Real& m_1_1,
	core::Real& m_2_2,
	core::Real& m_3_3,
	core::Real& m_4_3,
	core::Real& m_7_4
);

core::Real
CA_gdtmm(
	core::pose::Pose const& pose1,
	core::pose::Pose const& pose2,
	std::list<Size> residue_selection,
	core::Real& m_1_1,
	core::Real& m_2_2,
	core::Real& m_3_3,
	core::Real& m_4_3,
	core::Real& m_7_4
);

/// @brief Calculate gdtmm score based on the C-alpha positions in pose1 and pose2.
core::Real
CA_gdtmm(
	core::pose::Pose const& pose1,
	core::pose::Pose const& pose2,
	std::list<Size> residue_selection
);

/// @brief  Superimpose mod_pose onto ref_pose using the mapping of atoms from
/// mod_pose to ref_pose given by atom_map
Real
superimpose_pose(
	pose::Pose & mod_pose,
	pose::Pose const & ref_pose,
	id::AtomID_Map< id::AtomID > const & atom_map // from mod_pose to ref_pose
);

/// @brief Superimpose two poses by their calpha coordinates.  Ignores residues
/// that do not have atoms named "CA."
Real
calpha_superimpose_pose(
	pose::Pose & mod_pose,
	pose::Pose const & ref_pose
);

core::Real
CA_rmsd_symmetric(
  const core::pose::Pose & pose1,
  const core::pose::Pose & pose2
);

core::Real
rms_at_corresponding_heavy_atoms(
	pose::Pose const & mod_pose,
	pose::Pose const & ref_pose );

void
create_shuffle_map_recursive_rms(
	std::vector<int> sequence,
	int const N,
	std::vector< std::vector<int> > & map
);

Real
rms_at_corresponding_atoms(
	pose::Pose const & mod_pose,
	pose::Pose const & ref_pose,
	std::map< core::id::AtomID, core::id::AtomID > atom_id_map
													 );

void
setup_matching_heavy_atoms( core::pose::Pose const & pose1, core::pose::Pose const & pose2, 	std::map< core::id::AtomID, core::id::AtomID > & atom_id_map );

Real
rms_at_corresponding_heavy_atoms(
	pose::Pose const & mod_pose,
	pose::Pose const & ref_pose );

/// @brief utility function to calculate per-residue sidechain rmsd without superposition
core::Real
residue_sc_rmsd_no_super(
	core::conformation::ResidueCOP res1,
	core::conformation::ResidueCOP res2,
	bool const fxnal_group_only=false );


} // end namespace scoring
} // end namespace core

#endif
