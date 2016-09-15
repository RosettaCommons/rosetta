// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/rms/rms_util.hh
/// @brief  RMS stuff from rosetta++
/// @author James Thompson
/// @author Ian Davis
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_core_scoring_rms_util_HH
#define INCLUDED_core_scoring_rms_util_HH

// C/C++ headers
#include <list>
#include <map>
#include <string>

// External headers
#include <boost/unordered/unordered_map.hpp>

// Utility headers
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <utility>
#include <utility/vector1_bool.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>

// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/MiniPose.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {

typedef std::list< core::Size > ResidueSelection;
typedef utility::vector1< core::Size > ResidueSelectionVector;

void invert_exclude_residues( core::Size nres, utility::vector1<int> const& exclude_list, ResidueSelection& );

ResidueSelection invert_exclude_residues( core::Size nres, utility::vector1<int> const& exclude_list );

extern core::Real native_CA_rmsd( const core::pose::Pose &native_pose ,  const core::pose::Pose &pose );

extern core::Real native_CA_gdtmm( const core::pose::Pose &native_pose ,  const core::pose::Pose &pose );

/// @brief Returns a single, Global Distance Test-like value that measures the
/// extent to which the functional ends of a model's sidechains agree with their
/// counterparts in a given reference structure.
///
/// @detail Instead of comparing residue positions on the basis of CAs, gdtsc
/// uses a characteristic atom near the end of each sidechain type for the
/// evaluation of residue-residue distance deviations.
///
/// The traditional GDT score is a weighted sum of the fraction of residues
/// superimposed within limits of 1, 2, 4, and 8Å. For gdtsc, the backbone
/// superposition is used to calculate fractions of corresponding model-ref
/// sidechain atom pairs that fit under 10 distance-limit values from 0.5A
/// to 5A. Ambiguity in Asp or Glu terminal oxygen naming is not currently
/// considered.
///
/// Reference:
/// Keedy, DA. The other 90% of the protein. Proteins. 2009; 77 Suppl 9:29-49.
core::Real gdtsc(const core::pose::Pose& ref,
	const core::pose::Pose& model,
	const std::map<core::Size, core::Size>& residues);

/// @brief Returns the average fraction of residues superimposable under a
/// series of distance thresholds-- 0.5, 1.0, 2.0, and 4.0 Angstroms.
core::Real gdtha(const core::pose::Pose& ref,
	const core::pose::Pose& model,
	const std::map<core::Size, core::Size>& residues);

/// @brief Computes the RMSD of the jump residues between <model> and <native>,
/// storing the results in a map keyed by jump_id.
void compute_jump_rmsd(const core::pose::Pose& reference,
	const core::pose::Pose& model,
	boost::unordered_map<core::Size, core::Real>* rmsds);

/// @brief RMSD between residues, accounting for automorphisms
/// (symmetries).  For example if you have something like a tyrosine,
/// you won't get a higher rmsd just because you flipped the ring 180 degrees (Rocco).
/// Does NOT include H atoms -- they add lots of extra symmetries.
core::Real
automorphic_rmsd(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	bool superimpose
);

/// @brief  Compute the CA RMSD between two poses.
core::Real CA_rmsd(const core::pose::Pose& pose1,
	const core::pose::Pose& pose2,
	const std::map<core::Size, core::Size>& residues);

/// @brief  Compute the CA RMSD between two poses.
core::Real CA_gdtmm(const core::pose::Pose& pose1,
	const core::pose::Pose& pose2,
	const std::map<core::Size, core::Size>& residues);

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
is_protein_CA_or_equiv(
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

/// @brief Return true if the pose residues and atoms specified are non-peptide heavy atoms.
bool is_non_peptide_heavy_atom(
	core::pose::Pose const & pose1,
	core::pose::Pose const & /* pose2 */,
	core::uint const resno,
	core::uint const atomno );

bool
is_heavyatom(
	core::pose::Pose const & pose1,
	core::pose::Pose const & ,//pose2,
	core::Size resno,
	core::Size atomno
);

bool
is_scatom(
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

/////////////////////////////////////////////
// Predicate classes for more complex control

class Predicate: public utility::pointer::ReferenceCount {
public:
	Predicate() {};
	~Predicate() override = default;
	virtual bool operator()(
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2,
		core::Size resno,
		core::Size atomno) const = 0;
};

typedef utility::pointer::shared_ptr< Predicate > PredicateOP;
typedef utility::pointer::shared_ptr< Predicate const > PredicateCOP;

class IsProteinCAPredicate: public Predicate {
public:
	IsProteinCAPredicate() {}
	~IsProteinCAPredicate() override = default;
	bool operator()(
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2,
		core::Size resno,
		core::Size atomno) const override { return is_protein_CA(pose1, pose2, resno, atomno); }
};

class IsMainAtomPredicate: public Predicate {
public:
	IsMainAtomPredicate() {}
	~IsMainAtomPredicate() override = default;
	bool operator()(
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2,
		core::Size resno,
		core::Size atomno) const override { return is_protein_CA_or_equiv(pose1, pose2, resno, atomno); }
};

// (Fill in others as needed.)

//////////////////////
// Combining Predicates

class ResRangePredicate: public Predicate {
public:
	ResRangePredicate( core::Size start, core::Size end, PredicateCOP predicate ) : start_(start), end_(end), pred_(std::move(predicate)) {}
	~ResRangePredicate() override = default;
	bool operator()(
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2,
		core::Size resno,
		core::Size atomno) const override;
private:
	core::Size start_;
	core::Size end_;
	PredicateCOP pred_;
};

class SelectedResPredicate: public Predicate {
public:
	SelectedResPredicate( std::list< core::Size >  selected, PredicateCOP predicate ) : selected_(std::move(selected)), pred_(std::move(predicate)) {}
	~SelectedResPredicate() override = default;
	bool operator()(
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2,
		core::Size resno,
		core::Size atomno) const override;
private:
	std::list< core::Size > selected_;
	PredicateCOP pred_;
};

class ExcludedResPredicate: public Predicate {
public:
	ExcludedResPredicate( utility::vector1< Size > const & excluded, PredicateCOP predicate ) : excluded_(excluded), pred_(std::move(predicate)) {}
	~ExcludedResPredicate() override = default;
	bool operator()(
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2,
		core::Size resno,
		core::Size atomno) const override;
private:
	utility::vector1< Size > const & excluded_;
	PredicateCOP pred_;
};


//////////////////////////////////////////////////////////////////////////////

/// @brief Return the RMSD of the non-peptide heavy atoms of two poses.
core::DistanceSquared non_peptide_heavy_atom_RMSD( core::pose::Pose const & pose1, core::pose::Pose const & pose2 );

/// @brief Compute rmsd for residues between start and end.
/// If start and end aren't specified, use the entire pose.
core::Real
CA_rmsd(
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2,
	Size start = 1,
	Size end = 0
);

/// @brief Compute rmsd for residues between start and end.
/// If start and end aren't specified, use the entire pose.
core::Real
CA_or_equiv_rmsd(
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2,
	Size start = 1,
	Size end = 0
);

/// @brief Compute rmsd for only backbone residues (excluding carboxyl oxygen)
core::Real
bb_rmsd(
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2
);

/// @brief Compute rmsd for only backbone residues (including carboxyl oxygen)
core::Real
bb_rmsd_including_O(
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2
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
all_scatom_rmsd_nosuper(
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2
);

core::Real
all_atom_rmsd_nosuper(
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

/*
void
fill_rmsd_coordinates(
int & natoms,
ObjexxFCL::FArray2D< core::Real > & p1a,
ObjexxFCL::FArray2D< core::Real > & p2a,
const core::pose::Pose & pose1,
const core::pose::Pose & pose2,
std::string atom_name
); */

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
	ObjexxFCL::FArray2D< core::Real > p1a,
	ObjexxFCL::FArray2D< core::Real > p2a,
	int natoms
);

/// @brief Calculate gdtmm based on the given sets of xyz coordinates in p1a and p2a.
core::Real
xyz_gdtmm(
	ObjexxFCL::FArray2D< core::Real > p1a,
	ObjexxFCL::FArray2D< core::Real > p2a,
	core::Real& m_1_1,
	core::Real& m_2_2,
	core::Real& m_3_3,
	core::Real& m_4_3,
	core::Real& m_7_4
);

/// @brief Calculate gdtmm based on the given sets of xyz coordinates in p1a and p2a.
core::Real
xyz_gdtmm(
	ObjexxFCL::FArray2D< core::Real > p1a,
	ObjexxFCL::FArray2D< core::Real > p2a
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
/// mod_pose to ref_pose given by map< AtomID, AtomID >
/// @details The rms_calc_offset_val is a small constant value used by the numerical machinery to ensure a nonzero determinant.  This defaults to 1.0e-7.  Realign determines whether this is subtracted off again (default false).
Real
superimpose_pose(
	pose::Pose & mod_pose,
	pose::Pose const & ref_pose,
	std::map< id::AtomID, id::AtomID > const & atom_map, // from mod_pose to ref_pose
	core::Real const & rms_calc_offset_val = 1.0e-7,
	bool const realign=false
);

/// @brief Calculate gdttm score based on the C-alpha positions in pose1 and pose2.
void
CA_gdttm(
	core::pose::Pose const& pose1,
	core::pose::Pose const& pose2,
	core::Real &gdttm_score,
	core::Real &gdtha_score,
	std::list< Size > residue_selection //the std::list can be sorted! -- note std::sort can be applied to vectors
);

void
CA_gdttm(const core::pose::Pose& pose1,
	const core::pose::Pose& pose2,
	core::Real &gdttm_score,
	core::Real &gdtha_score,
	const std::map<core::Size, core::Size>& residues);

void
CA_gdttm(
	core::pose::Pose const& pose1,
	core::pose::Pose const& pose2,
	core::Real &gdttm_score,
	core::Real &gdtha_score
);

void
xyz_gdttm(
	ObjexxFCL::FArray2D< core::Real > p1a,
	ObjexxFCL::FArray2D< core::Real > p2a,
	core::Real &gdttm_score,
	core::Real &gdtha_score
);

/// @brief  Superimpose mod_pose onto ref_pose using the mapping of atoms from
/// mod_pose to ref_pose given by atom_map
/// @details The rms_calc_offset_val is a small constant value used by the numerical machinery to ensure a nonzero determinant.  This defaults to 1.0e-7. Realign determines whether this is subtracted off again (default false).
Real
superimpose_pose(
	pose::Pose & mod_pose,
	pose::Pose const & ref_pose,
	id::AtomID_Map< id::AtomID > const & atom_map, // from mod_pose to ref_pose
	core::Real const & rms_calc_offset_val = 1.0e-7,
	bool const realign=false
);

/// @brief  Superimpose mod_pose onto ref_pose using the mapping of atoms from
/// mod_pose to ref_pose given by atom_map
/// @details The rms_calc_offset_val is a small constant value used by the numerical machinery to ensure a nonzero determinant.  This defaults to 1.0e-7.  Realign determines whether this is subtracted off again (default false).
Real
superimpose_pose(
	pose::Pose & mod_pose,
	pose::MiniPose const & ref_pose,
	id::AtomID_Map< id::AtomID > const & atom_map, // from mod_pose to ref_pose
	core::Real const & rms_calc_offset_val = 1.0e-7,
	bool const realign=false
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
	std::map< core::id::AtomID, core::id::AtomID > const & atom_id_map
);

Real
rms_at_all_corresponding_atoms(
	pose::Pose const & mod_pose,
	pose::Pose const & ref_pose,
	std::map< core::id::AtomID, core::id::AtomID > const & atom_id_map
);

Real
rms_at_corresponding_atoms(
	pose::Pose const & mod_pose,
	pose::Pose const & ref_pose,
	std::map< core::id::AtomID, core::id::AtomID > const & atom_id_map,
	utility::vector1< Size > const & calc_rms_res
);

Real
rms_at_corresponding_atoms_no_super(
	pose::Pose const & mod_pose,
	pose::Pose const & ref_pose,
	std::map< core::id::AtomID, core::id::AtomID > const & atom_id_map );

Real
rms_at_corresponding_atoms_no_super(
	pose::Pose const & mod_pose,
	pose::Pose const & ref_pose,
	std::map< core::id::AtomID, core::id::AtomID > const & atom_id_map,
	utility::vector1< Size > const & calc_rms_res
);

void
setup_matching_heavy_atoms( core::pose::Pose const & pose1, core::pose::Pose const & pose2,  std::map< core::id::AtomID, core::id::AtomID > & atom_id_map );

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

void
setup_matching_CA_atoms( core::pose::Pose const & pose1, core::pose::Pose const & pose2,  std::map< core::id::AtomID, core::id::AtomID > & atom_id_map );

void
setup_matching_protein_backbone_heavy_atoms( core::pose::Pose const & pose1, core::pose::Pose const & pose2,
	std::map< core::id::AtomID, core::id::AtomID > & atom_id_map );

void
setup_matching_atoms_with_given_names( core::pose::Pose const & pose1, core::pose::Pose const & pose2,
	utility::vector1< std::string > const & atom_names_to_find,
	std::map< core::id::AtomID, core::id::AtomID > & atom_id_map );


} // end namespace scoring
} // end namespace core

#endif
