// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/energy_based_clustering/EnergyBasedClusteringProtocol.hh
/// @brief Performs the work done by the energy_based_clustering app.  Uses an energy-biased cookie-cutter approach to
/// cluster a large number of structures without generating an all-by-all RMSD matrix.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


#ifndef INCLUDED_protocols_cluster_energy_based_clustering_EnergyBasedClusteringProtocol_hh
#define INCLUDED_protocols_cluster_energy_based_clustering_EnergyBasedClusteringProtocol_hh

// Forward declarations include
#include <protocols/energy_based_clustering/EnergyBasedClusteringProtocol.fwd.hh>

// Associated headers
#include <protocols/energy_based_clustering/EnergyBasedClusteringOptions.hh>
#include <protocols/energy_based_clustering/EnergyBasedClusteringTests.fwd.hh>

// Core headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

// Protocols headers
#include <protocols/cyclic_peptide/CycpepSymmetryFilter.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.fwd.hh>

//Numeric headers
#include <numeric/xyzVector.fwd.hh>

namespace protocols {
namespace energy_based_clustering {

///@brief Performs the work done by the energy_based_clustering app.  Uses an energy-biased cookie-cutter approach to cluster a large number of structures without generating an all-by-all RMSD matrix.
class EnergyBasedClusteringProtocol : public utility::pointer::ReferenceCount {

	friend class ::EnergyBasedClusteringTests;

public:

	/// @brief Default constructor -- initializes from global options system.
	EnergyBasedClusteringProtocol();

	/// @brief Options constructor -- avoids access to global options system.
	EnergyBasedClusteringProtocol( EnergyBasedClusteringOptions const & options );

	virtual ~EnergyBasedClusteringProtocol();

	EnergyBasedClusteringProtocolOP
	clone() const;

	/// @brief Indicate which options-system options are relvant.
	static void register_options();

public:

	/// @brief Perform the clustering, based on options set in the options object.
	void go();

	/// @brief Get the number of clusters returned from the last call to go().
	/// @details Will be zero if go() has not been called.
	inline core::Size n_clusters_from_last_run() const { return n_clusters_from_last_run_; }

private: //Functions

	/// @brief Function to determine whether a value is in a list
	bool is_in_list ( core::Size const val, utility::vector1 < core::Size > const &vallist ) const;

	/// @brief Align one pose to another with an offset in the residue count.
	void align_with_offset (
		core::pose::Pose &pose1,
		core::pose::Pose const &pose2, //the target pose -- doesn't change.
		core::Size const offset,
		utility::vector1 <core::id::NamedAtomID> const &extra_atom_list
	) const;

	/// @brief Reconstruct a pose from posedata.
	/// @details The reconstruction_data vector is unused for anything except dihedral-based building whe rebuild_all_in_dihedral_mode is true.  In that case,
	/// this must be a full vector of coordinates of atoms.  Otherwise, it can be empty.
	void
	pose_from_posedata (
		core::pose::Pose const &inputpose,
		core::pose::Pose &outputpose,
		EBC_ClusterType const clustermode,
		utility::vector1<core::Real> const &posedata,
		utility::vector1<core::Real> const &reconstruction_data,
		bool const rebuild_all_in_dihedral_mode=false
	) const;

	/// @brief Sort the list of states in a cluster by energies.
	/// @details Inefficient selection sort used.  That's okay -- this sorts small lists.
	void sort_cluster_list( utility::vector1 <core::Size> &statelist, utility::vector1 <core::Real> &poseenergies ) const;

	/// @brief Function to calculate the RMSD between two poses, based on whatever is in "posedata" for the second only.
	/// This assumes that a pose already exists for the first, and that this is provided in refpose.
	core::Real
	calc_dist(
		utility::vector1 <core::Real> const &vect1, //Backbone dihedral vector 1
		utility::vector1 <core::Real> const &vect2, //Backbone dihedral vector 2
		EBC_ClusterType const clustmode,
		utility::vector1 < numeric::xyzVector < core::Real > > const &alignmentvect1, //Alignment vector 1 (for Cartesian clustering)
		utility::vector1 < numeric::xyzVector < core::Real > > const &alignmentvect2, //Alignment vector 2 (for Cartesian clustering)
		core::Size const nresidues,
		core::pose::Pose const &firstpose //Used for reference only!
	) const;

	/// @brief Overloaded form of calc_dist for cyclic permutations and for homooligomer swapping.
	core::Real
	calc_dist(
		utility::vector1 <core::Real> const &vect1,
		utility::vector1 <core::Real> const &vect2,
		EBC_ClusterType const clustmode,
		utility::vector1< numeric::xyzVector < core::Real > > const &alignmentvect1,
		utility::vector1< numeric::xyzVector < core::Real > > const &alignmentvect2,
		core::Size const nresidues,
		core::pose::Pose const &firstpose, //Used for reference only
		core::Size &offset, //Used only for cyclic permutations (OUTPUT)
		core::Size &permutation //Used only for homooligomer permutations (OUTPUT)
	) const;

	/// @brief Swap chains in an alignment vector (vector of x,y,z atom coordinates).
	/// @details Can only be called if options_.cluster_by_ == EBC_bb_cartesian.
	void swap_alignment_vector (
		utility::vector1 < numeric::xyzVector <core::Real > > const &parentvect,
		utility::vector1 < numeric::xyzVector <core::Real > >  &swappedvect,
		core::pose::Pose const &refpose,
		core::Size const permutation_number
	) const;

	/// @brief Returns the number of alignment atoms in a given residue.
	/// @details Some special-case logic here, since oxygens are counted in peptide-like building-blocks.
	core::Size alignment_atoms_in_res( core::pose::Pose const &pose, core::Size const position ) const;

	/// @brief Returns the number of alignment torsion angles stored for a given residue.
	core::Size alignment_torsions_in_res( core::pose::Pose const &pose, core::Size const position ) const;

	/// @brief Function containing the special-case logic for alpha, beta, and gamma-amino acids.
	bool use_this_atom( core::conformation::Residue const &res1, core::conformation::Residue const &res2, core::Size const atomno, std::string const &atname1, std::string const &atname2="" ) const;

	/// @brief Given a residue number, and atom number, and a pose, determine whether that atom gets used in the RMSD calculation.
	bool use_in_rmsd( core::pose::Pose const & pose1, core::pose::Pose const & pose2, core::Size resno, core::Size atomno, utility::vector1 <core::id::NamedAtomID> const &extra_atom_list ) const;

	/// @brief Given two poses, a residue number in the second pose, and atom number, and and offset for circular permutation, determine whether the atom is to be used in RMSD calculation.
	bool
	use_in_rmsd_offset(
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2,
		core::Size resno2,
		core::Size atomno,
		core::Size const pose1_offset,
		utility::vector1 <core::id::NamedAtomID> const &extra_atom_list
	) const;

	/// @brief Check that options have been set sensibly.
	void do_option_checks() const;

	/// @brief Get the global scorefunction, and set up constraints
	/// energies appropriately.
	core::scoring::ScoreFunctionOP set_up_scorefunction() const;

	/// @brief Given a vector of strings of format resnum:atomname, parse this out into a vector of NamedAtomIDs.
	void parse_extra_atom_list ( utility::vector1 < core::id::NamedAtomID > &extra_atom_list ) const;

	/// @brief Function to add user-specified constraints (specified with a CST file) to a pose.
	void add_user_constraints ( core::pose::Pose &mypose ) const;

	/// @brief Swap chains around in a pose (circular permutation of chain indices).
	void swap_chains(
		core::pose::Pose &currentpose, //Input and output
		core::Size const permutation_number //Input -- the chain perturbation number
	) const;

	/// @brief Get the start and end indices of chain chain_index, and store these in startaa and endaa.
	void
	get_start_and_end_aa(
		core::Size const chain_index,
		core::pose::Pose const &pose,
		core::Size & startaa,
		core::Size & endaa
	) const;

	/// @brief Get the start and end indices of chain chain_index, and store these in startaa and endaa.
	void
	get_start_and_end_indices(
		core::Size const chain_index,
		core::pose::Pose const &pose,
		core::Size &startindex,
		core::Size &endindex
	) const;

	/// @brief Function to form the disulfides based on user-specified disulfide positions (if provided), or based
	/// on automatic detection (if not).
	void make_disulfides ( core::pose::Pose &mypose ) const;

	/// @brief Function to mutate a pose to a chain of alanines.
	/// @details This does not mutate cysteine residues involved in disulfide bonds.  Note that this
	/// assumes that disulfide bonds have already been built in the pose.
	void mutate_to_alanine(core::pose::Pose &mypose) const;

	/// @brief Function to check whether two poses have matching classes of residues at matching positions:
	void check_backbones_match ( core::pose::Pose const &pose1, core::pose::Pose const &pose2 ) const;

	/// @brief Given an input pose, store only the relevant data needed for clustering.
	/// @details The "alignmentdata" array is only used for Cartesian clustering.
	void storeposedata(
		core::pose::Pose const &pose,
		utility::vector1 < core::Real > &posedata,
		utility::vector1< numeric::xyzVector< core::Real > > &alignmentdata, //Only for Cartesian clustering: x,y,z coordinates of atoms to be used for alignment.
		utility::vector1 < core::Real > &dihedral_mode_reconstruction_data, //Only for reconstructing pose in dihedral mode.
		EBC_ClusterType const clustermode,
		utility::vector1 < core::id::NamedAtomID > const &extra_atom_list
	) const;

	/// @brief Import all structures and score them, setting up derived data.
	/// @details Performs disk access.
	void do_initial_import_and_scoring( core::Size &count, core::Real &lowestE, core::Size &lowestE_index,
		core::pose::Pose &firstpose, protocols::cyclic_peptide::CycpepSymmetryFilter const &symmfilter,
		utility::vector1 <core::Real> &poseenergies, utility::vector1 < utility::vector1 <core::Real> > &posedata,
		utility::vector1< utility::vector1< numeric::xyzVector< core::Real > > > &alignmentdata,
		utility::vector1 < utility::vector1 <core::Real> > &dihedral_reconstruction_data,
		utility::vector1 < core::Size > &cluster_assignments, utility::vector1 < core::Size > &cluster_offsets,
		utility::vector1 < core::Size > &cluster_oligomer_permutations, core::scoring::ScoreFunctionOP sfxn,
		utility::vector1 < core::id::NamedAtomID > const &extra_atom_list
	) const;

	/// @brief Function to add cyclic constraints to a pose.
	void add_cyclic_constraints ( core::pose::Pose &mypose ) const;

private: //Data

	/// @brief Clustering options, either read from global options system or
	/// passed directly to the protocol.
	EnergyBasedClusteringOptions options_;

	/// @brief The contents of an optional constraints file.
	std::string constraints_file_contents_;

	/// @brief The number of clusters produced by the last run.
	core::Size n_clusters_from_last_run_;

};

} //energy_based_clustering
} //protocols

#endif //INCLUDED_protocols_cluster_energy_based_clustering_EnergyBasedClusteringProtocol_hh
