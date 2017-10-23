// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/energy_based_clustering/EnergyBasedClusteringOptions.hh
/// @brief A container for the options used by the EnergyBasedClusteringProtocol.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


#ifndef INCLUDED_protocols_cluster_energy_based_clustering_EnergyBasedClusteringOptions_hh
#define INCLUDED_protocols_cluster_energy_based_clustering_EnergyBasedClusteringOptions_hh

#include <protocols/energy_based_clustering/EnergyBasedClusteringOptions.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// Core headers
#include <core/types.hh>

namespace protocols {
namespace energy_based_clustering {

/// @brief The type of clustering to perform.
enum EBC_ClusterType {
	EBC_bb_cartesian=1, //Keep first
	EBC_bb_dihedral, //Keep second-to-last
	EBC_end_of_list = EBC_bb_dihedral //Keep last
};

/// @brief A container for the options used by the EnergyBasedClusteringProtocol.
/// @details All data are public in this class.  It's just a dumb bag for holding data.
class EnergyBasedClusteringOptions : public utility::pointer::ReferenceCount {

public:

	/// @brief Default constructor.
	/// @details Initializes from options system by default.  Can be disabled
	/// by passing 'false'.
	EnergyBasedClusteringOptions( bool const initialize_from_options=true );

	/// @brief Destructor.
	virtual ~EnergyBasedClusteringOptions();

	EnergyBasedClusteringOptionsOP
	clone() const;

	/// @brief Initialize this option from the global options system.
	/// @details Called by default constructor.
	void initialize_from_global_options();

public:

	/// @brief Should imported structures be subjected to a round of fast relaxation?  Default false.
	bool prerelax_;

	/// @brief The number of fastrelax rounds to apply if the prerelax_ option is used.  Default 1.
	core::Size relax_rounds_;

	/// @brief What should I use as the basis for clustering?
	EBC_ClusterType cluster_by_;

	/// @brief If clustering by backbone Cartesian coordinates, should beta carbons be included?  Default false.
	/// Note that if this option is used, none of the input structures can contain glycine.
	bool use_CB_;

	/// @brief The radius for clustering, in Angstroms for Cartesian clustering and degrees for dihedral clustering.
	/// Default 1.0.
	core::Real cluster_radius_;

	// @brief Should structures be weighted by exp(-E/(k_B*T)) when calculating cluster centers?  True by default.
	//bool weight_by_energy_;

	// @brief The value of k_b*T to use when weighting by exp(-E/(k_B*T)).  0.62 by default.
	//core::Real kbt_;

	/// @brief List of residues to ignore in alignments for clustering.  Default empty list.
	utility::vector1< core::Size > residues_to_ignore_;

	/// @brief List of chains to ignore in alignments for clustering.  Default empty list.
	utility::vector1< core::Size > chains_to_ignore_;

	/// @brief Maximum number of structures to output per cluster.  Default no limit (0).
	core::Size limit_structures_per_cluster_;

	/// @brief Maximum number of clusters to output.  Default no limit (0).
	core::Size limit_clusters_;

	/// @brief If true, constraints are added to make a peptide bond between the N- and C-termini.  If false (default),
	/// the termini are free.  Default false.
	bool cyclic_;

	/// @brief If provided, structures that do not have the desired symmetry are filtered out.  Set to 2 for C2
	/// or S2 symmetry, 3 for C3 symmetry, 4 for C4 or S4 symmetry, etc.  Unused (0) if not specified.  Can only
	/// be used with the cyclic_ option.
	core::Size cyclic_symmetry_;

	/// @brief If true, then SN symmetry is used instead of CN.  Unused if not specified.  Can only be used with
	/// the cyclic_ and cyclic_symmetry_ options.
	bool cyclic_symmetry_mirroring_;

	/// @brief The angle threshold, in degrees, for determining whether a cyclic peptide is symmetric.  Can only
	/// be used with the cyclic_ and cyclic_symmetry_ flags.  Defaults to 10.0 degrees.
	core::Real cyclic_symmetry_threshold_;

	/// @brief Has the user specified a cyclic symmetry threshold?
	bool cyclic_symmetry_threshold_specified_;

	/// @brief If true, all cyclic permutations are tried when comparing two structures for clustering.  Requires
	/// cyclic_ true.  Default false.
	bool cluster_cyclic_permutations_;

	/// @brief 1 by default, meaning that every cyclic permutation is clustered if cluster_cyclic_permutations_
	/// is true.  Values X > 1 mean that cyclic permutations shifted by X residues will be clustered.
	core::Size cyclic_permutation_offset_;

	/// @brief If true, the input structures will be converted to a chain of alanines (L- or D-) before scoring.
	/// Default false.
	bool mutate_to_ala_;

	/// @brief A space-separated list of positions that are disulfide-bonded.  For example,
	/// disulfide_positions_ 3 8 6 23 would mean that residues 3 and 8 are disulfide-bonded,
	/// as are residues 6 and 23.  Default empty list.
	utility::vector1< core::Size > disulfide_positions_;

	/// @brief If the structures contain multiple chains with identical sequence, setting this to true
	/// will test all permutations of chains when clustering.  Default false.
	bool homooligomer_swap_;

	/// @brief Write output to a silent file instead of to separate PDBs.  This will create two files:
	/// one that only contains the first member of each cluster, and one that contains everything.
	/// Default false.
	bool silent_output_;

	/// @brief A user-specified set of one or more constraints file.  Default unused.
	utility::vector1< std::string > cst_files_;

	/// @brief A list of additional atoms to use in the RMSD calculation,
	/// each in the format residue:atomname separated by whitespace.  For
	/// example, -v_extra_rms_atoms 7:SG 12:CG 12:CD 12:CE 12:NZ 14:OG.
	/// Default empty list.
	utility::vector1< std::string > extra_rms_atoms_;

	/// @brief In dihedral clustering mode, do we rebuild everything for output?
	/// @details Just rebuilds the mainchain if false.  True by default.
	bool rebuild_all_in_dihedral_mode_;

	/// @brief Prefix to prepend on output files.
	/// @details By default, read from options system.
	std::string output_prefix_;


};

} //energy_based_clustering
} //protocols

#endif //INCLUDED_protocols_cluster_energy_based_clustering_EnergyBasedClusteringOptions_hh
