// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/lddt.hh
/// @brief Utilities for calcuating lDDT.
///
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_scoring_lddt_hh
#define INCLUDED_core_scoring_lddt_hh

#include <core/pose/Pose.fwd.hh>

#include <core/conformation/Residue.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/id/AtomID.hh>

#include <core/types.hh>
#include <core/scoring/rms_util.hh> // For predicate

#include <utility/vector1.hh>

#include <map>


namespace core {
namespace scoring {

// Convenience functions

/// @brief Calculate the lDDT between two structures (1:1 residue correspondence)
/// across all heavy atoms.
/// If consider_alt is false, use the quicker algorithm which doesn't consider
/// automorphic changes (e.g. PHE ring flips)
core::Real
lddt(
	core::pose::Pose const & ref,
	core::pose::Pose const & model,
	bool consider_alt = true
);

/// @brief Calculate the per-residue lDDT between two structures (1:1 residue correspondence)
/// across all heavy atoms.
/// If consider_alt is false, use the quicker algorithm which doesn't consider
/// automorphic changes (e.g. PHE ring flips)
utility::vector1< core::Real >
per_res_lddt(
	core::pose::Pose const & ref,
	core::pose::Pose const & model,
	bool consider_alt = true
);

/// @brief Calculate the lDDT between the reference and the model.
/// https://dx.doi.org/10.1093/bioinformatics/btt473
/// Uses the same defaults as https://swissmodel.expasy.org/lddt/
/// Not implemented is any consideration for stereochemical/bond/angle deviation penalties
class lDDT_Calculator: public utility::VirtualBase {
public:

	/// @brief Create an lDDR_Calculator
	/// The predicate indicates which atoms to operate over.
	/// The default is all heavy atoms.
	lDDT_Calculator( PredicateCOP predicate = nullptr );

	/// @brief Set the predicate for which atoms to consider
	/// The default is all heavy atoms.
	/// Note the default of considering alternate chemical states is
	/// (currently) incompatible with including hydrogen atoms,
	/// so it is recommended to turn that off with any predicate
	/// which includes (potentially equivalent) hydrogens.
	void predicate( PredicateCOP predicate );

	/// @brief The distance similarity threasholds to use.
	/// The reported value will be an average of the lDDT for each individual threshold.
	void thresholds( utility::vector1< core::Real > const & setting ) { thresholds_ = setting; }

	/// @brief The inclusion radius: the cutoff to determine which distances are "local"
	void R0( core::Real val ) { R0_ = val; }

	/// @brief Only consider pairs if the polymeric sequence separation is seqsep or greater
	/// (different chains are infinite)
	/// This is different (one more) from the reference, to allow 0 to include distances within the same residue
	void seqsep( core::Size ss ) { seqsep_ = ss; }

	/// @brief The default is to ignore any interactions with the OXT atom.
	void ignore_oxt( bool setting ) { ignore_oxt_ = setting; }

	/// @brief When doing the calculation, use exact matching of atoms (false),
	/// Or consider chemically-equivaleint (automorphic) states to get the best value.
	/// Turning this off will speed calculations where alternate states aren't relevant
	/// (e.g. Calpha or backbone-only)
	void consider_alt_states( bool setting ) { consider_alt_ = setting; }

	//////// Calculators

	/// @brief Get the summary lDDT
	/// Assumes reference and model are the same length and match 1:1
	core::Real operator() (
		core::pose::Pose const & ref,
		core::pose::Pose const & model
	) const;

	/// @brief Get the summary lDDT, for the pairing specified in res_map.
	/// Only reference residues in res_map are considered for pairings.
	/// To include residues but set them unmatched, map them to zero.
	core::Real operator() (
		core::pose::Pose const & ref,
		core::pose::Pose const & model,
		std::map< core::Size, core::Size > const & res_map
	) const;

	/// @brief Get the residue lDDTs, for all atoms in the pose
	/// Assumes reference and model are the same length and match 1:1
	utility::vector1< core::Real > residue_lDDT(
		core::pose::Pose const & ref,
		core::pose::Pose const & model
	) const;

	/// @brief Get the residue lDDTs, for the pairing specified in res_map.
	/// Only reference residues in res_map are considered for pairings.
	/// To include residues but set them unmatched, map them to zero.
	/// The indicies for the returned residue_lddt are based on the *reference* structure.
	void residue_lDDT(
		core::pose::Pose const & ref,
		core::pose::Pose const & model,
		std::map< core::Size, core::Size > const & res_map,
		std::map< core::Size, core::Real > & residue_lddt // return by reference
	) const;

	/// @brief Get the atom lDDTs for all atoms in the poses
	/// Assumes reference and model are the same length and match 1:1
	/// The atom IDs for the returned value are based on the *reference* structure.
	std::map< core::id::AtomID, core::Real >
	atom_lDDT(
		core::pose::Pose const & ref,
		core::pose::Pose const & model
	) const;

	/// @brief Get the atom lDDTs, for the residue pairings specified in res_map.
	/// Only reference residues in res_map are considered for pairings.
	/// To include residues but set them unmatched, map them to zero.
	/// The atom IDs for the returned atom_lddt value are based on the *reference* structure.
	void atom_lDDT(
		core::pose::Pose const & ref,
		core::pose::Pose const & model,
		std::map< core::Size, core::Size > const & res_map,
		std::map< core::id::AtomID, core::Real > & atom_lddt // return by reference
	) const;

	/// @brief Get all lDDT statistics (global, residue and atom)
	/// Assumes reference and model are the same length and match 1:1
	/// The residue/atom IDs for the returned values are based on the *reference* structure.
	core::Real all_lDDT(
		core::pose::Pose const & ref,
		core::pose::Pose const & model,
		std::map< core::Size, core::Real > & residue_lddt, // return by reference
		std::map< core::id::AtomID, core::Real > & atom_lddt // return by reference
	) const;

	/// @brief Get all lDDT statistics (global, residue and atom)
	/// Only reference residues in res_map are considered for pairings.
	/// To include residues but set them unmatched, map them to zero.
	/// The residue/atom IDs for the returned values are based on the *reference* structure.
	core::Real all_lDDT(
		core::pose::Pose const & ref,
		core::pose::Pose const & model,
		std::map< core::Size, core::Size > const & res_map,
		std::map< core::Size, core::Real > & residue_lddt, // return by reference
		std::map< core::id::AtomID, core::Real > & atom_lddt // return by reference
	) const;

private:

	/// @brief A class to handle caching the calculations for the lddt, to speed things
	class lDDT_Cache {
	public:

		lDDT_Cache(PredicateCOP predicate, bool ignore_oxt) :
			predicate_( predicate ),
			ignore_oxt_( ignore_oxt )
		{}

		/// @brief Return (and cache) a list of the atom numbers for atoms matching the supplied predicate
		utility::vector1< core::Size > const &
		atoms_matching_predicate(
			core::pose::Pose const & pose,
			core::Size resno
		);

		/// @brief Return the number of states which the particular ResidueType has
		/// (Will return 1 for Residues without multiple mappings.)
		core::Size
		get_n_maps(
			core::chemical::ResidueType const & type
		);

		/// @brief Return (and cache) a specific chemically equivalent mapping
		/// (the equivalent mappings being a vector of atom index -> atom index)
		/// If the state_no is zero, or otherwise outside the range, will return the identity map.
		utility::vector1< core::Size > const &
		get_mapping_state(
			core::chemical::ResidueType const & type,
			core::Size state_no
		);

		/// @brief Return (and cache) a list of chemically equivalent mappings
		/// (the equivalent mappings being a vector of atom index -> atom index)
		utility::vector1< utility::vector1< core::Size > > const &
		residue_mapping_states(
			core::chemical::ResidueType const & type
		);

	private:

		/// @brief Generate the alternative mapping states for the particular ResidueType
		/// (TODO: Figure out how to auto-extract info about the predicate to inform state generation)
		void
		make_mapping_states(
			core::chemical::ResidueType const & type
		);

		/// @brief Utility function for make_mapping_states -- Flips the indexes for the specified atom on the passed vector
		/// Returns false on failure
		bool
		swap_indexes(
			utility::vector1< core::Size > & vec,
			core::chemical::ResidueType const & type,
			char const * const atm1, // String Literal
			char const * const atm2
		);

		/// @brief Utility function for make_mapping_states -- changes states to be from the automorphism iterator.
		void
		add_states_from_automorphisms(
			utility::vector1< utility::vector1< core::Size > > & states,
			core::chemical::ResidueType const & type
		);

	private:
		friend lDDT_Calculator;

		static utility::vector1< core::Size > const identity_map_;

	private:
		PredicateCOP predicate_;
		bool ignore_oxt_;

		// Data cache for atoms_matching_predicate()
		// Non-owning raw pointer, with limited lifetime
		std::map< core::conformation::Residue const *, utility::vector1< core::Size > > predicate_atoms_;

		// Data cache for a mapping of chemically equivalent states
		// Note that for space reasons, we don't consider hydrogen atoms.
		std::map< core::chemical::ResidueType const *, utility::vector1< utility::vector1< core::Size > > > mapping_states_;

	};

	/// @brief A class to store information about the various statistic counts.
	struct lDDT_Data {
		lDDT_Data( bool do_res = false, bool do_atom = false ):
			do_residue( do_res ),
			do_atomistic( do_atom )
		{}

		bool do_residue = false;
		bool do_atomistic = false;

		core::Size global_ndist = 0;
		core::Size global_nmatch = 0;

		// The residue-level statistics
		std::map< core::Size, core::Size > res_ndist;
		std::map< core::Size, core::Size > res_nmatch;

		// The atom-level statistics
		std::map< core::id::AtomID, core::Size > atom_ndist;
		std::map< core::id::AtomID, core::Size > atom_nmatch;

		// Safely calculate the global_lddt();
		core::Real global_lddt() const;

		void operator+=( lDDT_Data const & rhs );
		void operator-=( lDDT_Data const & rhs );

		/// @brief Helper function for operator+=
		template< typename T >
		static void add_inplace( std::map< T, core::Size> & lhs, std::map< T, core::Size > const & rhs );
		/// @brief Helper function for operator-=
		template< typename T >
		static void subtract_inplace( std::map< T, core::Size> & lhs, std::map< T, core::Size > const & rhs );
	};

	using lDDT_DataOP = utility::pointer::shared_ptr< lDDT_Data >;

	/// @brief Get the statistics, dispatching based on parameters.
	core::Real
	get_stats(
		core::pose::Pose const & ref,
		core::pose::Pose const & model,
		std::map< core::Size, core::Size > const & res_map,
		bool do_residue,
		std::map< core::Size, core::Real > & residue_stats, // return by reference
		bool do_atomistic,
		std::map< core::id::AtomID, core::Real > & atom_stats // return by reference
	) const;

	/// @brief Determine which states (from the cache) for each residue give the best overall lDDT.
	utility::vector1< core::Size >
	determine_alt_states(
		core::pose::Pose const & ref,
		core::pose::Pose const & model,
		std::map< core::Size, core::Size > const & res_map,
		lDDT_Cache & cache
	) const;

	core::Real
	get_stats_for_state(
		core::pose::Pose const & ref,
		core::pose::Pose const & model,
		std::map< core::Size, core::Size > const & res_map,
		utility::vector1< core::Size > const & state_assignment,
		lDDT_Cache & cache,
		lDDT_Data & data
	) const;

	void
	residue_pair_stats(
		core::pose::Pose const & ref,
		core::pose::Pose const & model,
		core::Size rres1_no,
		core::Size rres2_no,
		core::Size mres1_no,
		core::Size mres2_no,
		lDDT_Data & data,
		lDDT_Cache & cache,
		utility::vector1< core::Size > const & mres1_map = lDDT_Cache::identity_map_, // Mapping of residue atoms to alternate states.
		utility::vector1< core::Size > const & mres2_map = lDDT_Cache::identity_map_
	) const;

	///  @brief Turn the counts in lDDT_Data structure into fractions
	core::Real
	do_division(
		lDDT_Data & data,
		std::map< core::Size, core::Real > & residue_stats,
		std::map< core::id::AtomID, core::Real > & atom_stats
	) const;

	std::map< core::Size, core::Size >
	make_identity_map( core::pose::Pose const & pose ) const;

	///@brief Get the matching atom from two residues.
	/// Will return 0 if there is no matching atom.
	core::Size
	get_matching_atom(
		core::conformation::Residue const & rres,
		core::conformation::Residue const & mres,
		core::Size ratm,
		utility::vector1< core::Size > const & mres_map
	) const;

private:

	/// @brief The distance similarity threasholds to use.
	utility::vector1< core::Real > thresholds_ = { 0.5, 1, 2, 4 };

	/// @brief inclusion radius: the cutoff to determine which distances are "local"
	core::Real R0_ = 15.0;

	/// @brief Only consider pairs if the polymeric sequence separation is seqsep or greater
	/// (different chains are infinite)
	/// This is different (one more) from the reference, to allow 0 to include distances within the same residue
	core::Size seqsep_ = 1; // 1 is "everything except same residue"

	/// @brief Ignore OXT on terminal residues.
	bool ignore_oxt_ = true;

	/// @brief Consider automorphic residue states.
	bool consider_alt_ = true;

	/// @brief The predicate for which atoms to include
	PredicateCOP predicate_ = nullptr;
};

} // namespace scoring
} // namespace core

#endif
