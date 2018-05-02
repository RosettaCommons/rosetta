// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/guidance_scoreterms/buried_unsat_penalty/graph/BuriedUnsatPenaltyGraph.hh
/// @brief Headers for the BuriedUnsatPenaltyGraph class and its related Node and Edge classes.
/// @author Vikram K. Mulligan (vmullig@uw.edu).

#ifndef INCLUDED_core_pack_guidance_scoreterms_buried_unsat_penalty_graph_BuriedUnsatPenaltyGraph_HH
#define INCLUDED_core_pack_guidance_scoreterms_buried_unsat_penalty_graph_BuriedUnsatPenaltyGraph_HH

#include <core/pack/guidance_scoreterms/buried_unsat_penalty/graph/BuriedUnsatPenaltyGraph.fwd.hh>
#include <core/pack/guidance_scoreterms/buried_unsat_penalty/graph/BuriedUnsatPenaltyGraphOptions.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/hbonds/HBondOptions.fwd.hh>
#include <core/scoring/hbonds/HBondDatabase.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>

#include <numeric/xyzVector.hh>

#include <utility/graph/Graph.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/types.hh>

// Forward declarations for friendship:
#include <core/pack/guidance_scoreterms/buried_unsat_penalty/BuriedUnsatPenaltyTests.fwd.hh>
#include <core/pack/guidance_scoreterms/buried_unsat_penalty/graph/BuriedUnsatPenaltyGraphTests.fwd.hh>

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace buried_unsat_penalty {
namespace graph {

/// @brief A hydrogen bond between two atoms in two residues.  Since this is tied to a BuriedUnsatPenaltyEdge (which
/// in turn is associated with two residues), it stores only donor group index and acceptor group index.  The group
/// indices match the indexing in the corresponding node.
/// @details A "donor group" can be a single atom or a collection of related atoms (e.g. the pair of NH2 protons in
/// an aspartate side-chain).  This allows us to impose a little bit of prior knowledge: sometimes, you're satisfied if
/// one of the protons in a group is making a hydrogen bond.
class BuriedUnsatPenaltyGraphHbond : public utility::pointer::ReferenceCount {

	friend class ::BuriedUnsatPenaltyGraphTests; //To allow unit tests to interrogate inner workings of this class.
	friend class ::BuriedUnsatPenaltyGraphSymmetricTests; //To allow unit tests to interrogate inner workings of this class.
	friend class ::BuriedUnsatPenaltyTests; //To allow unit tests to interrogate inner workings of this class.
	friend class ::BuriedUnsatPenaltySymmetricTests; //To allow unit tests to interrogate inner workings of this class.

public:

	/// @brief Default constructor.
	BuriedUnsatPenaltyGraphHbond();

	/// @brief Initialization constructor.
	BuriedUnsatPenaltyGraphHbond( bool const first_node_is_the_acceptor, core::Size const donor_group_index, core::Size const acceptor_group_index, core::Real const energy, core::Size const lower_numbered_node_symmetry_copy_index, core::Size const higher_numbered_node_symmetry_copy_index );

	/// @brief Copy constructor.
	BuriedUnsatPenaltyGraphHbond( BuriedUnsatPenaltyGraphHbond const &src );

	/// @brief Destructor.
	virtual ~BuriedUnsatPenaltyGraphHbond();

	/// @brief Assignment operator.
	BuriedUnsatPenaltyGraphHbond & operator = ( BuriedUnsatPenaltyGraphHbond const & src );

public: //Accessors

	/// @brief Get the index of the Hbond donor group in the donor residue.
	inline core::Size donor_group() const { return donor_group_index_; }

	/// @brief Get the index of the Hbond acceptor group in the donor residue.
	inline core::Size acceptor_group() const { return acceptor_group_index_; }

	/// @brief Is the first node the acceptor?
	inline bool first_node_is_the_acceptor() const { return first_node_is_the_acceptor_; }

	/// @brief Get the symmetry copy index of the lower-numbered node.
	/// @details Returns 1 in the asymmetric case.
	inline core::Size lower_numbered_node_symmetry_copy_index() const { return lower_numbered_node_symmetry_copy_index_; }

	/// @brief Get the symmetry copy index of the higher-numbered node.
	/// @details Returns 1 in the asymmetric case.
	inline core::Size higher_numbered_node_symmetry_copy_index() const { return higher_numbered_node_symmetry_copy_index_; }

	/// @brief Get the symmetry copy index of the node that's the acceptor.
	/// @details Returns 1 in the asymmetric case.
	core::Size acceptor_symmetry_copy_index() const;

	/// @brief Get the symmetry copy index of the node that's the donor.
	/// @details Returns 1 in the asymmetric case.
	core::Size donor_symmetry_copy_index() const;

private:

	/// @brief Is the first node the acceptor?
	bool first_node_is_the_acceptor_;

	/// @brief The index of the Hbond donor group in the donor residue.
	core::Size donor_group_index_;

	/// @brief The index of the Hbond acceptor group in the acceptor residue.
	core::Size acceptor_group_index_;

	/// @brief The energy of the hydrogen bond.
	core::Real energy_;

	/// @brief The symmetry copy index of the lower-numbered node.
	/// @details Always 1 in the asymmetric case.
	core::Size lower_numbered_node_symmetry_copy_index_;

	/// @brief The symmetry copy index of the higher-numbered node.
	/// @details Always 1 in the asymmetric case.
	core::Size higher_numbered_node_symmetry_copy_index_;

}; //class BuriedUnsatPenaltyGraphHbond


/// @brief A class for a hydrogen bond donor group or acceptor group.
/// @details An "acceptor group" is (usually) just an acceptor atom (with the exception noted below).  A "donor group" is a donor proton and its
/// parent heavyatom, plus any other protons on that parent.
/// @note In rare cases (e.g. a hydroxyl), a group can be both donor and acceptor.
class BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup : public utility::pointer::ReferenceCount {

	friend class ::BuriedUnsatPenaltyGraphTests; //To allow unit tests to interrogate inner workings of this class.
	friend class ::BuriedUnsatPenaltyGraphSymmetricTests; //To allow unit tests to interrogate inner workings of this class.
	friend class ::BuriedUnsatPenaltyTests; //To allow unit tests to interrogate inner workings of this class.
	friend class ::BuriedUnsatPenaltySymmetricTests; //To allow unit tests to interrogate inner workings of this class.

public:

	/// @brief Default constructor -- explicitly deleted.
	BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup() = delete;

	/// @brief Options constructor.
	/// @details Note that protons are initialized to an empty list, and must be added later.
	BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup( bool const is_acceptor, bool const is_counted, core::Size const heavyatom_index, core::conformation::Residue const &residue );

	/// @brief Copy constructor -- explicitly deleted.
	BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup( BuriedUnsatPenaltyGraphHbond const &src ) = delete;

	/// @brief Destructor.
	virtual ~BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup();

	/// @brief Assignment operator.
	BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup & operator=( BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup const &src );

public: //Accessor functions

	/// @brief Is this group a donor?
	inline bool is_donor() const { return !proton_indices().empty(); }

	/// @brief Is this group an acceptor?
	inline bool is_acceptor() const { return is_acceptor_; }

	/// @brief Is this group a group to count (e.g. a buried group)?
	inline bool is_counted() const { return is_counted_; }

	/// @brief Get the atom index in the relevant residue of the heavyatom for this group.
	inline core::Size heavyatom_index() const { return heavyatom_index_; }

	/// @brief Get the vector of atom indices in the relevant residue of the protons attached to the heavyatom.
	inline utility::vector1< core::Size > const & proton_indices() const { return proton_indices_; }

	/// @brief Get the number of protons in this group.
	/// @detais Synonymous with max_donated_hbond_count().
	inline core::Size n_protons() const { return proton_indices_.size(); }

	/// @brief Add a proton to the list of protons in this group.
	void add_proton_index( core::Size const index_in );

	/// @brief Determine the maximum number of hydrogen bonds that a group can accept, based on the identity of the heavyatom.
	/// @details Oxygens can accept 2; nitrogens can accept 1.  For now, the rule is that simple.  At some point, we might look up the information from a database
	/// lookup table based on the details of the oxygen type or nitrogen type or whatnot.
	/// @note Returns 0 if "is_acceptor" is false.  This is a static function.
	static core::Size determine_max_accepted_hbond_count( bool const is_acceptor, core::Size const heavyatom_index, core::conformation::Residue const & residue );

	/// @brief Get the maximum donated hydrogen bond count.
	/// @details Synonym for n_protons().
	inline core::Size max_donated_hbond_count() const { return proton_indices_.size(); }

	/// @brief Get the maximum accepted hydrogen bond count.
	inline core::Size max_accepted_hbond_count() const { return max_accepted_hbond_count_; }

private:

	/// @brief Is this group an acceptor?
	bool is_acceptor_;

	/// @brief Is this a group to count (e.g. a buried group)?
	bool is_counted_;

	/// @brief The atom index in the relevant residue of the heavyatom for this group.
	core::Size heavyatom_index_;

	/// @brief The atom indices in the relevant residue of the protons attached to the heavyatom.
	utility::vector1< core::Size > proton_indices_;

	/// @brief The maximum number of hbonds that this group can accept.
	core::Size max_accepted_hbond_count_;

}; //class BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup


/// @brief Data stored inside a BuriedUnsatPenaltyNode.
/// @details This is a separate class so that nodes can rapidly be "copied", just by changing the stored pointer
/// to the data object.
class BuriedUnsatPenaltyNodeData : public utility::pointer::ReferenceCount {

	friend class ::BuriedUnsatPenaltyGraphTests; //To allow unit tests to interrogate inner workings of this class.
	friend class ::BuriedUnsatPenaltyGraphSymmetricTests; //To allow unit tests to interrogate inner workings of this class.
	friend class ::BuriedUnsatPenaltyTests; //To allow unit tests to interrogate inner workings of this class.
	friend class ::BuriedUnsatPenaltySymmetricTests; //To allow unit tests to interrogate inner workings of this class.

public:

	/// @brief Default constructor.
	BuriedUnsatPenaltyNodeData();

	/// @brief Copy constructor (default).
	BuriedUnsatPenaltyNodeData( BuriedUnsatPenaltyNodeData const & ) = default;

	/// @brief Destructor.
	virtual ~BuriedUnsatPenaltyNodeData();

	/// @brief Create a copy by owning pointer and return the owning pointer to the copy.
	virtual BuriedUnsatPenaltyNodeDataOP clone() const;

public: //Setters and getters

	/// @brief Set the residue position.
	void residue_position( core::Size const &setting );

	/// @brief Get the residue position.
	inline core::Size residue_position() const { return residue_position_; }

	/// @brief Set the residue position.
	void rotamer_index( core::Size const &setting );

	/// @brief Get the residue position.
	inline core::Size rotamer_index() const { return rotamer_index_; }

	/// @brief Set the options pointer.
	/// @details This copies the input pointer; it doesn't clone the object.
	void options( BuriedUnsatPenaltyGraphOptionsCOP options_in );

	/// @brief Get the options.
	inline BuriedUnsatPenaltyGraphOptionsCOP options() const { return options_; }

	/// @brief Set the hbonds options pointer.
	/// @details This copies the input pointer; it doesn't clone the object.
	void hbond_options( core::scoring::hbonds::HBondOptionsCOP options_in );

	/// @brief Get the options.
	inline core::scoring::hbonds::HBondOptionsCOP hbond_options() const { return hbond_options_; }

	/// @brief Non-const access to the donor/acceptor groups.
	inline utility::vector1< BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup > & donor_acceptor_groups() { return donor_acceptor_groups_; }

	/// @brief Const access to the donor/acceptor groups.
	inline utility::vector1< BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup > const & donor_acceptor_groups() const { return donor_acceptor_groups_; }

	/// @brief Access the group_heavyatom_index_to_group_index_ map (non-const).
	inline std::map< core::Size, core::Size > & group_heavyatom_index_to_group_index() { return group_heavyatom_index_to_group_index_; }

	/// @brief Access the group_heavyatom_index_to_group_index_ map (const).
	inline std::map< core::Size, core::Size > const & group_heavyatom_index_to_group_index() const { return group_heavyatom_index_to_group_index_; }

	/// @brief Access the donor_acceptor_groups_intrares_hbonds_donated_ vector (non-const).
	inline utility::vector1< core::Size > & donor_acceptor_groups_intrares_hbonds_donated() { return donor_acceptor_groups_intrares_hbonds_donated_; }

	/// @brief Access the donor_acceptor_groups_intrares_hbonds_donated_ vector (const).
	inline utility::vector1< core::Size > const & donor_acceptor_groups_intrares_hbonds_donated() const { return donor_acceptor_groups_intrares_hbonds_donated_; }

	/// @brief Access the donor_acceptor_groups_intrares_hbonds_accepted_ vector (non-const).
	inline utility::vector1< core::Size > & donor_acceptor_groups_intrares_hbonds_accepted() { return donor_acceptor_groups_intrares_hbonds_accepted_; }

	/// @brief Access the donor_acceptor_groups_intrares_hbonds_accepted_ vector (const).
	inline utility::vector1< core::Size > const & donor_acceptor_groups_intrares_hbonds_accepted() const { return donor_acceptor_groups_intrares_hbonds_accepted_; }

	/// @brief Access the donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_donated_ vector (non-const).
	inline utility::vector1< core::Size > & donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_donated() { return donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_donated_; }

	/// @brief Access the donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_donated_ vector (const).
	inline utility::vector1< core::Size > const & donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_donated() const { return donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_donated_; }

	/// @brief Access the donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_accepted_ vector (non-const).
	inline utility::vector1< core::Size > & donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_accepted() { return donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_accepted_; }

	/// @brief Access the donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_accepted_ vector (const).
	inline utility::vector1< core::Size > const & donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_accepted() const { return donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_accepted_; }

private:

	/// @brief The residue index (in the pose) of the residue represented by this node.
	core::Size residue_position_;

	/// @brief The rotmaer index represented by this node.  Will always be 1 if this is a non-packable position.
	core::Size rotamer_index_;

	/// @brief The donor and acceptor groups in the residue.
	utility::vector1< BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup > donor_acceptor_groups_;

	/// @brief Given the group heavyatom index (atom numbering for the residue), get the group index (in the node).
	std::map< core::Size, core::Size > group_heavyatom_index_to_group_index_;

	/// @brief The base count for the number of intraresidue hydrogen bonds for which each donor/acceptor group is a DONOR.
	utility::vector1< core::Size > donor_acceptor_groups_intrares_hbonds_donated_;

	/// @brief The base count for the number of intraresidue hydrogen bonds for which each donor/acceptor group is an ACCEPTOR.
	utility::vector1< core::Size > donor_acceptor_groups_intrares_hbonds_accepted_;

	/// @brief The base count for the number of interresidue hydrogen bonds to symmetric copies of this residue for which each donor/acceptor group is a DONOR.
	utility::vector1< core::Size > donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_donated_;

	/// @brief The base count for the number of interresidue hydrogen bonds to symmetric copies of this residue for which each donor/acceptor group is an ACCEPTOR.
	utility::vector1< core::Size > donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_accepted_;

	/// @brief A const owning pointer to an options container, containing settings for the BuriedUnsatPenaltyGraph.
	BuriedUnsatPenaltyGraphOptionsCOP options_;

	/// @brief A const owning pointer to an options container for hydrogen bond calculations.
	core::scoring::hbonds::HBondOptionsCOP hbond_options_;

}; //class BuriedUnsatPenaltyNodeData

/// @brief Each BuriedUnsatPenaltyNode represents a rotamer from the RotamerSets object, or a residue (if we're doing a simple scoring pass).
/// @details Internally, this stores a list of hydrogen bond donor and acceptor atoms that must be satisfied.
class BuriedUnsatPenaltyNode : public utility::graph::Node {

	friend class ::BuriedUnsatPenaltyGraphTests; //To allow unit tests to interrogate inner workings of this class.
	friend class ::BuriedUnsatPenaltyGraphSymmetricTests; //To allow unit tests to interrogate inner workings of this class.
	friend class ::BuriedUnsatPenaltyTests; //To allow unit tests to interrogate inner workings of this class.
	friend class ::BuriedUnsatPenaltySymmetricTests; //To allow unit tests to interrogate inner workings of this class.

public:

	/// @brief Default constructor -- explicitly deleted.
	BuriedUnsatPenaltyNode() = delete;

	/// @brief Node constructor.
	BuriedUnsatPenaltyNode( utility::graph::Graph * owner, platform::Size const node_id );

	/// @brief Copy constructor -- explicitly deleted.
	BuriedUnsatPenaltyNode( BuriedUnsatPenaltyNode const &src ) = delete;

	/// @brief Destructor.
	virtual ~BuriedUnsatPenaltyNode();

public: //These functions are public, but are only intended to be called by BuriedUnsatPenaltyGraph class:

	/// @brief Initialize this node, setting its residue_position, rotamer_index,
	/// @details The options owning pointer is stored directly; the options are not cloned.
	void initialize_node( core::Size const residue_position, core::Size const rotamer_index, core::conformation::ResidueCOP residue, core::pose::Pose const &pose, BuriedUnsatPenaltyGraphOptionsCOP options, core::scoring::hbonds::HBondOptionsCOP hbond_options, core::scoring::hbonds::HBondDatabase const &hbond_database, bool const is_symmetric );

	/// @brief Given another node, set this node to copy that one.
	/// @details Note: this does NOT copy the rotamer index, but instead sets it to 1.  The donated_hbond_count_
	/// and accepted_hbond_count_ vars are initialized to the intra-residue donated and intra-residue accepted
	/// counts, respectively (i.e. inter-residue donated and inter-residue accepted counts have *not* yet been
	/// added in).
	void copy_from( utility::graph::Node const * other_node ) override;

	/// @brief Given the index of a donor or acceptor heavyatom, get the donor/acceptor group index.
	inline core::Size get_donor_acceptor_group_from_heavyatom_index( core::Size const heavyatom_index ) const { return stored_data_->group_heavyatom_index_to_group_index().at(heavyatom_index); }

	/// @brief Get a donor or acceptor group, by donor or acceptor group index.
	/// @details This version provides const access.
	inline BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup const & donor_acceptor_group( core::Size const donor_acceptor_group_index ) const {
		debug_assert( donor_acceptor_group_index > 0 && donor_acceptor_group_index <= stored_data_->donor_acceptor_groups().size() );
		return stored_data_->donor_acceptor_groups()[donor_acceptor_group_index];
	}

	/// @brief For every donor/acceptor group, reset the counts of hydrogen bonds to/from that group.
	/// @details The donated_hbond_count_ and accepted_hbond_count_ vars are initialized to the
	/// intra-residue donated and intra-residue accepted counts, respectively (i.e. inter-residue
	/// donated and inter-residue accepted counts have *not* yet been added in).
	void clear_hbond_counts();

	/// @brief Add one to the number of accepted hydrogen bonds for the Nth donor/acceptor group in this node.
	void increment_accepted_hbond_count_for_group( core::Size const group_index );

	/// @brief Add one to the number of donated hydrogen bonds for the Nth donor/acceptor group in this node.
	void increment_donated_hbond_count_for_group( core::Size const group_index );

	/// @brief After hydrogen bonds from edges have been counted, report the number of donor/acceptor groups in each of several categories.
	/// @details The integer instances passed in are incremented appropriately by this function.  (So if there is an incoming count of 5 unsaturated acceptors, and
	/// the current node has 3 more, the unsat_acceptor_count will end up being 8).
	/// @param[out] unsat_acceptor_count The number of acceptor (and not donor) groups that are unsatisfied.
	/// @param[out] unsat_donor_count The number of donor (and not acceptor) groups that are unsatisfied.
	/// @param[out] unsat_acceptor_and_donor_count The number of groups that are both donors and acceptors (e.g. hydroxyls) that are unsatisfied (i.e. lack either a donated hbond or an accepted hbond).
	/// @param[out] oversat_acceptor_count The number of acceptor (and not donor) groups that are oversatisfied.
	/// @param[out] oversat_donor_count The number of donor (and not acceptor) groups that are oversatisfied.  (This generally doesn't happen).
	/// @param[out] oversat_acceptor_and_donor_count The number of groups that are both donors and acceptors (e.g. hydroxyls) that are unsatisfied (i.e. have either too many donated hbonds or too many accepted hbonds).
	void increment_counts(core::Size &unsat_acceptor_count, core::Size &unsat_donor_count, core::Size &unsat_acceptor_and_donor_count, core::Size &oversat_acceptor_count, core::Size &oversat_donor_count, core::Size &oversat_acceptor_and_donor_count ) const;

	/// @brief After hydrogen bonds from edges have been counted, report the number of donor/acceptor groups in each of several categories.
	/// @details The integer instances passed in are DECREMENTED appropriately by this function.  (So if there is an incoming count of 5 unsaturated acceptors, and
	/// the current node accounts for 3 unsaturated acceptors, the unsat_acceptor_count will drop to 2).
	/// @param[out] unsat_acceptor_count The number of acceptor (and not donor) groups that are unsatisfied.
	/// @param[out] unsat_donor_count The number of donor (and not acceptor) groups that are unsatisfied.
	/// @param[out] unsat_acceptor_and_donor_count The number of groups that are both donors and acceptors (e.g. hydroxyls) that are unsatisfied (i.e. lack either a donated hbond or an accepted hbond).
	/// @param[out] oversat_acceptor_count The number of acceptor (and not donor) groups that are oversatisfied.
	/// @param[out] oversat_donor_count The number of donor (and not acceptor) groups that are oversatisfied.  (This generally doesn't happen).
	/// @param[out] oversat_acceptor_and_donor_count The number of groups that are both donors and acceptors (e.g. hydroxyls) that are unsatisfied (i.e. have either too many donated hbonds or too many accepted hbonds).
	/// @param[in] data The BuriedUnsatPenaltyNodeData to use.  Since this has often been replaced by the time this function is called, we need to pass it in separately.
	void decrement_counts(core::Size &unsat_acceptor_count, core::Size &unsat_donor_count, core::Size &unsat_acceptor_and_donor_count, core::Size &oversat_acceptor_count, core::Size &oversat_donor_count, core::Size &oversat_acceptor_and_donor_count, graph::BuriedUnsatPenaltyNodeDataCOP data ) const;

	/// @brief Get the residue index (position in the pose) for this node.
	inline core::Size residue_position() const { return stored_data_->residue_position(); }

	/// @brief Get the rotamer index for this node.
	inline core::Size rotamer_index() const { return stored_data_->rotamer_index(); }

	/// @brief Get the number of donor/acceptor groups.
	inline core::Size num_donor_acceptor_groups() const { return stored_data_->donor_acceptor_groups().size(); }

	/// @brief Const access to the stored data.
	inline BuriedUnsatPenaltyNodeDataCOP stored_data() const { return stored_data_; }

private: //Private functions

	/// @brief Determine whether an atom is buried by the method of sidechain neighbors.
	bool is_buried( core::Size const atom_index, core::conformation::Residue const & residue, core::pose::Pose const &pose ) const;

	/// @brief Given a residue, detect intramolecular hydrogen bonds.
	/// @details Populates the donor_acceptor_groups_intrares_hbonds_donated_ and donor_acceptor_groups_intrares_hbonds_accepted_ lists.
	/// @note This version is for the ASYMMETRIC and SYMMETRIC cases.  It only detects the intraresidue hbonds.  In the symmetric case,
	/// the detect_intra_residue_hydrogen_bonds_symmetric() function should be called immediately after this to detect the interresidue hbonds
	/// between this residue and its symmetry mates.
	void detect_intra_residue_hydrogen_bonds( core::conformation::Residue const &residue, core::scoring::hbonds::HBondDatabase const &hbond_database, BuriedUnsatPenaltyNodeDataOP new_stored_data );

	/// @brief Given a residue, detect intramolecular hydrogen bonds.
	/// @details Populates the donor_acceptor_groups_intrares_hbonds_donated_ and donor_acceptor_groups_intrares_hbonds_accepted_ lists.
	/// @note This version is for the SYMMETRIC case.  It detects interresidue hbonds between this residue and its symmetry mates.  This should be
	/// called after detect_intra_residue_hydrogen_bonds().
	void detect_intra_residue_hydrogen_bonds_symmetric( core::Size const res_index, core::conformation::ResidueCOP rotamer, core::pose::Pose const &pose, core::conformation::symmetry::SymmetryInfoCOP const symminfo, core::scoring::hbonds::HBondDatabase const &hbond_database, core::scoring::hbonds::HBondOptions const &hbond_options, BuriedUnsatPenaltyNodeDataOP new_stored_data );


private:

	/// @brief Owning pointer to the data object.
	BuriedUnsatPenaltyNodeDataCOP stored_data_;

	/// @brief The number of hbonds that each group is donating.
	utility::vector1< core::Size > donated_hbond_count_;

	/// @brief The number of hbonds that each group is accepting.
	utility::vector1< core::Size > accepted_hbond_count_;

}; //BuriedUnsatPenaltyNode


class BuriedUnsatPenaltyEdgeData : public utility::pointer::ReferenceCount {

	friend class ::BuriedUnsatPenaltyGraphTests; //To allow unit tests to interrogate inner workings of this class.
	friend class ::BuriedUnsatPenaltyGraphSymmetricTests; //To allow unit tests to interrogate inner workings of this class.
	friend class ::BuriedUnsatPenaltyTests; //To allow unit tests to interrogate inner workings of this class.
	friend class ::BuriedUnsatPenaltySymmetricTests; //To allow unit tests to interrogate inner workings of this class.

public:

	/// @brief Constructor.
	BuriedUnsatPenaltyEdgeData();

	/// @brief Copy constructor -- explicitly deleted.
	BuriedUnsatPenaltyEdgeData( BuriedUnsatPenaltyEdgeData const & ) = delete;

	/// @brief Destructor.
	virtual ~BuriedUnsatPenaltyEdgeData() = default;

public:

	/// @brief Access the hbonds list (const access).
	inline utility::vector1< BuriedUnsatPenaltyGraphHbond > const &hbonds_list() const { return hbonds_list_; }

	/// @brief Add a hydrogen bond to a newly-created edge data object.
	/// @details Note that acceptor_group and donor_group are group indices in the respective nodes, not atom indices in the respective residues.
	/// @note The symmetry copy indices should both be 1 in the asymmetric case.
	void add_hbond( bool const lower_numbered_node_is_acceptor, core::Size const acceptor_group, core::Size const donor_group, core::Real const hbond_energy, core::Size const lower_numbered_node_symmetry_copy_index, core::Size const higher_numbered_node_symmetry_copy_index);

private:

	/// @brief A list of all of the hydrogen bonds between two residues.
	utility::vector1< BuriedUnsatPenaltyGraphHbond > hbonds_list_;

};

/// @brief Each BuriedUnsatPenaltyEdge represents a hydrogen bonding interaction between two residues, and stores
/// information about (a) the number of hydrogen bonds, and (b) the atoms involved.
class BuriedUnsatPenaltyEdge : public utility::graph::Edge {

	friend class ::BuriedUnsatPenaltyGraphTests; //To allow unit tests to interrogate inner workings of this class.
	friend class ::BuriedUnsatPenaltyGraphSymmetricTests; //To allow unit tests to interrogate inner workings of this class.
	friend class ::BuriedUnsatPenaltyTests; //To allow unit tests to interrogate inner workings of this class.
	friend class ::BuriedUnsatPenaltySymmetricTests; //To allow unit tests to interrogate inner workings of this class.

public:

	/// @brief Constructor -- explicitly deleted.
	BuriedUnsatPenaltyEdge() = delete;

	/// @brief Edge constructor.
	BuriedUnsatPenaltyEdge( utility::graph::Graph * owner, platform::Size const first_node_ind, platform::Size const second_node_ind );

	/// @brief Edge copy-like constructor.
	BuriedUnsatPenaltyEdge( utility::graph::Graph * owner, BuriedUnsatPenaltyEdge const &src );

	/// @brief Copy constructor -- explicitly deleted.
	BuriedUnsatPenaltyEdge( BuriedUnsatPenaltyEdge const &src ) = delete;

	/// @brief Destructor.
	virtual ~BuriedUnsatPenaltyEdge();

public:

	/// @brief Initialize this edge from another.
	void copy_from( utility::graph::Edge const * src) override;

	/// @brief Get the number of hbonds in this edge.
	inline core::Size n_hbonds() const { return edge_data_->hbonds_list().size(); }

	/// @brief Access a particular hbond in this edge.
	inline BuriedUnsatPenaltyGraphHbond const & hbond( core::Size const hbond_index ) const {
		debug_assert( hbond_index > 0 && hbond_index <= edge_data_->hbonds_list().size() );
		return edge_data_->hbonds_list()[hbond_index];
	}

	/// @brief Set the data object for this edge.
	/// @details Data object pointer is copied; object is NOT cloned.
	void set_edge_data( BuriedUnsatPenaltyEdgeDataCOP edge_data_in );

private:

	/// @brief A container for the hydrogen bond data, which can be shared with other edges.
	BuriedUnsatPenaltyEdgeDataCOP edge_data_;

}; //BuriedUnsatPenaltyEdge class


/// @brief The BuriedUnsatPenaltyGraph consists of nodes representing resiudues (or rotamers in packing mode) and edges representing hydrogen-bonding
/// interactions.  Each node stores a list of hydrogen bond donor and acceptor atoms that we're seeking to satisfy; each edge stores the donors and acceptors
/// that are connected by hydrogen bonds when two residues (or rotamers) interact.
class BuriedUnsatPenaltyGraph : public utility::graph::Graph {

	friend class ::BuriedUnsatPenaltyGraphTests; //To allow unit tests to interrogate inner workings of this class.
	friend class ::BuriedUnsatPenaltyGraphSymmetricTests; //To allow unit tests to interrogate inner workings of this class.
	friend class ::BuriedUnsatPenaltyTests; //To allow unit tests to interrogate inner workings of this class.
	friend class ::BuriedUnsatPenaltySymmetricTests; //To allow unit tests to interrogate inner workings of this class.

public:

	typedef utility::graph::Graph parent;

public:

	/// @brief Default constructor -- explicitly deleted.
	BuriedUnsatPenaltyGraph() = delete;

	/// @brief Options constructor.
	/// @details Note: this stores the owning pointer to the options; it doesn't clone them.  The hbond options *are* cloned and modified, though.
	BuriedUnsatPenaltyGraph( BuriedUnsatPenaltyGraphOptionsCOP options, core::scoring::hbonds::HBondOptionsCOP hbond_options );

	/// @brief Nodecount constructor with options.
	/// @details Note: this stores the owning pointer to the options; it doesn't clone them.  The hbond options *are* cloned and modified, though.
	BuriedUnsatPenaltyGraph( platform::Size const num_nodes, BuriedUnsatPenaltyGraphOptionsCOP options, core::scoring::hbonds::HBondOptionsCOP hbond_options );

	/// @brief Copy constructor.
	BuriedUnsatPenaltyGraph( BuriedUnsatPenaltyGraph const &src );

	/// @brief Destructor.
	~BuriedUnsatPenaltyGraph() override;

public: //Static helper functions

	/// @brief Set up options for hbond detection.
	static void configure_hbond_options( core::scoring::hbonds::HBondOptions &hbondoptions );

public: //Public member functions

	/// @brief Needed override from base class.
	void delete_edge( utility::graph::Edge * edge ) override;

	/// @brief Needed override from base class.
	void delete_node( utility::graph::Node * node ) override;

	/// @brief Set whether this is a graph that just stores one rotamer per position.
	void set_always_rotamer_one( bool const setting );

	/// @brief Provide Pymol commands to colour the pose grey, non-buried donor and acceptor groups cyan, and buried acceptor
	/// and donor groups orange.  Useful for debugging degree of burial.
	/// @details To use, pass in a pose.  If this graph contains residues corresponding to those in the pose, commands for colouring
	/// them will be written out.
	void provide_pymol_commands_to_show_groups( std::ostream &out, core::pose::Pose const &pose ) const;

	/// @brief Initialize a BuriedUnsatPenaltyGraph from a pose, for scoring.
	void initialize_graph_for_scoring( core::pose::Pose const &pose );

	/// @brief Initialize a BuriedUnsatPenaltyGraph from a pose and a residue set, for packing.
	void initialize_graph_for_packing( core::pose::Pose const &pose, core::pack::rotamer_set::RotamerSets const &rotamersets, bool const only_scoring=false );

	/// @brief Given this BuriedUnsatPenaltyGraph with some number of nodes, iterate through each node and update the
	/// internally-stored counts for unsats and oversats based on the edges connected to that node.
	/// @details Calls compute_unsats_for_node().
	void compute_unsats_all_nodes();

	/// @brief Given two lists (one of changed nodes, one of their partners), update the internally-stored counts for unsats and oversats for
	/// those nodes only.
	/// @details calls compute_unsats_for_node().
	void compute_unsats_changed_nodes( utility::vector1< core::Size > const & changed_node_indices, utility::vector1< core::Size > const & changed_node_partners );

	/// @brief Given this BuriedUnsatPenaltyGraph with some number of nodes and the index of a node, update the
	/// internally-stored counts for unsats and oversats based on the edges connected to that node.
	void compute_unsats_for_node( core::Size const node_index );

	/// @brief Given an index of a node in this graph, an owning pointer to another graph, and a node index in the other graph, copy the node
	/// from the other graph to the node in this graph, flush the edges that were connected to the node in this graph, and copy those edges from
	/// the other graph that can be connected to nodes in this graph.
	/// @details Note that the logic for determining whether an edge from the other graph can be related to this graph is based on ResidueCOPs.
	void copy_node_and_connected_edges( core::Size const node_index_in_this_graph, BuriedUnsatPenaltyGraph const & other_graph, core::Size const node_index_in_other_graph );

	/// @brief Given the sequence position and rotamer index of a residue, get the corresponding node index.
	inline core::Size get_node_index( core::Size const seqpos, core::Size const rotamer_index ) const {
		std::pair< core::Size, core::Size > address( seqpos, always_rotamer_one_ ? 1 : rotamer_index );
		runtime_assert( residuepos_rotamerindex_to_nodeindex_.find(address) != residuepos_rotamerindex_to_nodeindex_.end() );
		return residuepos_rotamerindex_to_nodeindex_.at(address);
	}

	/// @brief Given a const owning pointer to a residue, determine whether a corresponding node exists.
	inline bool has_node_corresponding_to_residue( core::conformation::ResidueCOP residue ) const {
		return (residue_memory_address_to_nodeindex_.find( residue ) != residue_memory_address_to_nodeindex_.end() );
	}

	/// @brief Given a const owning pointer to a residue, get the corresponding node index.
	/// @details Will throw an error if the owning pointer corresponds to no node in this object.
	inline core::Size get_node_index( core::conformation::ResidueCOP residue ) const {
		return residue_memory_address_to_nodeindex_.at( residue );
	}

	/// @brief Get the memory address of the residue corresponding to a particular node.
	inline core::conformation::ResidueCOP nodeindex_to_residue_memory_address( core::Size const nodeindex ) const {
		return nodeindex_to_residue_memory_address_.at(nodeindex);
	}

private: //Private functions

	/// @brief Add an edge, representing an interresidue hydrogen bonded interaction (consisting of one or more interresidue hydrogen bonds), to the graph.
	/// @details Returns a pointer to the newly-created edge.
	/// @note In the graph base class, the first node index is always numbered lower than the second node index for an edge.  The swap happens automatically.
	BuriedUnsatPenaltyEdge * add_edge( core::Size const seqpos1, core::Size const rotamer_index1, core::Size const seqpos2, core::Size const rotamer_index2 );

	/// @brief Add an edge, representing an interresidue hydrogen bonded interaction (consisting of one or more interresidue hydrogen bonds), to the graph.
	/// @details Returns a pointer to the newly-created edge.  This version copies an edge from another graph.
	/// @note In the graph base class, the first node index is always numbered lower than the second node index for an edge.  The swap happens automatically.
	BuriedUnsatPenaltyEdge * add_edge( core::Size const node_index1, core::Size const node_index2, BuriedUnsatPenaltyEdge const & other_edge );

	/// @brief Determine whether an edge exists in the graph, by seqpos and rotamer index.
	bool get_edge_exists( core::Size const seqpos1, core::Size const rotamer_index1, core::Size const seqpos2, core::Size const rotamer_index2 ) const;

	/// @brief Retrieve an edge by seqpos and rotamer index.
	/// @details Nonconst version.
	BuriedUnsatPenaltyEdge * find_edge( core::Size const seqpos1, core::Size const rotamer_index1, core::Size const seqpos2, core::Size const rotamer_index2 );

	/// @brief Retrieve an edge by seqpos and rotamer index.
	/// @details Const version.
	BuriedUnsatPenaltyEdge const * find_edge( core::Size const seqpos1, core::Size const rotamer_index1, core::Size const seqpos2, core::Size const rotamer_index2 ) const;

protected:

	/// @brief Factory method for node creation, defined by derived graph
	/// classes, called by the base class
	utility::graph::Node* create_new_node( platform::Size node_index ) override;

	/// @brief Factory method for edge creation, defined by derived graph
	/// classes, called by the base class
	utility::graph::Edge* create_new_edge( platform::Size index1, platform::Size index2 ) override;

	/// @brief This is also needed for edge creation, when copying graphs.
	utility::graph::Edge * create_new_edge( utility::graph::Edge const * example_edge ) override;

private:

	/// @brief Initialize a node to represent a given residue position and rotamer index, and store in it all of the relevant hydrogen bond donors and acceptors.
	/// @param[in] node_index The index of the node to initialize.
	/// @param[in] residue_position The index of the residue that this node represents in the pose.
	/// @param[in] rotamer_index The index of the rotamer for this residue that this node represents.
	/// @param[in] residue The residue object itself, for extracting hbond donor/acceptor information.
	/// @param[in] pose The pose, for context.
	/// @param[in] hb_data The hydrogen bonding database object.
	/// @param[in] is_symmetric Is this a symmetric pose?  (We figure this out once and only once, to avoid repeated dynamic_casts.)
	void initialize_node( core::Size const node_index, core::Size const residue_position, core::Size const rotamer_index, core::conformation::ResidueCOP residue, core::pose::Pose const &pose, core::scoring::hbonds::HBondDatabase const &hb_data, bool const is_symmetric );

private: //Data members:

	/// @brief In some cases, the graph should only store one rotamer per position.
	bool always_rotamer_one_;

	/// @brief A map of (residue position, rotamer index)-->(node index).
	std::map< std::pair< core::Size, core::Size >, core::Size > residuepos_rotamerindex_to_nodeindex_;

	/// @brief A map of (Residue memory address)-->(node index).
	std::map< core::conformation::ResidueCOP, core::Size > residue_memory_address_to_nodeindex_;

	/// @brief A map of (node index)-->(Residue memory address).
	std::map< core::Size, core::conformation::ResidueCOP > nodeindex_to_residue_memory_address_;

	/// @brief A const owning pointer to an options container, containing settings for the BuriedUnsatPenaltyGraph.
	BuriedUnsatPenaltyGraphOptionsCOP options_;

	/// @brief A const owning pointer to an options container for hydrogen bonds.
	core::scoring::hbonds::HBondOptionsCOP hbond_options_;

	/// @brief Storage for the edges of this graph.
	/// @details The Graph base class is weird and nasty.  If your derived edge class differs
	/// from the base edge class, the Graph base clas fails to manage the edges correctly on
	/// destruction, and this creates a memory leak.  The workaround is to manage edges yourself
	/// in the derived Graph class, which is silly.
	boost::unordered_object_pool< BuriedUnsatPenaltyEdge > * bunsat_edge_pool_;

}; //BuriedUnsatPenaltyGraph class


} //graph
} //buried_unsat_penalty_graph
} //guidance_scoreterms
} //pack
} //core

#endif
