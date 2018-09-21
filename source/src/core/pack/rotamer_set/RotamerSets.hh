// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/RotamerSet/RotamerSets.hh
/// @brief  RotamerSets class declaration
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_pack_rotamer_set_RotamerSets_hh
#define INCLUDED_core_pack_rotamer_set_RotamerSets_hh

// Unit Headers
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
//#include <core/pack/rotamer_set/rotamer_building_functions.hh> // WHY IS THIS IN THE .HH?

// Package Headers

#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/pack/rotamer_set/FixbbRotamerSets.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>
#include <core/pack/interaction_graph/PrecomputedPairEnergiesInteractionGraph.fwd.hh>
#include <core/pack/interaction_graph/OnTheFlyInteractionGraph.fwd.hh>

#ifdef WIN32 //VC++ needs full class declaration
#include <core/pack/rotamer_set/RotamerSet.hh> // WIN32 INCLUDE
#endif

// Project Headers
#include <core/conformation/Residue.fwd.hh>
#include <utility/graph/Graph.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility Headers
#include <core/types.hh>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace rotamer_set {

class RotamerSets : public FixbbRotamerSets
{
public:
	typedef task::PackerTaskCOP PackerTaskCOP;

public:
	RotamerSets();
	~RotamerSets();

	void
	set_task( task::PackerTaskCOP task);

	void
	build_rotamers(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scfxn,
		utility::graph::GraphCOP packer_neighbor_graph
	);

	//fd add explicit rotamers at a position.  Deletes current rotamers at this position!
	void
	set_explicit_rotamers(
		core::Size moltenresid,
		RotamerSetOP rotamers
	);

	void
	dump_pdb( pose::Pose const & pose, std::string const & filename ) const;

	virtual uint nrotamers() const;
	virtual uint nrotamers_for_moltenres( uint ) const;

	virtual uint nmoltenres() const;
	virtual uint total_residue() const;

	virtual
	uint
	moltenres_2_resid( uint ) const;

	virtual
	uint
	resid_2_moltenres( uint ) const;

	virtual
	uint
	moltenres_for_rotamer( uint ) const;

	virtual
	uint
	res_for_rotamer( uint ) const;

	virtual
	core::conformation::ResidueCOP
	rotamer( uint ) const;

	virtual
	core::conformation::ResidueCOP
	rotamer_for_moltenres( uint moltenres_id, uint rotamerid ) const;

	virtual
	uint
	nrotamer_offset_for_moltenres( uint ) const;

	/// @brief Does this RotamerSets object store a rotamer set for a residue at position resid
	/// in the pose?
	/// @details Rotamer sets for non-packable residues aren't generated, but could conceivably
	/// be queried by protocols or energy functions.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	bool has_rotamer_set_for_residue( uint const resid ) const;

	virtual
	RotamerSetCOP
	rotamer_set_for_residue( uint resid ) const;

	virtual
	RotamerSetOP
	rotamer_set_for_residue( uint resid );

	virtual
	RotamerSetCOP
	rotamer_set_for_moltenresidue( uint moltenresid ) const;

	virtual
	RotamerSetOP
	rotamer_set_for_moltenresidue( uint moltenresid );

	virtual
	RotamerSetVector::const_iterator begin()
	{ return set_of_rotamer_sets_.begin(); }

	virtual
	RotamerSetVector::const_iterator end()
	{ return set_of_rotamer_sets_.end(); }

	/// convert rotid in full rotamer enumeration into rotamer id on its source residue
	virtual
	uint
	rotid_on_moltenresidue( uint rotid ) const;

	/// convert moltenres rotid to id in full rotamer enumeration
	virtual
	uint
	moltenres_rotid_2_rotid( uint moltenres, uint moltenresrotid ) const;

	/// access to packer_task_
	PackerTaskCOP
	task() const;

	void prepare_sets_for_packing( pose::Pose const & pose, scoring::ScoreFunction const &);

	/// @brief Give the pose a chance to stash any data needed by the _rotset_
	///        need nonconst access to pose
	virtual
	void
	initialize_pose_for_rotsets_creation(
		pose::Pose & /*pose*/
	) const {}

	///fd: make this public
	///    the reason?  we now can update rotamersets interactively from code, so this should be exposed
	void update_offset_data();

	/// @brief Used as a debug assert to ensure that RotamerSetsOperations call update_offset_data()
	/// @details Without this check, rotamer sets operations could potentially introduce very-hard-to-debug
	///          glitches if they don't call update_offset_data().
	bool offset_data_up_to_date();

private:
	void copy_residue_conenctions_and_variants(pose::Pose const & pose, conformation::ResidueOP cloneRes, Size seqpos, Size asym_length);


public:

	/// @brief Precompute all rotamer pair and rotamer one-body energies, populating
	/// the given interaction graph.
	virtual
	void
	compute_energies(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scfxn,
		utility::graph::GraphCOP packer_neighbor_graph,
		interaction_graph::InteractionGraphBaseOP ig
	);

	/// @brief Precompute all rotamer one-body energies, populating the given
	/// interaction graph.
	virtual
	void
	compute_one_body_energies(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scfxn,
		utility::graph::GraphCOP packer_neighbor_graph,
		interaction_graph::InteractionGraphBaseOP ig
	);

	/// @brief Precompute all rotamer pair energies between neighboring RotamerSets (residues)
	/// populating the given interaction graph.
	///
	/// Public so it can be used by the GreenPacker.
	virtual
	void
	precompute_two_body_energies(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scfxn,
		utility::graph::GraphCOP packer_neighbor_graph,
		interaction_graph::PrecomputedPairEnergiesInteractionGraphOP pig,
		bool const finalize_edges = true
	);

private:
	/// @brief Marks all protein vertices in the on-the-fly interaction graph as ones
	/// that should distinguish between the backbones and sidechains.  Then, adds edges to
	/// the on-the-fly interaction graph between neighboring RotamerSets,
	/// and figures out, for those edges, which pairs of rotamer groups (e.g. ala1/arg2, ser1/phe2)
	/// are close enough to interact for their energies to need calculation.
	void
	prepare_otf_graph(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scfxn,
		utility::graph::GraphCOP packer_neighbor_graph,
		interaction_graph::OnTheFlyInteractionGraphOP otfig
	);

	/// @brief computes one body energies for the on-the-fly graph, calculating proline-correction
	/// terms for protien-residues that allow prolines and storing them on the otf edges
	void
	compute_proline_correction_energies_for_otf_graph(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scfxn,
		utility::graph::GraphCOP packer_neighbor_graph,
		interaction_graph::OnTheFlyInteractionGraphOP otfig
	);

public:

	static
	core::PackerEnergy
	get_bb_bbE(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scfxn,
		conformation::Residue const & res1,
		conformation::Residue const & res2
	);

	static
	core::PackerEnergy
	get_sc_bbE(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scfxn,
		conformation::Residue const & res1,
		conformation::Residue const & res2
	);

public: // WHY? G++, WHY?!

	virtual
	utility::vector1< uint > const &
	resid_2_moltenres_vector() const {
		return resid_2_moltenres_;
	}

	virtual
	utility::vector1< uint > const &
	moltenres_2_resid_vector() const {
		return moltenres_2_resid_;
	}

	virtual
	void
	show( std::ostream & out ) const;

private:
	uint nmoltenres_ = 0;
	uint total_residue_ = 0;

	uint nrotamers_ = 0;

	RotamerSetVector set_of_rotamer_sets_;
	utility::vector1< uint > resid_2_moltenres_;
	utility::vector1< uint > moltenres_2_resid_;
	utility::vector1< uint > nrotamer_offsets_;

	// originating moltenres for a particular rotamer in the enumeration of all rotamers
	utility::vector1< uint > moltenres_for_rotamer_;
	utility::vector1< uint > nrotamers_for_moltenres_;

	PackerTaskCOP task_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // namespace rotamer_set
} // namespace pack
} // namespace core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pack_rotamer_set_RotamerSets )
#endif // SERIALIZATION


#endif // INCLUDED_core_pack_RotamerSet_RotamerSets_HH
