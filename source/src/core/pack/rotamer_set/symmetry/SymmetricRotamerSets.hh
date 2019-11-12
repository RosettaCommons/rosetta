// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/RotamerSet/RotamerSets.hh
/// @brief  RotamerSets class declaration, for symmetric packing
/// @author Ingemar Andre


#ifndef INCLUDED_core_pack_rotamer_set_symmetry_SymmetricRotamerSets_hh
#define INCLUDED_core_pack_rotamer_set_symmetry_SymmetricRotamerSets_hh

// Unit Headers
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.fwd.hh>

// Package Headers
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/interaction_graph/SymmOnTheFlyInteractionGraph.fwd.hh>

// Project Headers
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>

// Utility headers
#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace rotamer_set {
namespace symmetry {

typedef utility::vector1< RotamerSetOP > RotamerSetVector;

class SymmetricRotamerSets : public RotamerSets
{
public:
	typedef task::PackerTaskCOP PackerTaskCOP;
	typedef conformation::symmetry::SymmetricConformation SymmetricConformation;
	typedef conformation::symmetry::SymmetryInfoCOP SymmetryInfoCOP;

public:
	SymmetricRotamerSets();
	~SymmetricRotamerSets() override;

	void
	compute_energies(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scfxn,
		utility::graph::GraphCOP packer_neighbor_graph,
		interaction_graph::InteractionGraphBaseOP ig,
		core::Size const threads_to_request
	) override;

public:

	/// @brief Append to a vector of calculations that already contains one-body energy calculations additonal work units
	/// for two-body energy calculations, for subsequent multi-threaded evaluation.
	/// @details Each individual calculation is for the interaction of all possible rotamer pairs at two positions.  Again,
	/// not as granular as conceivably possible, but easier to set up.
	/// @note The work_vector vector is extended by this operation.  This version is for the symmetric case.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitue.org).
	void
	append_two_body_energy_computations_to_work_vector(
		core::pose::Pose const & pose,
		core::scoring::ScoreFunction const & scfxn,
		utility::graph::GraphCOP packer_neighbor_graph,
		interaction_graph::PrecomputedPairEnergiesInteractionGraphOP pig,
		utility::vector1< basic::thread_manager::RosettaThreadFunction > & work_vector,
		basic::thread_manager::RosettaThreadAssignmentInfo const & thread_assignment_info
	) const override;

public:

	//fpd function to set some pose data needed SymmetricRotamerSets
	void
	initialize_pose_for_rotsets_creation(
		pose::Pose & pose
	) const override;


private:
	void
	prepare_symm_otf_interaction_graph(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scfxn,
		utility::graph::GraphCOP packer_neighbor_graph,
		interaction_graph::SymmOnTheFlyInteractionGraphOP ig
	);

	void
	compute_proline_correction_energies_for_otf_graph(
		pose::Pose const & pose,
		conformation::symmetry::SymmetryInfoCOP symm_info,
		scoring::ScoreFunction const & scfxn,
		utility::graph::GraphCOP packer_neighbor_graph,
		interaction_graph::SymmOnTheFlyInteractionGraphOP otfig
	);

public:

	/// @brief Generate a new rotamer set oriented onto a reference rotamer set by cloning the reference set and reorienting
	/// using symmetry information.  Return an owning pointer to the new, reoriented rotamer set.
	/// @details Orients all rotamers in a rotamer_set to a different (symmetrical) position
	/// @note If set_up_mirror_symmetry_if_has_mirror_symmetry_ is true, then residues in mirrored subunits have their
	/// ResidueTypes switched to the types of the opposite chirality.  If false (the default), then they keep the same
	/// types, and only have their geometry mirrored.  The default is false, which is computationally less expensive
	/// at the expense of having the incorrect types in mirrored subunits.
	static
	RotamerSetOP
	orient_rotamer_set_to_symmetric_partner(
		pose::Pose const & pose,
		RotamerSetCOP rotset_for_seqpos,
		uint const & symmpos,
		bool const set_up_mirror_types_if_has_mirror_symmetry=false
	);

private:

	bool
	final_visit_to_edge(
		pose::Pose const & pose,
		utility::graph::GraphCOP packer_neighbor_graph,
		uint ii_resid,
		uint jj_resid
	);

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

}
} // namespace rotamer_set
} // namespace pack
} // namespace core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pack_rotamer_set_symmetry_SymmetricRotamerSets )
#endif // SERIALIZATION


#endif // INCLUDED_core_pack_RotamerSet_RotamerSets_HH
