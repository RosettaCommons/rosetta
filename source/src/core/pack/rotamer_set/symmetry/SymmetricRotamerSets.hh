// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
	~SymmetricRotamerSets();

	virtual
	void
	compute_energies(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scfxn,
		graph::GraphCOP packer_neighbor_graph,
		interaction_graph::InteractionGraphBaseOP ig
	);

	virtual
	void
	compute_one_body_energies(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scfxn,
		graph::GraphCOP packer_neighbor_graph,
		interaction_graph::InteractionGraphBaseOP ig
	);

	/// @brief precomputes all rotamer pair energies between neighboring RotamerSets( residues )
	/// and stores those energies in an intereaction graph capable of storing them
	/// public so it can be used by the GreenPacker.
	virtual
	void
	precompute_two_body_energies(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scfxn,
		graph::GraphCOP packer_neighbor_graph,
		interaction_graph::PrecomputedPairEnergiesInteractionGraphOP pig,
		bool const finalize_edges = true
	);

private:
	void
	prepare_symm_otf_interaction_graph(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scfxn,
		graph::GraphCOP packer_neighbor_graph,
		interaction_graph::SymmOnTheFlyInteractionGraphOP ig
	);

	void
	compute_proline_correction_energies_for_otf_graph(
		pose::Pose const & pose,
		conformation::symmetry::SymmetryInfoCOP symm_info,
		scoring::ScoreFunction const & scfxn,
		graph::GraphCOP packer_neighbor_graph,
		interaction_graph::SymmOnTheFlyInteractionGraphOP otfig
	);

	RotamerSetOP
	orient_rotamer_set_to_symmetric_partner(
		pose::Pose const & pose,
		uint const & setpos,
		uint const & symmpos
	);

	bool
	final_visit_to_edge(
		pose::Pose const & pose,
		graph::GraphCOP packer_neighbor_graph,
		uint ii_resid,
		uint jj_resid
	);

};

}
} // namespace rotamer_set
} // namespace pack
} // namespace core


#endif // INCLUDED_core_pack_RotamerSet_RotamerSets_HH
