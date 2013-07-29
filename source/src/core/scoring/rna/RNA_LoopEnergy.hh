// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RNA_LoopEnergy.hh
/// @brief  Radius of gyration score for RNA, to match Rosetta++
/// @author Rhiju Das


#ifndef INCLUDED_core_scoring_rna_RNA_LoopEnergy_hh
#define INCLUDED_core_scoring_rna_RNA_LoopEnergy_hh


// Package headers
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace rna {


class RNA_LoopEnergy : public methods::WholeStructureEnergy  {
public:
	typedef methods::WholeStructureEnergy  parent;

public:

	/// @brief Defines a loop closure energy based on FullModelInfo and pose geometry.
	RNA_LoopEnergy();

	/// clone
	virtual
	methods::EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	virtual
	void
	setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const;

	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & totals
	) const;

/////////////////////////////////
	void
	eval_atom_derivative(
	  id::AtomID const & atom_id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		ScoreFunction const &,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2 ) const;

	void
	indicate_required_context_graphs(
		utility::vector1< bool > & /*context_graphs_required*/
	) const {}

private:

	void
	update_rna_loop_atoms_and_lengths( pose::Pose & pose ) const;

	Real
	get_loop_distance2( Vector const & v_takeoff,  Vector const & v_landing ) const;

	Vector
	get_loop_distance2_deriv( Vector const & v_takeoff,  Vector const & v_landing ) const;

	virtual
	core::Size version() const;

	Real const rna_loop_fixed_cost_;
	Real const persistence_length2_;
	Real const kB_T_;

	// might be better to cache the following inside the pose...
	mutable utility::vector1< Size > rna_loop_lengths_;
	mutable utility::vector1< core::id::AtomID > loop_takeoff_atoms_;
	mutable utility::vector1< core::id::AtomID > loop_landing_atoms_;

};


}
}
}

#endif // INCLUDED_core_scoring_methods_RNA_LoopEnergy_HH
