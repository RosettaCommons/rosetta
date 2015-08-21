// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/method/EnergyMethod.hh
/// @brief  Base class for energy classes
/// @author Phil Bradley


#ifndef INCLUDED_core_scoring_methods_EnergyMethod_hh
#define INCLUDED_core_scoring_methods_EnergyMethod_hh

// Unit headers
#include <core/scoring/methods/EnergyMethod.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/scoring/methods/EnergyMethodCreator.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>

#include <core/pose/Pose.fwd.hh>

#include <core/kinematics/MinimizerMapBase.fwd.hh>

#include <core/id/AtomID.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

// Numeric headers


namespace core {
namespace scoring {
namespace methods {

/// @brief base class for the energy method hierarchy
class EnergyMethod : public utility::pointer::ReferenceCount  {
public:
	typedef utility::pointer::ReferenceCount  parent;

public:
	/// @brief Constructor with EnergyMethodCreator, which lists the score
	/// types that this energy method is responsible for.
	EnergyMethod( EnergyMethodCreatorOP creator );

	/// @brief Copy constructor copies over the score types of the source
	EnergyMethod( EnergyMethod const & src );

	virtual
	~EnergyMethod();

	virtual
	EnergyMethodOP
	clone() const = 0;

	/// @brief if an energy method needs to cache data in the Energies object,
	/// before packing begins, then it does so during this function. The packer
	/// must ensure this function is called. The default behavior is to do nothing.
	virtual
	void
	setup_for_packing( pose::Pose &, utility::vector1< bool > const &, utility::vector1< bool > const & ) const;

	/// @brief If an energy method needs to cache data in a packing::RotamerSet object before
	/// rotamer energies are calculated, it does so during this function. The packer
	/// must ensure this function is called. The default behavior is to do nothing.
	virtual
	void
	prepare_rotamers_for_packing(
		pose::Pose const &,
		conformation::RotamerSetBase &
	) const;


	/// @brief If the pose changes in the middle of a packing (as happens in rotamer trials) and if
	/// an energy method needs to cache data in the pose that corresponds to its current state,
	/// then the method must update that data when this function is called.  The packer must
	/// ensure this function gets called.  The default behavior is to do nothing.
	virtual
	void
	update_residue_for_packing(
		pose::Pose &,
		Size resid
	) const;

	/// @brief if an energy method needs to cache something in the pose (e.g. in pose.energies()),
	/// before scoring begins, it must do so in this method.  All long range energy
	/// functions must initialize their LREnergyContainers before scoring begins.
	/// The default is to do nothing.
	virtual
	void
	setup_for_scoring( pose::Pose &, ScoreFunction const & ) const;

	/// @brief Called at the beginning of atom tree minimization, this method
	/// allows the derived class the opportunity to initialize pertinent data
	/// that will be used during minimization.  During minimzation, the chemical
	/// structure of the pose is constant, so assumptions on the number of atoms
	/// per residue and their identities are safe so long as the pose's Energies
	/// object's "use_nblist()" method returns true.
	virtual
	void
	setup_for_minimizing(
		pose::Pose & ,
		ScoreFunction const & ,
		kinematics::MinimizerMapBase const &
	) const;

	/// @brief Called immediately before atom- and DOF-derivatives are calculated
	/// allowing the derived class a chance to prepare for future calls.
	virtual
	void
	setup_for_derivatives( pose::Pose & pose, ScoreFunction const & sfxn ) const;

	/// @brief called at the end of derivatives evaluation
	virtual
	void
	finalize_after_derivatives( pose::Pose &, ScoreFunction const &  ) const;

	/// @brief Should this EnergyMethod have score and derivative evaluation
	/// evaluated only in the context of the whole Pose, or can it be included
	/// in a decomposed manner for a residue or a set of residue-pairs that are
	/// not part of the Pose that's serving as their context?  The default
	/// method implemented in the base class returns true in order to grandfather
	/// in EnergyMethods that have not had their derivatives changed to take
	/// advantage of the new derivative-evaluation machinery.  Methods that return
	/// "true" will not have their residue-energy(-ext) / residue-pair-energy(-ext)
	/// methods invoked by the ScoreFunction during its traversal of the
	/// MinimizationGraph, and instead will be asked to perform all their work
	/// during finalize_total_energies().  Similarly, they will be expected to
	/// perform all their work during eval_atom_deriv() instead of during the
	/// ScoreFunction's traversal of the MinimizationGraph for derivative evaluation.
	/// IMPORTANT: Methods that return "true" cannot be included in RTMin.
	virtual
	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const;

	/// @brief Should this EnergyMethod have score and derivative evaluation
	/// evaluated both in the context of the whole Pose and in the context
	/// of residue or residue-pairs?  This covers scoring terms like env-smooth
	/// wherein the CBeta's get derivatives for increasing the neighbor counts
	/// for surrounding residues, and terms like constraints, which are definable
	/// on arbitrary number of residues (e.g. more than 2); both of these terms
	/// could be used in RTMin, and both should use the residue and residue-pair
	/// evaluation scheme with the MinimizationGraph for the majority of the
	/// work they do.  (Now, high-order constraints (3-body or above) will not
	/// be properly evaluated within RTMin.).  The default implementation
	/// returns "false".
	virtual
	bool
	defines_high_order_terms( pose::Pose const & ) const;

	/// @brief Evaluate the XYZ derivative for an atom in the pose.
	/// Called during the atomtree derivative calculation, atom_tree_minimize.cc,
	/// through the ScoreFunction::eval_atom_derivative intermediary.
	/// F1 and F2 should not zeroed, rather, this class should accumulate its contribution
	/// from this atom's XYZ derivative
	///
	/// @details The derivative scheme is based on that of Abe, Braun, Noguti and Go (1984)
	/// "Rapid Calculation of First and Second Derivatives of Conformational Energy with
	/// Respect to Dihedral Angles for Proteins. General Recurrent Equations"
	/// Computers & Chemistry 8(4) pp. 239-247. F1 and F2 correspond roughly to Fa and Ga,
	/// respectively, of equations 7a & 7b in that paper.
	virtual
	void
	eval_atom_derivative(
		id::AtomID const & id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		ScoreFunction const & sfxn,
		EnergyMap const & emap,
		Vector & F1,
		Vector & F2
	) const;


	/* /// @brief Evaluate the derivative for a DOF in the atom tree.
	virtual
	Real
	eval_dof_derivative(
	id::DOF_ID const &,
	id::TorsionID const &,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap const &
	) const;*/

	/// @brief called by the ScoreFunction at the end of energy evaluation.
	/// The derived class has the opportunity to accumulate a score
	/// into the pose's total_energy EnergyMap.  WholeStructure energies
	/// operate within this method; any method using a NeighborList during
	/// minimization would also operate within this function call.
	virtual
	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const & sfxn,
		EnergyMap & total_energy
	) const;

	/// @brief Returns the score types that this energy method computes.
	ScoreTypes const &
	score_types() const
	{
		return score_types_;
	}

	/// @brief Return one of the 7 kinds of energy methods that exist:
	/// e.g. context-dependent-one-body vs whole-structure.
	virtual
	EnergyMethodType
	method_type() const = 0;

	/// @brief Indicate in the context-graphs-required list which
	/// context-graphs this energy method requires that the Pose
	/// maintain when doing neighbor evaluation.  Context graphs are
	/// allowed
	virtual
	void indicate_required_context_graphs(
		utility::vector1< bool > & context_graphs_required
	) const = 0;

	/// @brief Return the version of the energy method
	virtual
	core::Size version() const = 0;

protected:

	void
	set_score_types( EnergyMethodCreatorOP creator );

private:
	ScoreTypes score_types_;


};

} // methods
} // scoring
} // core

#endif // INCLUDED_core_scoring_ScoreFunction_HH
