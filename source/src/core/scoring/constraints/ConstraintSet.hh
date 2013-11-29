// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


#ifndef INCLUDED_core_scoring_constraints_ConstraintSet_hh
#define INCLUDED_core_scoring_constraints_ConstraintSet_hh

// Unit headers
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/func/Func.hh> /// cant get to compile w/ Func.fwd.hh


// Package headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/MinimizationData.fwd.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/Constraints.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/constraints/DOF_Constraint.hh>
#include <core/scoring/ScoreType.hh>

// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
// AUTO-REMOVED #include <core/id/TorsionID.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/conformation/signals/LengthEvent.fwd.hh> //for observing conformation length changes
#include <core/conformation/signals/ConnectionEvent.fwd.hh>
// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

//Utility Headers

// C++ Headers
#include <map>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace constraints {

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
/// silly helper class, a wrapper for std::map so we can hold in owning_ptr
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

class ResidueConstraints : public utility::pointer::ReferenceCount {
public:
	///@brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~ResidueConstraints();
	typedef std::map< Size, ConstraintsOP > Map;
	typedef Map::const_iterator const_iterator;
	typedef Map::iterator iterator;

public:
	const_iterator
	begin() const
	{
		return map_.begin();
	}

	const_iterator
	end() const
	{
		return map_.end();
	}

	const_iterator
	find( Size const seqpos ) const
	{
		return map_.find( seqpos );
	}

	iterator
	find( Size const seqpos )
	{
		return map_.find( seqpos );
	}

	void
	erase( Size const seqpos )
	{
		map_.erase( seqpos );
	}

	bool
	has( Size const seqpos )
	{
		return ( map_.find( seqpos ) != map_.end() );
	}

	void
	insert( Size const seqpos, ConstraintsOP cst )
	{
		map_.insert( std::make_pair( seqpos, cst ) );
	}

	Size
	size() const
	{
		return map_.size();
	}

	void
	clear()
	{
		map_.clear();
	}

private:
	Map map_;
};



////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
// ConstraintSet
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

class ConstraintSet : public utility::pointer::ReferenceCount {
public:
	typedef id::AtomID AtomID;
	typedef id::DOF_ID DOF_ID;

	typedef conformation::Residue Residue;
	typedef pose::Pose Pose;

	typedef utility::vector1< ResidueConstraintsOP > ResiduePairConstraints;
	typedef ResidueConstraints::const_iterator ResiduePairConstraintsIterator;

	//	typedef utility::vector1< DOF_ConstraintOP > DOF_Constraints;

public:
	ConstraintSet();
	ConstraintSet( ConstraintSet const & other );
	ConstraintSet(
		ConstraintSet const & other,
		core::Size start_residue,
		core::Size end_residue
	);
	/// @brief Destructor so far only detaches from conformation
	virtual ~ConstraintSet() { this->detach_from_conformation(); }

	virtual ConstraintSetOP clone() const;

	/// @brief Copies the data from this ConstraintSet into a new object and
	/// returns its OP; atoms are mapped to atoms with the same name in dest pose
	/// ( e.g. for switch from centroid to fullatom ) if a sequence_mapping is
	/// present it is used to map residue numbers .. NULL = identity mapping to
	/// the new object. This will really clone all constraints since they have to
	/// change their atom-numbers and residue-numbers
	virtual ConstraintSetOP
	remapped_clone(
		pose::Pose const& src,
		pose::Pose const& dest,
		id::SequenceMappingCOP smap = NULL
	) const;

	/// @brief  like remapped_clone, but constraints also steal_def from src-pose
	/// use, e.g., to get a new set of CoordinateConstraints for given xyz
	/// coordinates in src-pose
	virtual ConstraintSetOP steal_def_clone(
		pose::Pose const& src,
		pose::Pose const& dest,
		id::SequenceMappingCOP smap = NULL
	) const;


	/// @brief remaps the constraints in this particular constraint set according
	/// to brief the passed in sequence mapping --- redundant with
	/// remapped_clone!!!
	void
	remap_residue_positions(
		id::SequenceMapping const & smap
	);

	/// @brief Cache the ConstraintsCOP for a particular residue
	/// in the res_data_cache, under the "cst_res_data" element of
	/// the data-cache's CachableDataOP array.  Derived ConstraintSet
	/// classes should decide HERE whether any constraints should be evaluated
	/// for this residue, since, once a ConstraintsOP is added
	/// to the minimization graph, the derived class looses the chance to
	/// veto the evaluation of that constraint.
	virtual
	void
	setup_for_minimizing_for_residue(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		kinematics::MinimizerMapBase const & minmap,
		ResSingleMinimizationData & res_data_cache
	) const;


	/// @brief Cache the ConstraintsCOP for a particular residue pair
	/// in the respair_data_cache, under the "cst_respair_data" element of
	/// the data-cache's CachableDataOP array.  Derived ConstraintSet
	/// classes should decide HERE whether any constraints should be evaluated
	/// between this pair of residues, since, once a ConstraintsOP is added
	/// to the minimization graph, the derived class looses the chance to
	/// veto the evaluation of that constraint.
	virtual
	void
	setup_for_minimizing_for_residue_pair(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		kinematics::MinimizerMapBase const & minmap,
		ResSingleMinimizationData const & res1_data_cache,
		ResSingleMinimizationData const & res2_data_cache,
		ResPairMinimizationData & respair_data_cache
	) const;


 	virtual void
 	setup_for_scoring( pose::Pose & pose, ScoreFunction const & scfxn ) const;

 	virtual void
 	setup_for_derivatives( pose::Pose & pose, ScoreFunction const & scfxn ) const;

	///
	virtual
	void
	residue_pair_energy(
		Residue const & rsd1,
		Residue const & rsd2,
		Pose const & pose,
		ScoreFunction const & scorefxn,
		EnergyMap & emap
	) const;

	/// @brief Switching over to a pairwise decomposable eval-atom-deriv system for
	/// RTMin means deprecating the old "evaluate an atom's derivative wrt the entire structure"
	/// This function is preserved (for now) for use by the RNA_TorsionEnergy
 	virtual
	void
 	deprecated_eval_atom_derivative(
 		id::AtomID const & atom_id,
 		pose::Pose const & pose,
 		ScoreFunction const &,
 		EnergyMap const & weights,
 		Vector & F1,
 		Vector & F2
 	) const;

	/// @brief evaluate the derivatives for an atom that contains 3- or higher-body
	/// constraints.  Such derivatives cannot be evalauated in an extra-posal context
	/// (e.g. such as in RTMin).
 	virtual
	void
 	eval_multibody_atom_derivative(
 		id::AtomID const & atom_id,
 		pose::Pose const & pose,
 		ScoreFunction const &,
 		EnergyMap const & weights,
 		Vector & F1,
 		Vector & F2
 	) const;

	/// uses the dof constraints
	/*
	Real
	eval_dof_derivative(
		id::DOF_ID const & id,
		id::TorsionID const & tor,
		pose::Pose const & pose,
		ScoreFunction const & scorefxn,
		EnergyMap const & weights
	) const;*/


	///
	virtual bool
	residue_pair_constraint_exists( int const pos1, int const pos2 ) const
	{
		return ( residue_pair_constraints_.size() >= Size(pos1) &&
			residue_pair_constraints_[ pos1 ] &&
			residue_pair_constraints_[ pos1 ]->has( pos2 ) );
	}

	///
	virtual bool
	residue_pair_constraints_exists( Size const pos ) const {
		return ( (residue_pair_constraints_.size() >= pos ) &&
							residue_pair_constraints_[ pos ]
		);
	}

	///
	virtual void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	///
	void
	eval_intrares_energy(
											 conformation::Residue const & rsd,
											 EnergyMap & emap
											 ) const;

	/// Does *NOT* zero the emap values, just adds the additional contribution to the
	/// existing emap energies (so can be called inside finalize_total_energies)
	virtual void
	eval_non_residue_pair_energy(
		Pose const & pose,
 		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	///
	void
	add_constraint( ConstraintCOP cst );

	void
	add_constraints( ConstraintCOPs cst_list );

	///@brief add another constraint set to this constraint set
	void
	add_constraints( ConstraintSetCOP const cst_set );

	/// @brief Returns true if the constraint was successfully found and removed.
	/// if object comparison is set to true, the constraint to be removed is found
	/// through the Constraint::== operator and not through pointer comparison
	bool
	remove_constraint(
		ConstraintCOP cst,
		bool object_comparison );

	// @brief returns true if all the constraints in the list were successfully
	// found and removed.
	bool
	remove_constraints(
		ConstraintCOPs cst_list,
		bool object_comparison );

	/// @brief  Note -- still hacky. Will not be included in packing, just scoring
	/// and minimization
	void
	add_dof_constraint(
		DOF_ID const & id,
		FuncOP func,
		ScoreType const & t = dof_constraint
	);

	/// @brief Returns all constraints in the set as a flat list, regardless of
	/// type.
	ConstraintCOPs
	get_all_constraints() const;

	ResiduePairConstraintsIterator
	residue_pair_constraints_begin( Size resid ) const;

	ResiduePairConstraintsIterator
	residue_pair_constraints_end( Size resid ) const;

	void
	on_length_change( conformation::signals::LengthEvent const & event );

	void
	on_connection_change(
		core::conformation::signals::ConnectionEvent const & event
	);

	void
	attach_to_conformation( core::conformation::Conformation * conformation );

	void
	detach_from_conformation();

	Size
	revision_id() const;

	// output statements, useful for debugging
	virtual void
	show( std::ostream& out ) const;

	// output statements, useful for debugging
	virtual void
	show_definition( std::ostream& out, core::pose::Pose const& ) const;

	// this is the worst function name ever. This might as well be print_string.
	virtual void
	show_numbers( std::ostream& out ) const;

	virtual Size
	show_violations(
		std::ostream & out,
		pose::Pose &,
		Size verbose_level,
		Real threshold = 1
	) const;

	bool
	has_residue_pair_constraints() const {
		return residue_pair_constraints_.size() > 0;
	}

	bool
	has_intra_residue_constraints() const {
		return intra_residue_constraints_.size() > 0;
	}

	bool
	has_dof_constraints() const {
		return dof_constraints_.size() > 0;
	}

	bool
	has_non_residue_pair_constraints() const {
		return non_residue_pair_constraints_.size() > 0;
	}

	bool
	has_constraints() const {
		return has_residue_pair_constraints()
			|| has_intra_residue_constraints()
			|| has_dof_constraints()
			|| has_non_residue_pair_constraints();
	}

	void
	clear();

	bool
	is_empty() const;

protected:


	virtual
	void
	deprecated_eval_atom_derivative_for_residue_pairs(
		id::AtomID const & atom_id,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;

	void
	mark_revision_id_expired();

private:
	void
	add_residue_pair_constraint(
		Size const pos1,
		Size const pos2,
		ConstraintCOP cst
	);

	/// @brief Returns true iff the constraint was successfully found and removed.
	bool
	remove_residue_pair_constraint(
		Size const pos1,
		Size const pos2,
		ConstraintCOP cst,
		bool object_comparison
	);

protected:
	ResiduePairConstraints const &
	residue_pair_constraints() const;

	Constraints const& non_residue_pair_constraints() const {
		return non_residue_pair_constraints_;
	}

private:

	/// constraints are added symmetrically
	ResiduePairConstraints residue_pair_constraints_;

	ResidueConstraints intra_residue_constraints_;

	Constraints non_residue_pair_constraints_;

	// do not put anything in this -- used to return iterators to empty containers
	ResidueConstraints empty_rsdcst_;

	// vector of DOF constraints
	DOF_ConstraintOPs dof_constraints_;

	// a constraint set's "state identifier" for rapidly determining if any constraints
	// have changed since the last time the constraint energies were evaluated
	mutable Size revision_id_;
	mutable bool revision_id_current_;

	core::conformation::Conformation * conformation_pt_;
};

std::ostream & operator << (std::ostream & os, ConstraintSet const & set);

} // constraints
} // scoring
} // core

#endif
