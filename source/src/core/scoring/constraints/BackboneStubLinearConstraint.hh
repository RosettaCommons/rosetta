// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constraints/BackboneStubLinearConstraint.hh
///
/// @brief
/// @author Lei Shi


#ifndef INCLUDED_core_scoring_constraints_BackboneStubLinearConstraint_hh
#define INCLUDED_core_scoring_constraints_BackboneStubLinearConstraint_hh

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/BackboneStubLinearConstraint.fwd.hh>

#include <utility/pointer/owning_ptr.hh>

#include <core/conformation/Residue.fwd.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {


/// @brief This constraint favors the backbone landing on a "stub" backbone, which puts the sidechain in a pre-determined desirable location
///
class BackboneStubLinearConstraint : public Constraint
{
public:
	virtual std::string type() const {
		return "BackboneStubLinear";
	}

	BackboneStubLinearConstraint(
		pose::Pose const & pose,
		Size const seqpos,
		AtomID const & fixed_atom_id,
		conformation::Residue const & target_rsd,
		core::Real const & superposition_bonus,
		core::Real const & CB_force_constant
	);

	virtual ~BackboneStubLinearConstraint() {};

	virtual Size natoms() const { return atom_ids_.size(); };

	virtual AtomID const & atom( Size const index ) const { return atom_ids_[index]; };

	/// @brief possibility to compare constraint according to data
	/// and not just pointers
	virtual bool operator == ( Constraint const & other ) const;

	virtual bool same_type_as_me( Constraint const & other ) const;

	// Needed to get the base class overloads
	using Constraint::score;
	using Constraint::dist;

	virtual
	void
	score( func::XYZ_Func const & xyz_func, EnergyMap const & weights, EnergyMap & emap ) const;

	/// @details "Distance" for BackboneStubLinearConstraint isn't all that simple
	virtual
	core::Real
	dist( core::scoring::func::XYZ_Func const & ) const { return 0; }

	virtual
	void
	fill_f1_f2(
		AtomID const & atom,
		func::XYZ_Func const & xyz,
		Vector & F1,
		Vector & F2,
		EnergyMap const & weights
	) const;

	virtual void show( std::ostream& out ) const;

	/// @brief returns the private member seqpos_
	core::Size seqpos() const;
	virtual
	ConstraintOP clone() const;

	/// @brief Copies the data from this Constraint into a new object and returns an OP
	/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
	/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
	/// to the new object. Intended to be implemented by derived classes.
	virtual ConstraintOP remapped_clone(
		pose::Pose const & src,
		pose::Pose const & dest,
		id::SequenceMappingCOP map = NULL
	) const;


	/*virtual
	ConstraintOP
	remap_resid( core::id::SequenceMapping const & seqmap ) const;
	*/

private:

	core::Real superposition_bonus_;
	core::Real CB_force_constant_;

	Size seqpos_;
	AtomID CB_atom_id_, CA_atom_id_, C_atom_id_, N_atom_id_;
	utility::vector1< AtomID > atom_ids_;

	core::Vector CB_target_, CA_target_, C_target_, N_target_, CB_CA_target_, C_N_target_;

	AtomID fixed_atom_id_;
	core::Vector fixed_reference_point_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	BackboneStubLinearConstraint();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // BackboneStubLinearConstraint


} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_constraints_BackboneStubLinearConstraint )
#endif // SERIALIZATION


#endif // INCLUDED_core_scoring_constraints_BackboneStubLinearConstraint_HH
