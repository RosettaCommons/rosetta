// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constraints/BackboneStubConstraint.hh
///
/// @brief
/// @author John Karanicolas, Sarel Fleishman


#ifndef INCLUDED_core_scoring_constraints_BackboneStubConstraint_hh
#define INCLUDED_core_scoring_constraints_BackboneStubConstraint_hh

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/BackboneStubConstraint.fwd.hh>

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
class BackboneStubConstraint : public Constraint
{
public:
	virtual std::string type() const {
		return "BackboneStub";
	}

	BackboneStubConstraint(
		pose::Pose const & pose,
		Size const seqpos,
		AtomID const & fixed_atom_id,
		pose::Pose const & target_pose,
		Size const target_seqpos,
		core::Real const & superposition_bonus,
		core::Real const & force_constant
	);


	BackboneStubConstraint(
		pose::Pose const & pose,
		Size const seqpos,
		AtomID const & fixed_atom_id,
		conformation::Residue const & target_rsd,
		core::Real const & superposition_bonus,
		core::Real const & force_constant,
		std::string const & sidechain_atom_name,
		std::string const & primary_backbone_atom_name,
		std::string const & backbone_atom1_name,
		std::string const & backbone_atom2_name
	);

	BackboneStubConstraint(
		pose::Pose const & pose,
		Size const seqpos,
		AtomID const & fixed_atom_id,
		conformation::Residue const & target_rsd,
		core::Real const & superposition_bonus,
		core::Real const & force_constant
	);

	/// @brief Copy constructor -- performs a deep copy
	BackboneStubConstraint( BackboneStubConstraint const & src );

	virtual ~BackboneStubConstraint() {};

	virtual Size natoms() const { return atom_ids_.size(); };

	virtual AtomID const & atom( Size const index ) const { return atom_ids_[index]; };

	/// @brief possibility to compare constraint according to data
	/// and not just pointers
	bool operator == ( Constraint const & other ) const;

	bool same_type_as_me( Constraint const & other ) const;

	virtual
	void
	score( func::XYZ_Func const & xyz_func, EnergyMap const & weights, EnergyMap & emap ) const;

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
	core::Real force_constant_;

	Size seqpos_;
	AtomID sidechain_atom_id_, primary_atom_id_, backbone_atom2_id_, backbone_atom1_id_;
	utility::vector1< AtomID > atom_ids_;

	core::Vector sidechain_target_, primary_backbone_target_, backbone2_target_, backbone1_target_, primary_sidechain_target_, backbone1_2_target_;

	AtomID fixed_atom_id_;
	core::Vector fixed_reference_point_;

	utility::pointer::shared_ptr< AngleConstraint > ang_cst_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	BackboneStubConstraint();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // BackboneStubConstraint


} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_constraints_BackboneStubConstraint )
#endif // SERIALIZATION


#endif // INCLUDED_core_scoring_constraints_BackboneStubConstraint_HH
