// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief

#ifndef INCLUDED_core_scoring_constraints_NamedAngleConstraint_hh
#define INCLUDED_core_scoring_constraints_NamedAngleConstraint_hh

//Unit headers
#include <core/scoring/constraints/NamedAngleConstraint.fwd.hh>
#include <core/scoring/constraints/AngleConstraint.hh>

//Core headers
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>

//Utility headers
#include <utility/vector1.hh>

// C++ Headers

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {

class NamedAngleConstraint : public AngleConstraint {
public:
	NamedAngleConstraint(
		id::NamedAtomID const & a1,
		id::NamedAtomID const & a2,
		id::NamedAtomID const & a3,
		func::FuncOP func,
		ScoreType scoretype = angle_constraint );

	virtual std::string type() const override;

	virtual ConstraintOP clone() const override;

	/// @brief Copies the data from this Constraint into a new object and returns an OP
	/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
	/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
	/// to the new object. Intended to be implemented by derived classes.
	virtual ConstraintOP remapped_clone(
		pose::Pose const & src,
		pose::Pose const & dest,
		id::SequenceMappingCOP map = NULL ) const override;

	/// @brief This overrride updates the sequence numbering but not the atom names.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	virtual ConstraintOP remap_resid( core::id::SequenceMapping const &seqmap ) const override;

	virtual bool operator == ( Constraint const & rhs ) const override;
	virtual bool same_type_as_me( Constraint const & other ) const override;

	virtual void show_def( std::ostream & out, pose::Pose const & pose ) const override;
	void show_def_nopose( std::ostream & out ) const;

	virtual void read_def(
		std::istream & in,
		pose::Pose const & pose,
		func::FuncFactory const & func_factory ) override;

	// Needed to get the base class overloads
	using Constraint::score;
	using Constraint::dist;

	virtual void score( func::XYZ_Func const & xyz, EnergyMap const &, EnergyMap & emap ) const override;

	//fpd use the same machinery as NamedAtomPairConstraint where names are mapped to indices here
	virtual void setup_for_scoring( func::XYZ_Func const &, ScoreFunction const & ) const override;


	virtual
	core::Real
	dist( core::scoring::func::XYZ_Func const & xyz ) const override;

protected:
	NamedAngleConstraint( NamedAngleConstraint const & src );


private:
	id::NamedAtomID named_atom1_;
	id::NamedAtomID named_atom2_;
	id::NamedAtomID named_atom3_;

	core::Size type1_id_;
	core::Size type2_id_;
	core::Size type3_id_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	NamedAngleConstraint();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


}
}
}

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_constraints_NamedAngleConstraint )
#endif // SERIALIZATION


#endif
