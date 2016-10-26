// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/constraints/AmbiguousNMRConstraint.hh
/// @brief contains declarations for a type of constraint that holds multiple
/// other constrains that belong to each other and are all evaluate at once
/// @author Florian Richter (floric@u.washington.edu, march 2008)


#ifndef INCLUDED_core_scoring_constraints_AmbiguousNMRConstraint_hh
#define INCLUDED_core_scoring_constraints_AmbiguousNMRConstraint_hh

#include <core/scoring/constraints/AmbiguousNMRConstraint.fwd.hh>

// Unit header
#include <core/scoring/constraints/MultiConstraint.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/ScoreType.hh>

#include <core/id/AtomID.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/pose/Pose.fwd.hh>

//Utility Headers
#include <numeric/xyzVector.fwd.hh>

// STL Headers
#include <map>

#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {


class AmbiguousNMRConstraint : public MultiConstraint {
public:

	/// @brief default Constructor
	AmbiguousNMRConstraint( func::FuncOP func = NULL );

	/// @brief Constructor
	AmbiguousNMRConstraint( ConstraintCOPs const & cst_in, func::FuncOP func );

	virtual
	ConstraintOP clone() const;

	virtual
	ConstraintOP clone( func::FuncOP func ) const;

	virtual
	MultiConstraintOP empty_clone() const;

	virtual std::string type() const;

	virtual bool operator == ( Constraint const & other ) const;
	virtual bool same_type_as_me( Constraint const & other ) const;

	// Needed to get the base class overloads
	using Constraint::score;
	using Constraint::dist;

	/// @brief compute score
	virtual
	void
	score( func::XYZ_Func const & xyz_func, EnergyMap const & weights, EnergyMap & emap ) const;

	virtual
	core::Real
	dist( func::XYZ_Func const & xyz ) const;

	/// @brief add individual constraint into AmbiguousNMRConstraint
	virtual
	void
	add_individual_constraint( ConstraintCOP cst_in );

	virtual
	ConstraintOP
	remap_resid( core::id::SequenceMapping const &seqmap ) const;

	/// @brief compute atom deriv
	virtual
	void
	fill_f1_f2(
		AtomID const & atom,
		func::XYZ_Func const & xyz,
		Vector & F1,
		Vector & F2,
		EnergyMap const & weights
	) const;

	/// @brief Returns the func::Func object associated with this Constraint object.
	virtual
	func::Func const & get_func() const;

	//  virtual
	//  void show( std::ostream& out ) const;

	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// @details read definition of a multiconstraint. Since a MultiConstraint is essentially a vector of
	virtual void
	read_def(
		std::istream& data,
		core::pose::Pose const& pose,
		func::FuncFactory const & func_factory
	);

	virtual
	void show_def( std::ostream& out, pose::Pose const& pose ) const;

	virtual
	Size show_violations( std::ostream & out, pose::Pose const & pose, Size verbose_level, Real threshold = 1.0 ) const;

protected:
	/// @brief Explicit copy constructor so that derived classes will recieve a deep copy
	/// of the Func this class contains.
	AmbiguousNMRConstraint( AmbiguousNMRConstraint const & src );

private:
	func::FuncOP func_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //AmbiguousNMRConstraint

} //constraints
} //scoring
} //core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_constraints_AmbiguousNMRConstraint )
#endif // SERIALIZATION


#endif
