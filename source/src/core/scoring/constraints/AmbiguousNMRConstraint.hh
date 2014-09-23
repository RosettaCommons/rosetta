// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
// AUTO-REMOVED #include <core/scoring/func/XYZ_Func.hh>
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


namespace core {
namespace scoring {
namespace constraints {


class AmbiguousNMRConstraint : public MultiConstraint {
public:

	/// @brief default Constructor
	AmbiguousNMRConstraint( func::FuncOP func = NULL );

	/// @brief Constructor
	AmbiguousNMRConstraint( const ConstraintCOPs & cst_in, func::FuncOP func );

// 	/// @brief
// 	void
// 	init_cst_score_types();

	///
	virtual
	ConstraintOP clone() const {
		if ( size() > 0 ) {
			return ConstraintOP( new AmbiguousNMRConstraint( member_constraints_, func_ ) );
		} else {
			return ConstraintOP( new AmbiguousNMRConstraint( func_ ) );
		}
	}

	///
	virtual
	ConstraintOP clone( func::FuncOP func ) const {
		if ( size() > 0 ) {
			return ConstraintOP( new AmbiguousNMRConstraint( member_constraints_, func ) );
		} else {
			return ConstraintOP( new AmbiguousNMRConstraint( func ) );
		}
	}

	virtual
	MultiConstraintOP empty_clone() const {
		return MultiConstraintOP( new AmbiguousNMRConstraint( get_func().clone() ) );
	}

	virtual std::string type() const {
		return "AmbiguousNMRConstraint";
	}


	/// @brief compute score
	virtual
	void
	score( func::XYZ_Func const & xyz_func, EnergyMap const & weights, EnergyMap & emap ) const;


	virtual
	core::Real
	dist( core::pose::Pose const & pose ) const;

	virtual
	core::Real
	dist( func::XYZ_Func const & xyz ) const;

	///@brief add individual constraint into AmbiguousNMRConstraint
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
	func::Func const & get_func() const {
		runtime_assert( func_ != 0 );
		return *func_;
	}


	virtual Real score( pose::Pose const& pose ) const {
		return func_->func( dist( pose ) );
	}


// 	virtual
// 	void show( std::ostream& out ) const;

////////////////////////////////////////////////////////////////////////////////////////////////////
///@details read definition of a multiconstraint. Since a MultiConstraint is essentially a vector of
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

private:
	func::FuncOP func_;

}; //AmbiguousNMRConstraint

} //constraints
} //scoring
} //core

#endif
