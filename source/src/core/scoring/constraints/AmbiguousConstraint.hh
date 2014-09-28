// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief contains declarations for a type of constraint that holds multiple other
/// @brief constrains where only the one with  lowest energy is considered
/// @brief
/// @author Florian Richter (floric@u.washington.edu, march 2008)


#ifndef INCLUDED_core_scoring_constraints_AmbiguousConstraint_hh
#define INCLUDED_core_scoring_constraints_AmbiguousConstraint_hh

// Unit header
#include <core/scoring/constraints/MultiConstraint.hh>
#include <core/scoring/constraints/AmbiguousConstraint.fwd.hh>
// AUTO-REMOVED #include <core/scoring/func/XYZ_Func.hh>
#include <core/scoring/func/XYZ_Func.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/func/FuncFactory.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>

#include <core/id/AtomID.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/pose/Pose.fwd.hh>


//Utility Headers
#include <numeric/xyzVector.fwd.hh>

#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>



namespace core {
namespace scoring {
namespace constraints {

///@brief Nested constraint  where only the one with  lowest energy is considered.
class AmbiguousConstraint : public MultiConstraint {
public:

	/// @brief Constructor
	AmbiguousConstraint();

	/// @brief Constructor
	AmbiguousConstraint( ConstraintCOPs & cst_in ) ;

	///
	virtual
	ConstraintOP clone() const {
		return ConstraintOP( new AmbiguousConstraint(*this) );
	}

	virtual
	MultiConstraintOP empty_clone() const {
		return MultiConstraintOP( new AmbiguousConstraint );
	}

	/// @brief
	void
	init_cst_score_types();

	std::string type() const {
		return "AmbiguousConstraint";
	}

	/// @brief read in constraint defiinition
	void
	read_def( std::istream& data, pose::Pose const& pose,func::FuncFactory const& func_factory );

	/// @brief possibility to compare constraint according to data
	/// and not just pointers
	bool operator == ( Constraint const & other ) const;

	/// @brief compute score
	void
	score( func::XYZ_Func const & xyz_func, EnergyMap const & weights, EnergyMap & emap ) const;

	core::Real
	calculate_total_cst_score( EnergyMap const & weights, EnergyMap & emap) const;

	virtual
	ConstraintOP
	remap_resid( core::id::SequenceMapping const &seqmap ) const;

	/// @brief compute atom deriv
	void
	fill_f1_f2(
		AtomID const & atom,
		func::XYZ_Func const & xyz,
		Vector & F1,
		Vector & F2,
		EnergyMap const & weights
	) const;

	void show( std::ostream& out) const;

	ConstraintCOP active_constraint() const;

	//	void read_def( std::istream& in, pose::Pose const& pose,func::FuncFactory const& func_factory );

	Size show_violations( std::ostream& out, pose::Pose const& pose, Size verbose_level, Real threshold = 1.0 ) const;



private:
	mutable ConstraintCOP active_constraint_;
	mutable core::Real low_total_cst_score_;
	mutable EnergyMap low_EMap_;
	mutable EnergyMap temp_EMap_;
	ScoreTypes cst_score_types_;

}; //AmbiguousConstraint

} //constraints
} //scoring
} //core

#endif
