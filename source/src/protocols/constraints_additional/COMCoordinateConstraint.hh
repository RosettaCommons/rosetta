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

#ifndef INCLUDED_protocols_constraints_additional_COMCoordinateConstraint_hh
#define INCLUDED_protocols_constraints_additional_COMCoordinateConstraint_hh

// Package headers
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/XYZ_Func.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/func/FuncFactory.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>

#include <utility/vector1.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>

namespace protocols {
namespace constraints_additional {

using namespace core;

class COMCoordinateConstraint : public scoring::constraints::Constraint {
public:

	/// null constructor
	COMCoordinateConstraint( ) :
		scoring::constraints::Constraint( scoring::coordinate_constraint )
	{ }

	/// ctor from atom list + input pose
	COMCoordinateConstraint(
	  utility::vector1< AtomID > const & atms,
		Vector const & COM_target,
		Real stdv,
		Real interval,
		scoring::ScoreType scoretype = scoring::coordinate_constraint
	);

	COMCoordinateConstraint(
	  utility::vector1< AtomID > const & atms,
		Vector const & COM_target,
		Real stdv,
		scoring::ScoreType scoretype = scoring::coordinate_constraint
	);

	COMCoordinateConstraint(
	  utility::vector1< AtomID > const & atms,
		Vector const & COM_target,
		scoring::ScoreType scoretype = scoring::coordinate_constraint
	);

	virtual scoring::constraints::ConstraintOP clone() const {
		return scoring::constraints::ConstraintOP( new COMCoordinateConstraint( atms_, COM_target_, stdv_, interval_ ) );
	}

	virtual scoring::constraints::ConstraintOP remapped_clone(
		pose::Pose const& src,
	  pose::Pose const& dest,
	  id::SequenceMappingCOP smap
		) const ;


	void
	score( scoring::func::XYZ_Func const & xyz,
				 scoring::EnergyMap const &, scoring::EnergyMap & emap ) const;

	// do some pre-scoring calculations
	void setup_for_scoring( scoring::func::XYZ_Func const & xyz,
													scoring::ScoreFunction const &scfxn ) const;

	// call the setup_for_derivatives for each constraint
	void setup_for_derivatives(  scoring::func::XYZ_Func const & xyz, scoring::ScoreFunction const &scfxn ) const;

	// atom deriv
	virtual
	void
	fill_f1_f2(
		AtomID const & atom,
		scoring::func::XYZ_Func const & xyz,
		Vector & F1,
		Vector & F2,
		scoring::EnergyMap const & weights
	) const;

	std::string type() const;


	Size
	natoms() const;

	virtual
	scoring::constraints::ConstraintOP
	remap_resid( id::SequenceMapping const & ) const { return NULL; }


	AtomID const &
	atom( Size const n ) const;

	void show( std::ostream& out ) const;

	void show_def( std::ostream& out, pose::Pose const & pose ) const;
	void read_def( std::istream& in, pose::Pose const & pose, scoring::func::FuncFactory const & func_factory );

	Size show_violations( std::ostream & out, pose::Pose const & pose, Size verbose_level, Real threshold = 1.0 ) const;

	//protected:
	//void init( pose::Pose const & start_pose );

private:
	// data
	mutable Vector dCOM_; // to change value at setup_for_scoring...is there any better way than this?
	Vector COM_target_;
	utility::vector1< AtomID > atms_;
	Real stdv_;
	Real interval_;
};

}
}

#endif
