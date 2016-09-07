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

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace protocols {
namespace constraints_additional {

class COMCoordinateConstraint : public core::scoring::constraints::Constraint {
public:

	/// null constructor
	COMCoordinateConstraint( ) :
		core::scoring::constraints::Constraint( core::scoring::coordinate_constraint )
	{ }

	/// ctor from atom list + input pose
	COMCoordinateConstraint(
		utility::vector1< AtomID > const & atms,
		core::Vector const & COM_target,
		core::Real stdv,
		core::Real interval,
		core::scoring::ScoreType scoretype = core::scoring::coordinate_constraint
	);

	COMCoordinateConstraint(
		utility::vector1< AtomID > const & atms,
		core::Vector const & COM_target,
		core::Real stdv,
		core::scoring::ScoreType scoretype = core::scoring::coordinate_constraint
	);

	COMCoordinateConstraint(
		utility::vector1< AtomID > const & atms,
		core::Vector const & COM_target,
		core::scoring::ScoreType scoretype = core::scoring::coordinate_constraint
	);

	core::scoring::constraints::ConstraintOP clone() const override {
		return core::scoring::constraints::ConstraintOP( new COMCoordinateConstraint( atms_, COM_target_, stdv_, interval_ ) );
	}

	core::scoring::constraints::ConstraintOP remapped_clone(
		core::pose::Pose const& src,
		core::pose::Pose const& dest,
		core::id::SequenceMappingCOP smap
	) const override ;

	bool operator == ( core::scoring::constraints::Constraint const & other ) const override;

	bool same_type_as_me( core::scoring::constraints::Constraint const & other ) const override;



	void
	score( core::scoring::func::XYZ_Func const & xyz,
		core::scoring::EnergyMap const &, core::scoring::EnergyMap & emap ) const override;

	// do some pre-scoring calculations
	void setup_for_scoring( core::scoring::func::XYZ_Func const & xyz,
		core::scoring::ScoreFunction const &scfxn ) const override;

	// call the setup_for_derivatives for each constraint
	void setup_for_derivatives(  core::scoring::func::XYZ_Func const & xyz, core::scoring::ScoreFunction const &scfxn ) const override;

	// atom deriv

	void
	fill_f1_f2(
		AtomID const & atom,
		core::scoring::func::XYZ_Func const & xyz,
		core::Vector & F1,
		core::Vector & F2,
		core::scoring::EnergyMap const & weights
	) const override;

	std::string type() const override;


	Size
	natoms() const override;


	core::scoring::constraints::ConstraintOP
	remap_resid( core::id::SequenceMapping const & ) const override { return nullptr; }


	AtomID const &
	atom( Size const n ) const override;

	void show( std::ostream& out ) const override;

	void show_def( std::ostream& out, core::pose::Pose const & pose ) const override;
	void read_def( std::istream& in, core::pose::Pose const & pose, core::scoring::func::FuncFactory const & func_factory ) override;

	Size show_violations( std::ostream & out, core::pose::Pose const & pose, Size verbose_level, core::Real threshold = 1.0 ) const override;

	//protected:
	//void init( pose::Pose const & start_pose );

private:
	// data
	mutable core::Vector dCOM_; // to change value at setup_for_scoring...is there any better way than this?
	core::Vector COM_target_;
	utility::vector1< AtomID > atms_;
	core::Real stdv_;
	core::Real interval_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

}
}

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_constraints_additional_COMCoordinateConstraint )
#endif // SERIALIZATION


#endif
