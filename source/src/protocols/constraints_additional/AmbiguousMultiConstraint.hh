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
/// @author Florian Richter (floric@u.washington.edu, april 2009)


#ifndef INCLUDED_protocols_constraints_additional_AmbiguousMultiConstraint_hh
#define INCLUDED_protocols_constraints_additional_AmbiguousMultiConstraint_hh

// Unit header
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <protocols/constraints_additional/AmbiguousMultiConstraint.fwd.hh>
#include <core/scoring/func/XYZ_Func.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/id/SequenceMapping.fwd.hh>

#include <core/id/AtomID.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/pose/Pose.fwd.hh>

//Utility Headers
#include <numeric/xyzVector.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace constraints_additional {


class AmbiguousMultiConstraint : public core::scoring::constraints::AmbiguousConstraint {
public:

  /// @brief Constructor
	AmbiguousMultiConstraint(
		core::Size num_act_csts);

	/// @brief Constructor
  AmbiguousMultiConstraint(
		core::Size num_act_csts,
		core::scoring::constraints::ConstraintCOPs & cst_in ) ;


	virtual
	core::scoring::constraints::ConstraintOP clone() const {
		return core::scoring::constraints::ConstraintOP( new AmbiguousMultiConstraint(*this) );
	}

	std::string type() const {
		return "AmbiguousMultiConstraint";
	}

	/// @brief read in constraint defiinition
	//void
	//read_def( std::istream& data, pose::Pose const& pose, FuncFactory const& func_factory );

  /// @brief compute score
  void
  score(
		core::scoring::func::XYZ_Func const & xyz_func,
		core::scoring::EnergyMap const & weights,
		core::scoring::EnergyMap & emap ) const;


	virtual
	core::scoring::constraints::ConstraintOP
	remap_resid( core::id::SequenceMapping const &seqmap ) const;

	/// @brief compute atom deriv
	void
	fill_f1_f2(
		core::id::AtomID const & atom,
		core::scoring::func::XYZ_Func const & xyz,
		core::Vector & F1,
		core::Vector & F2,
		core::scoring::EnergyMap const & weights
	) const;

  void show( std::ostream& out) const;

	//	void read_def( std::istream& in, pose::Pose const& pose, FuncFactory const& func_factory );

  Size show_violations( std::ostream& out, core::pose::Pose const& pose, core::Size verbose_level, core::Real threshold = 1.0 ) const;


private:
	core::Size num_active_constraints_;

	mutable core::scoring::constraints::ConstraintCOPs active_constraints_;

}; //AmbiguousMultiConstraint

}
}

#endif
