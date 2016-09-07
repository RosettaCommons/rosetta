// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


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



	core::scoring::constraints::ConstraintOP clone() const override {
		return core::scoring::constraints::ConstraintOP( new AmbiguousMultiConstraint(*this) );
	}

	bool operator == ( core::scoring::constraints::Constraint const & other ) const override;

	bool same_type_as_me( core::scoring::constraints::Constraint const & other ) const override;


	std::string type() const override {
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
		core::scoring::EnergyMap & emap ) const override;



	core::scoring::constraints::ConstraintOP
	remap_resid( core::id::SequenceMapping const &seqmap ) const override;

	/// @brief compute atom deriv
	void
	fill_f1_f2(
		core::id::AtomID const & atom,
		core::scoring::func::XYZ_Func const & xyz,
		core::Vector & F1,
		core::Vector & F2,
		core::scoring::EnergyMap const & weights
	) const override;

	void show( std::ostream& out) const override;

	// void read_def( std::istream& in, pose::Pose const& pose, FuncFactory const& func_factory );

	Size show_violations( std::ostream& out, core::pose::Pose const& pose, core::Size verbose_level, core::Real threshold = 1.0 ) const override;


private:
	core::Size num_active_constraints_;

	mutable core::scoring::constraints::ConstraintCOPs active_constraints_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	AmbiguousMultiConstraint();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //AmbiguousMultiConstraint

}
}

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_constraints_additional_AmbiguousMultiConstraint )
#endif // SERIALIZATION


#endif
