// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/constraint_generator/DistanceConstraintGenerator.hh
/// @brief Generates AtomPair constraints to enforce a given distance between two residue subsets
/// @author Tom Linsky (tlinsky@gmail.com)

#ifndef INCLUDED_protocols_constraint_generator_DistanceConstraintGenerator_hh
#define INCLUDED_protocols_constraint_generator_DistanceConstraintGenerator_hh

// Unit headers
#include <protocols/constraint_generator/DistanceConstraintGenerator.fwd.hh>
#include <protocols/constraint_generator/ConstraintGenerator.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace constraint_generator {

///@brief Generates AtomPair constraints to enforce a given distance between two residue subsets
class DistanceConstraintGenerator : public protocols::constraint_generator::ConstraintGenerator {
public:
	DistanceConstraintGenerator();

	virtual ~DistanceConstraintGenerator();

	static std::string
	class_name();

	protocols::constraint_generator::ConstraintGeneratorOP
	clone() const;

	virtual core::scoring::constraints::ConstraintCOPs
	apply( core::pose::Pose const & pose ) const;

protected:
	virtual void
	parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data );

public:
	void
	set_function( std::string const & func_def_str );

private:
	core::scoring::constraints::ConstraintOP
	create_constraint(
		core::pose::Pose const & pose,
		core::Size const resid1,
		core::Size const resid2 ) const;

private:
	core::scoring::func::FuncCOP func_;
	core::select::residue_selector::ResidueSelectorCOP selector1_;
	core::select::residue_selector::ResidueSelectorCOP selector2_;
	std::string atom_name1_;
	std::string atom_name2_;
};

} //protocols
} //constraint_generator

#endif //INCLUDED_protocols_constraint_generator_DistanceConstraintGenerator_hh

