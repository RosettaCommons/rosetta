// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/constraint_generator/DihedralConstraintGenerator.hh
/// @brief A Constraint Generator for dihedral angles.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_constraint_generator_DihedralConstraintGenerator_hh
#define INCLUDED_protocols_constraint_generator_DihedralConstraintGenerator_hh

// Unit headers
#include <protocols/constraint_generator/DihedralConstraintGenerator.fwd.hh>
#include <protocols/constraint_generator/ConstraintGenerator.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace constraint_generator {

///@brief A Constraint Generator for dihedral angles.
class DihedralConstraintGenerator : public protocols::constraint_generator::ConstraintGenerator {
public:
	DihedralConstraintGenerator();

	virtual ~DihedralConstraintGenerator();

	protocols::constraint_generator::ConstraintGeneratorOP
	clone() const;

	virtual core::scoring::constraints::ConstraintCOPs
	apply( core::pose::Pose const & pose ) const;

protected:
	virtual void
	parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data );

private:

};

} //protocols
} //constraint_generator

#endif //INCLUDED_protocols_constraint_generator_DihedralConstraintGenerator_hh

