// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/constraint_generator/DihedralConstraintGenerator.cc
/// @brief A Constraint Generator for dihedral angles.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <protocols/constraint_generator/DihedralConstraintGenerator.hh>
#include <protocols/constraint_generator/DihedralConstraintGeneratorCreator.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.constraint_generator.DihedralConstraintGenerator" );

namespace protocols {
namespace constraint_generator {

protocols::constraint_generator::ConstraintGeneratorOP
DihedralConstraintGeneratorCreator::create_constraint_generator() const
{
	return protocols::constraint_generator::ConstraintGeneratorOP( new DihedralConstraintGenerator );
}

std::string
DihedralConstraintGeneratorCreator::keyname() const
{
	return DihedralConstraintGeneratorCreator::constraint_generator_name();
}

std::string
DihedralConstraintGeneratorCreator::constraint_generator_name()
{
	return "DihedralConstraintGenerator";
}

DihedralConstraintGenerator::DihedralConstraintGenerator():
	protocols::constraint_generator::ConstraintGenerator( DihedralConstraintGeneratorCreator::constraint_generator_name() )
{
}

DihedralConstraintGenerator::~DihedralConstraintGenerator()
{}

protocols::constraint_generator::ConstraintGeneratorOP
DihedralConstraintGenerator::clone() const
{
	return protocols::constraint_generator::ConstraintGeneratorOP( new DihedralConstraintGenerator( *this ) );
}

void
DihedralConstraintGenerator::parse_tag( utility::tag::TagCOP , basic::datacache::DataMap & )
{
}

core::scoring::constraints::ConstraintCOPs
DihedralConstraintGenerator::apply( core::pose::Pose const & ) const
{
	return core::scoring::constraints::ConstraintCOPs();
}

} //protocols
} //constraint_generator






