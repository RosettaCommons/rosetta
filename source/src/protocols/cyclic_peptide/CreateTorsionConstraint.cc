// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
///
/// @file   protocols/cyclic_peptide/CreateTorsionConstraint.cc
/// @brief  Add torsion constraints to the current pose conformation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#include <protocols/cyclic_peptide/CreateTorsionConstraint.hh>
#include <protocols/cyclic_peptide/CreateTorsionConstraintCreator.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/kinematics/FoldTree.hh>

#include <protocols/loops/loops_main.hh>
#include <protocols/loops/util.hh>

#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.cyclic_peptide.CreateTorsionConstraint" );

namespace protocols {
namespace cyclic_peptide {

CreateTorsionConstraint::CreateTorsionConstraint() //:
{}
CreateTorsionConstraint::~CreateTorsionConstraint(){}

void CreateTorsionConstraint::apply( core::pose::Pose & /*pose*/ )
{
	//TODO
	return;
}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
CreateTorsionConstraint::parse_my_tag(
	TagCOP /*tag*/,
	basic::datacache::DataMap &,
	Filters_map const &,
	moves::Movers_map const &,
	Pose const &
)
{

	//TODO
	return;
}
	
moves::MoverOP CreateTorsionConstraint::clone() const { return new CreateTorsionConstraint( *this ); }
moves::MoverOP CreateTorsionConstraint::fresh_instance() const { return new CreateTorsionConstraint; }

protocols::moves::MoverOP
CreateTorsionConstraintCreator::create_mover() const {
	return new CreateTorsionConstraint;
}

std::string
CreateTorsionConstraintCreator::keyname() const
{
	return CreateTorsionConstraintCreator::mover_name();
}

std::string
CreateTorsionConstraintCreator::mover_name()
{
	return "CreateTorsionConstraint";
}

std::string
CreateTorsionConstraint::get_name() const {
	return "CreateTorsionConstraint";
}
	
} // moves
} // protocols
