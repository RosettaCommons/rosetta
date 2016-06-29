// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/constraint_generator/DistanceConstraintGenerator.cc
/// @brief Generates AtomPair constraints to enforce a given distance between two residue subsets
/// @author Tom Linsky (tlinsky@gmail.com)

// Unit headers
#include <protocols/constraint_generator/DistanceConstraintGenerator.hh>
#include <protocols/constraint_generator/DistanceConstraintGeneratorCreator.hh>

// Core headers
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/FuncFactory.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/ResidueVector.hh>
#include <core/select/residue_selector/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.constraint_generator.DistanceConstraintGenerator" );

namespace protocols {
namespace constraint_generator {

protocols::constraint_generator::ConstraintGeneratorOP
DistanceConstraintGeneratorCreator::create_constraint_generator() const
{
	return protocols::constraint_generator::ConstraintGeneratorOP( new DistanceConstraintGenerator );
}

std::string
DistanceConstraintGeneratorCreator::keyname() const
{
	return DistanceConstraintGenerator::class_name();
}

DistanceConstraintGenerator::DistanceConstraintGenerator():
	protocols::constraint_generator::ConstraintGenerator( DistanceConstraintGenerator::class_name() ),
	func_( new core::scoring::func::HarmonicFunc( 11.0, 1.0 ) ),
	selector_( new core::select::residue_selector::TrueResidueSelector )
{
}

DistanceConstraintGenerator::~DistanceConstraintGenerator()
{}

protocols::constraint_generator::ConstraintGeneratorOP
DistanceConstraintGenerator::clone() const
{
	return protocols::constraint_generator::ConstraintGeneratorOP( new DistanceConstraintGenerator( *this ) );
}

std::string
DistanceConstraintGenerator::class_name()
{
	return "DistanceConstraintGenerator";
}

void
DistanceConstraintGenerator::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data )
{
	selector_ = core::select::residue_selector::parse_residue_selector( tag, data );

	std::string const func_str = tag->getOption< std::string >( "function", "" );
	if ( !func_str.empty() ) set_function( func_str );
}

core::scoring::constraints::ConstraintCOPs
DistanceConstraintGenerator::apply( core::pose::Pose const & pose ) const
{
	if ( !func_ ) {
		std::stringstream msg;
		msg << class_name() << "::apply() Func must be set!" << std::endl;
		utility_exit_with_message( msg.str() );
	}
	if ( !selector_ ) {
		std::stringstream msg;
		msg << class_name() << "::apply() Selector must be set!" << std::endl;
		utility_exit_with_message( msg.str() );
	}
	core::scoring::constraints::ConstraintCOPs csts;

	core::select::residue_selector::ResidueVector const residues( selector_->apply( pose ) );

	for ( core::select::residue_selector::ResidueVector::const_iterator r=residues.begin(); r!=residues.end(); ++r ) {
		core::select::residue_selector::ResidueVector::const_iterator next = r;
		for ( core::select::residue_selector::ResidueVector::const_iterator r2=++next; r2!=residues.end(); ++r2 ) {
			csts.push_back( create_constraint( pose, *r, *r2 ) );
		}
	}

	return csts;
}

core::scoring::constraints::ConstraintOP
DistanceConstraintGenerator::create_constraint(
	core::pose::Pose const & pose,
	core::Size const resid1,
	core::Size const resid2 ) const
{
	core::id::AtomID const a1( pose.residue( resid1 ).nbr_atom(), resid1 );
	core::id::AtomID const a2( pose.residue( resid2 ).nbr_atom(), resid2 );
	return core::scoring::constraints::ConstraintOP( new core::scoring::constraints::AtomPairConstraint( a1, a2, func_->clone() ) );
}

core::scoring::func::FuncOP
parse_func( std::string const & func_str )
{
	std::istringstream func_data( func_str );
	std::string func_type;
	func_data >> func_type;

	core::scoring::func::FuncFactory func_factory;
	core::scoring::func::FuncOP func = func_factory.new_func( func_type );
	func->read_data( func_data );
	return func;
}

void
DistanceConstraintGenerator::set_function( std::string const & func_str )
{
	func_ = parse_func( func_str );
	TR << "Created func: " << func_str << std::endl;
}

} //protocols
} //constraint_generator

