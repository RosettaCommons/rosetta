// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constraint_generator/DistanceConstraintGenerator.cc
/// @brief Generates AtomPair constraints to enforce a given distance between two residue subsets
/// @author Tom Linsky (tlinsky@gmail.com)

// Unit headers
#include <protocols/constraint_generator/DistanceConstraintGenerator.hh>
#include <protocols/constraint_generator/DistanceConstraintGeneratorCreator.hh>

// Core headers
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/FuncFactory.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/ResidueVector.hh>
#include <core/select/residue_selector/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// Boost headers
#include <boost/assign.hpp>

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
	selector1_(),
	selector2_(),
	atom_name1_( "" ),
	atom_name2_( "" )
{
}

DistanceConstraintGenerator::~DistanceConstraintGenerator() = default;

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
	std::string const selector1_name = tag->getOption< std::string >( "residue_selector1", "" );
	if ( !selector1_name.empty() ) selector1_ = core::select::residue_selector::get_residue_selector( selector1_name, data );

	std::string const selector2_name = tag->getOption< std::string >( "residue_selector2", "" );
	if ( !selector2_name.empty() ) selector2_ = core::select::residue_selector::get_residue_selector( selector2_name, data );

	atom_name1_ = tag->getOption< std::string >( "atom_name1", atom_name1_ );
	atom_name2_ = tag->getOption< std::string >( "atom_name2", atom_name2_ );

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
	if ( !selector1_ ) {
		std::stringstream msg;
		msg << class_name() << "::apply() Selector 1 must be set!" << std::endl;
		utility_exit_with_message( msg.str() );
	}
	if ( !selector2_ ) {
		std::stringstream msg;
		msg << class_name() << "::apply() Selector 2 must be set!" << std::endl;
		utility_exit_with_message( msg.str() );
	}

	using core::select::residue_selector::ResidueVector;
	using core::scoring::constraints::ConstraintOP;

	core::scoring::constraints::ConstraintCOPs csts;
	ResidueVector const subset1 = selector1_->apply( pose );
	ResidueVector const subset2 = selector2_->apply( pose );

	for ( unsigned long r : subset1 ) {
		for ( unsigned long r2 : subset2 ) {
			// don't create dist cst to same residue
			if ( r == r2 ) continue;
			csts.push_back( create_constraint( pose, r, r2 ) );
		}
	}

	TR << "Created " << csts.size() << " constraints." << std::endl;

	if ( csts.size() <= 1 ) return csts;

	TR.Debug << "Rolling constraints into ambiguous constraint" << std::endl;
	ConstraintOP amb_cst( new core::scoring::constraints::AmbiguousConstraint( csts ) );
	return boost::assign::list_of (amb_cst);
}

core::scoring::constraints::ConstraintOP
DistanceConstraintGenerator::create_constraint(
	core::pose::Pose const & pose,
	core::Size const resid1,
	core::Size const resid2 ) const
{
	core::id::AtomID a1( pose.residue( resid1 ).nbr_atom(), resid1 );
	core::id::AtomID a2( pose.residue( resid2 ).nbr_atom(), resid2 );
	if ( !atom_name1_.empty() ) {
		a1.atomno() = pose.residue( resid1 ).type().atom_index( atom_name1_ );
	}
	if ( !atom_name2_.empty() ) {
		a2.atomno() = pose.residue( resid2 ).type().atom_index( atom_name2_ );
	}
	TR.Debug << "Creating distance constraint between " << a1 << " and " << a2 << std::endl;
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

