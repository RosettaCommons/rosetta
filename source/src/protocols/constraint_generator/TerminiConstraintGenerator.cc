// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constraint_generator/TerminiConstraintGenerator.cc
/// @brief Generates distance constraints between the upper and lower termini
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/constraint_generator/TerminiConstraintGenerator.hh>
#include <protocols/constraint_generator/TerminiConstraintGeneratorCreator.hh>

// Protocol headers
#include <protocols/constraint_generator/util.hh>

// Core headers
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>

// Boost headers
#include <boost/assign.hpp>

static THREAD_LOCAL basic::Tracer TR( "protocols.constraint_generator.TerminiConstraintGenerator" );

namespace protocols {
namespace constraint_generator {

protocols::constraint_generator::ConstraintGeneratorOP
TerminiConstraintGeneratorCreator::create_constraint_generator() const
{
	return protocols::constraint_generator::ConstraintGeneratorOP( new TerminiConstraintGenerator );
}

std::string
TerminiConstraintGeneratorCreator::keyname() const
{
	return TerminiConstraintGenerator::class_name();
}

TerminiConstraintGenerator::TerminiConstraintGenerator():
	protocols::constraint_generator::ConstraintGenerator( TerminiConstraintGenerator::class_name() ),
	min_distance_( 0.0 ),
	max_distance_( 11.0 ),
	sd_( 1.0 ),
	weight_( 1.0 )
{
}

TerminiConstraintGenerator::~TerminiConstraintGenerator() = default;

ConstraintGeneratorOP
TerminiConstraintGenerator::clone() const
{
	return ConstraintGeneratorOP( new TerminiConstraintGenerator( *this ) );
}

void
TerminiConstraintGenerator::parse_tag( utility::tag::TagCOP , basic::datacache::DataMap & )
{
}

core::Size
last_protein_residue( core::pose::Pose const & pose )
{
	for ( core::Size resid=pose.total_residue(); resid>=1; --resid ) {
		if ( pose.residue( resid ).is_protein() ) return resid;
	}
	return 0;
}

core::Size
first_protein_residue( core::pose::Pose const & pose )
{
	for ( core::Size resid=1; resid<=pose.total_residue(); ++resid ) {
		if ( pose.residue( resid ).is_protein() ) return resid;
	}
	return 0;
}

core::scoring::constraints::ConstraintCOPs
TerminiConstraintGenerator::apply( core::pose::Pose const & pose ) const
{
	static std::string const tag = "constraint_between_N_&_C_terminal_Calpha";
	core::scoring::func::FuncOP basefunc( new core::scoring::constraints::BoundFunc( min_distance_, max_distance_, sd_, tag ) );
	core::scoring::func::FuncOP cstfunc = scalar_weighted( basefunc, weight_ );

	core::Size const lower = first_protein_residue( pose );
	core::Size const upper = last_protein_residue( pose );

	if ( ( lower == 0 )  || ( upper == 0 ) ) {
		TR.Error << "TerminiConstraintGenerator::apply(): no protein residues were found in the pose. Not generating constraints" << std::endl;
		return core::scoring::constraints::ConstraintCOPs();
	}
	if ( !pose.residue( lower ).has( "CA" ) ) {
		TR.Error << "TerminiConstraintGenerator::apply(): residue " << lower << " does not have a CA atom, but it is a protein residue. Not generating constraints" << std::endl;
		return core::scoring::constraints::ConstraintCOPs();
	}
	if ( !pose.residue( upper ).has( "CA" ) ) {
		TR.Error << "TerminiConstraintGenerator::apply(): residue " << lower << " does not have a CA atom, but it is a protein residue. Not generating constraints" << std::endl;
		return core::scoring::constraints::ConstraintCOPs();
	}

	core::id::AtomID const atom1( pose.residue_type( lower ).atom_index( "CA" ), first_protein_residue( pose ) );
	core::id::AtomID const atom2( pose.residue_type( upper ).atom_index( "CA" ), last_protein_residue( pose ) );

	TR << "Constraining atoms " << atom1 << " and " << atom2 << ", min_distance=" << min_distance_
		<< " max_distance=" << max_distance_ << std::endl;
	core::scoring::constraints::ConstraintOP const cst( new core::scoring::constraints::AtomPairConstraint( atom1, atom2, cstfunc ) );

	return boost::assign::list_of( cst );
}

void
TerminiConstraintGenerator::set_weight( core::Real const weight )
{
	weight_ = weight;
}

void
TerminiConstraintGenerator::set_sd( core::Real const sd )
{
	sd_ = sd;
}

void
TerminiConstraintGenerator::set_min_distance( core::Real const dist )
{
	min_distance_ = dist;
}

void
TerminiConstraintGenerator::set_max_distance( core::Real const dist )
{
	max_distance_ = dist;
}

core::Real
TerminiConstraintGenerator::weight() const
{
	return weight_;
}

core::Real
TerminiConstraintGenerator::min_distance() const
{
	return min_distance_;
}

core::Real
TerminiConstraintGenerator::max_distance() const
{
	return max_distance_;
}

} //protocols
} //constraint_generator

