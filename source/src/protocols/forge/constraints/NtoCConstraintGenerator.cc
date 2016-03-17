// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/forge/constraints/NtoCConstraintGenerator.cc
///
/// @brief
/// @author Nobuyasu Koga( nobuyasu@uw.edu ) , October 2009
/// @modified Tom Linsky (tlinsky@uw.edu), Nov 2012

// Unit header
#include <protocols/forge/constraints/NtoCConstraintGenerator.hh>
#include <protocols/forge/constraints/NtoCConstraintGeneratorCreator.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/func/ScalarWeightedFunc.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/id/SequenceMapping.hh>
#include <protocols/toolbox/match_enzdes_util/util_functions.hh>

// Project headers
#include <core/chemical/ResidueType.hh>

// utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

// boost headers
#include <boost/assign.hpp>

static THREAD_LOCAL basic::Tracer TR( "protocols.forge.constraints.NtoCConstraintGenerator" );

namespace protocols {
namespace forge {
namespace constraints {

std::string
NtoCConstraintGeneratorCreator::keyname() const
{
	return NtoCConstraintGeneratorCreator::mover_name();
}

protocols::moves::MoverOP
NtoCConstraintGeneratorCreator::create_mover() const {
	return protocols::moves::MoverOP( new NtoCConstraintGenerator() );
}

std::string
NtoCConstraintGeneratorCreator::mover_name()
{
	return "NtoCConstraintGenerator";
}

/// @brief
NtoCConstraintGenerator::NtoCConstraintGenerator():
	RemodelConstraintGenerator(),
	dist_( 11.0 ),
	coef_( 1.0 )
{}

/// @brief
NtoCConstraintGenerator::NtoCConstraintGenerator( Real const dist, Real const coef ):
	RemodelConstraintGenerator(),
	dist_( dist ),
	coef_( coef )
{}

/// @brief
NtoCConstraintGenerator::~NtoCConstraintGenerator() {}

void
NtoCConstraintGenerator::parse_my_tag( TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose )
{
	RemodelConstraintGenerator::parse_my_tag( tag, data, filters, movers, pose );
	set_weight( tag->getOption< core::Real >( "weight", coef_ ) );
	set_distance( tag->getOption< core::Real >( "dist", dist_ ) );
}

std::string
NtoCConstraintGenerator::get_name() const
{
	return NtoCConstraintGeneratorCreator::mover_name();
}

protocols::moves::MoverOP
NtoCConstraintGenerator::fresh_instance() const
{
	return protocols::moves::MoverOP( new NtoCConstraintGenerator() );
}

protocols::moves::MoverOP
NtoCConstraintGenerator::clone() const
{
	return protocols::moves::MoverOP( new NtoCConstraintGenerator( *this ) );
}

/// @brief set weight
void
NtoCConstraintGenerator::set_weight( Real const coef )
{
	coef_ = coef;
}

/// @brief set distance of constraint
void
NtoCConstraintGenerator::set_distance( Real const dist )
{
	dist_ = dist;
}


/// @brief
core::scoring::constraints::ConstraintCOPs
NtoCConstraintGenerator::generate_constraints( Pose const & pose )
{
	using namespace core::scoring::constraints;

	std::string tag( "constraint_between_N_&_C_terminal_Calpha" );
	Real lb( 0.0 );
	Real ub( dist_ );
	Real sd( 1.0 );
	core::scoring::func::ScalarWeightedFuncOP cstfunc( new core::scoring::func::ScalarWeightedFunc( coef_, core::scoring::func::FuncOP( new BoundFunc( lb, ub, sd, tag ) ) ) );

	Size last_residue = protocols::toolbox::match_enzdes_util::get_last_protein_residue( pose );
	Size first_residue = protocols::toolbox::match_enzdes_util::get_first_protein_residue( pose );
	TR << "first residue in NtoC generation is:" << first_residue << " and last is:" << last_residue << " out of total=" << pose.total_residue() << std::endl;
	core::id::AtomID atom1( pose.residue_type( first_residue ).atom_index( "CA" ), first_residue );
	core::id::AtomID atom2( pose.residue_type( last_residue ).atom_index( "CA" ), last_residue );
	ConstraintOP const cst( new AtomPairConstraint( atom1, atom2, cstfunc ) );

	TR << "Constraints between N- and C- terminal: " << first_residue << "-" << last_residue << ", dist=" << dist_ << std::endl;

	return boost::assign::list_of( cst );
} //generate constraints


} //namespace constraints
} //namespace forge
} //namespace protocols
