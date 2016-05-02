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

// Project headers

// utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// boost headers

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
	cg_()
{
}

/// @brief
NtoCConstraintGenerator::NtoCConstraintGenerator( Real const dist, Real const coef ):
	RemodelConstraintGenerator(),
	cg_()
{
	cg_.set_weight( coef );
	cg_.set_max_distance( dist );
}

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
	cg_.set_weight( tag->getOption< core::Real >( "weight", cg_.weight() ) );
	cg_.set_max_distance( tag->getOption< core::Real >( "dist", cg_.max_distance() ) );
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

/*
/// @brief set weight
void
NtoCConstraintGenerator::set_weight( Real const coef )
{
	cg_.set_weight( coef );
}

/// @brief set distance of constraint
void
NtoCConstraintGenerator::set_distance( Real const dist )
{
	cg_.set_distance( dist );
}
*/


/// @brief
void
NtoCConstraintGenerator::generate_remodel_constraints( Pose const & pose )
{
	add_constraints( cg_.apply( pose ) );
} //generate constraints


} //namespace constraints
} //namespace forge
} //namespace protocols
