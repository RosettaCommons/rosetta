// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/fldsgn/SheetRemodelConstraintGenerator.cc
/// @brief Remodel constraint generator for adding sheet constraints
/// @author Tom Linsky (tlinsky@gmail.com)

// Unit header
#include <protocols/fldsgn/SheetRemodelConstraintGenerator.hh>
#include <protocols/fldsgn/SheetRemodelConstraintGeneratorCreator.hh>

// Protocol headers
#include <protocols/fldsgn/SheetConstraintGenerator.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.fldsgn.SheetRemodelConstraintGenerator" );

namespace protocols {
namespace fldsgn {

SheetRemodelConstraintGenerator::SheetRemodelConstraintGenerator():
	protocols::forge::remodel::RemodelConstraintGenerator(),
	cg_( new SheetConstraintGenerator )
{
}

SheetRemodelConstraintGenerator::~SheetRemodelConstraintGenerator()
{
}

std::string
SheetRemodelConstraintGenerator::get_name() const
{
	return SheetRemodelConstraintGenerator::class_name();
}

protocols::moves::MoverOP
SheetRemodelConstraintGenerator::clone() const
{
	return protocols::moves::MoverOP( new SheetRemodelConstraintGenerator( *this ) );
}

void
SheetRemodelConstraintGenerator::parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & )
{
	cg_->parse_my_tag( tag, data );
}

void
SheetRemodelConstraintGenerator::generate_remodel_constraints( core::pose::Pose const & pose )
{
	if ( !cg_ ) {
		std::stringstream msg;
		msg << class_name() << "::generate_remodel_constraints(): "
			<< "constraint generator not set!" << std::endl;
		utility_exit_with_message( msg.str() );
	}
	add_constraints( cg_->apply( pose ) );
}

protocols::moves::MoverOP
SheetRemodelConstraintGeneratorCreator::create_mover() const
{
	return protocols::moves::MoverOP( new SheetRemodelConstraintGenerator );
}

std::string
SheetRemodelConstraintGeneratorCreator::keyname() const
{
	return SheetRemodelConstraintGenerator::class_name();
}

} //protocols
} //fldsgn

