// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/constraint_generator/ConstraintGenerator.cc
///
/// @brief Abstract class for generating constraints
/// @author Tom Linsky ( tlinsky at uw dot edu )

// Unit Headers
#include <protocols/constraint_generator/ConstraintGenerator.hh>

// Protocol Headers
#include <protocols/constraint_generator/ConstraintsManager.hh>

// Core headers
#include <core/id/SequenceMapping.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/func/ScalarWeightedFunc.hh>

// Basic/Utility Headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.moves.ConstraintGenerator" );

namespace protocols {
namespace constraint_generator {

ConstraintGenerator::ConstraintGenerator( std::string const & class_name ):
	utility::pointer::ReferenceCount(),
	class_name_( class_name ),
	id_( "unnamed_constraint_generator" )
{}

ConstraintGenerator::~ConstraintGenerator()
{}

/// @brief This is called if this mover is instantiated from XML
void
ConstraintGenerator::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data )
{
	// if there are any options that we might want, they will go here...
	id_ = tag->getOption<std::string>( "name" );
	parse_tag( tag, data );
}

std::string const &
ConstraintGenerator::id() const
{
	return id_;
}

void
ConstraintGenerator::set_id( std::string const & id )
{
	id_ = id;
}

std::string const &
ConstraintGenerator::class_name() const
{
	return class_name_;
}

void
ConstraintGenerator::clear_stored_constraints() const
{
	if ( id_ == "" ) return;

	if ( ConstraintsManager::get_instance()->constraints_exist( id_ ) ) {
		ConstraintsManager::get_instance()->remove_constraints( id_ );
	}
}

void
ConstraintGenerator::store_constraints( core::scoring::constraints::ConstraintCOPs const & csts ) const
{
	if ( id_ == "" ) {
		TR.Warning << "ID is not set for this " << this->class_name() << " object. Constraints will not be removable by XML." << std::endl;
		return;
	}
	ConstraintsManager::get_instance()->store_constraints( id_, csts );

	// store the csts
	TR.Debug << "Stored constraints as " << id_ << std::endl;
}

core::scoring::constraints::ConstraintCOPs const &
ConstraintGenerator::retrieve_constraints() const
{
	return ConstraintsManager::get_instance()->retrieve_constraints( id_ );
}

} //namespace constraint_generator
} //namespace protocols
