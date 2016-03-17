// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/forge/constraints/RemoveConstraints.cc
///
/// @brief
/// @author Tom Linsky (tlinsky@uw.edu), Nov 2012

// Unit header
#include <protocols/forge/constraints/RemoveConstraints.hh>
#include <protocols/forge/constraints/RemoveConstraintsCreator.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <protocols/moves/ConstraintGenerator.hh>

// Project headers
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <protocols/rosetta_scripts/util.hh>

// utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer TR( "protocols.forge.constraints.RemoveConstraints" );

namespace protocols {
namespace forge {
namespace constraints {

std::string
RemoveConstraintsCreator::keyname() const
{
	return RemoveConstraintsCreator::mover_name();
}

protocols::moves::MoverOP
RemoveConstraintsCreator::create_mover() const {
	return protocols::moves::MoverOP( new RemoveConstraints() );
}

std::string
RemoveConstraintsCreator::mover_name()
{
	return "RemoveConstraints";
}

/// @brief
RemoveConstraints::RemoveConstraints()
: Mover(),
	generator_( /* NULL */ ),
	generator_id_( "" )
{}

/// @brief
RemoveConstraints::RemoveConstraints( protocols::moves::ConstraintGeneratorOP generator )
: Mover()
{
	set_generator( generator );
}

/// @brief
RemoveConstraints::~RemoveConstraints() {}

void
RemoveConstraints::parse_my_tag( TagCOP const tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & )
{
	generator_id_ = tag->getOption< std::string >( "generator", generator_id_ );
	if ( generator_id_ == "" ) {
		utility_exit_with_message( "No Cst Generator was specified in the XML for RemoveConstraints." );
	}
	protocols::moves::MoverOP mover = protocols::rosetta_scripts::parse_mover( generator_id_, movers );
	debug_assert( utility::pointer::dynamic_pointer_cast< protocols::moves::ConstraintGenerator >( mover ) );
	protocols::moves::ConstraintGeneratorOP cg;
	if ( ( cg = utility::pointer::static_pointer_cast< protocols::moves::ConstraintGenerator >( mover ) ) ) {
		set_generator( cg );
	} else {
		utility_exit_with_message( "Error parsing generator option to RemoveConstraints: the specified mover " + generator_id_ + " is not a constraint generator." );
	}
	TR << "Cst generator =" << generator_->id() << " with name=" << generator_id_ << std::endl;
}

std::string
RemoveConstraints::get_name() const
{
	return RemoveConstraintsCreator::mover_name();
}

protocols::moves::MoverOP
RemoveConstraints::fresh_instance() const
{
	return protocols::moves::MoverOP( new RemoveConstraints() );
}

protocols::moves::MoverOP
RemoveConstraints::clone() const
{
	return protocols::moves::MoverOP( new RemoveConstraints( *this ) );
}

void
RemoveConstraints::set_generator( protocols::moves::ConstraintGeneratorOP generator )
{
	runtime_assert( generator != 0 );
	generator_ = generator;
}

/// @brief find the constraint set added by the generator and remove them
void
RemoveConstraints::apply( core::pose::Pose & pose )
{
	runtime_assert( generator_ || ( generator_id_ != "" ) );
	// if the generator_id_ is set, look up constraints in the static map and remove them.
	core::Size const start_ncsts = pose.constraint_set()->get_all_constraints().size();
	if ( generator_id_ != "" ) {
		TR << "Before removing csts from " << generator_->get_name() << ", there were " << start_ncsts << " constraints in the pose." << std::endl;
		pose.remove_constraints( protocols::moves::ConstraintGenerator::lookup_stored_constraints( generator_id_ ) );
		TR << "There are " << pose.constraint_set()->get_all_constraints().size() << " constraints remaining in the pose." << std::endl;
	} else {
		// otherwise, use the pointer
		TR << "Before removing csts from " << generator_->get_name() << ", there were " << start_ncsts << " constraints in the pose." << std::endl;
		generator_->remove_constraints_from_pose( pose );
		TR << "There are " << pose.constraint_set()->get_all_constraints().size() << " constraints remaining in the pose." << std::endl;
	}
}


} //namespace constraints
} //namespace forge
} //namespace protocols
