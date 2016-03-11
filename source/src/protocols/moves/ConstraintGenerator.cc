// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/ConstraintGenerator.cc
///
/// @brief
/// @author Florian Richter, floric@u.washington.edu, april 2009
/// @modified Tom Linsky, tlinsky@uw.edu

#include <protocols/moves/ConstraintGenerator.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/id/SequenceMapping.hh>

#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.moves.ConstraintGenerator" );
namespace protocols {
namespace moves {

std::map<std::string,core::scoring::constraints::ConstraintCOPs> ConstraintGenerator::cst_map_;

/// @details Auto-generated virtual destructor
ConstraintGenerator::~ConstraintGenerator()
{}


ConstraintGenerator::ConstraintGenerator():
	Mover("ConstraintGenerator"),
	id_( "unnamed_constraint_generator" ),
	csts_()
{}

/// @details When called, generates constraints for the pose, and adds them to the pose
void
ConstraintGenerator::apply( core::pose::Pose & pose )
{
	init( pose );
	add_constraints_to_pose( pose );
}

/// @brief This is called if this mover is instantiated from XML
void
ConstraintGenerator::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	// if there are any options that we might want, they will go here...
	id_ = tag->getOption<std::string>( "name", id_ );
	TR << "Setting id for type " << this->get_name() << " = " << id_ << std::endl;
}

void
ConstraintGenerator::add_constraints_to_pose( core::pose::Pose & pose )
{
	csts_.clear();
	clear_stored_constraints();

	TR.Debug << "Generating remodel constraints" << std::endl;
	core::scoring::constraints::ConstraintCOPs csts = generate_constraints( pose );
	TR.Debug << "Done generating csts" << std::endl;

	//safeguard against an RCG not generating anything
	if ( csts.size() == 0 ) return;

	TR.Debug << this->get_name() << " generated " << csts.size() << " constraints." << std::endl;

	csts_ = pose.add_constraints( csts );
	store_constraints( csts_ );
}

void
ConstraintGenerator::remove_constraints_from_pose( core::pose::Pose & pose ) const
{
	core::scoring::constraints::ConstraintCOPs remodel_csts( csts_ );

	TR << this->get_name() << " is about to try to remove " << remodel_csts.size() << " constraints." << std::endl;

	if ( remodel_csts.size() == 0 ) return;

	//TR << this->get_name() << " is about to try to remove " << remodel_csts.size() << " constraints." << std::endl;

	if ( ! pose.remove_constraints( remodel_csts, true ) ) {
		throw EXCN_RemoveCstsFailed();
	}
}

core::scoring::constraints::ConstraintCOPs const &
ConstraintGenerator::constraints() const
{
	return csts_;
}

std::string
ConstraintGenerator::id() const
{
	return id_;
}

void
ConstraintGenerator::set_id( std::string const & id )
{
	id_ = id;
}

void
ConstraintGenerator::clear_stored_constraints()
{
	if ( id_ == "" ) return;
	// find id for this class in the map
	std::map< std::string, core::scoring::constraints::ConstraintCOPs >::iterator cst_it( cst_map_.find( id_ ) );
	// if the constraints are found, warn the user and erase the old data
	if ( cst_it != cst_map_.end() ) {
		TR << "Overwriting constraints for " << this->get_name() << " named " << id_ << std::endl;
		cst_map_.erase( cst_it );
	}
}

void
ConstraintGenerator::store_constraints( core::scoring::constraints::ConstraintCOPs const & csts )
{
	if ( id_ == "" ) {
		TR.Warning << "ID is not set for this " << this->get_name() << " object. Constraints will not be removable by XML." << std::endl;
		return;
	}
	// clear any stored csts
	clear_stored_constraints();

	// store the csts
	TR.Debug << "Storing constraints as " << id_ << std::endl;
	cst_map_.insert( std::pair< std::string, core::scoring::constraints::ConstraintCOPs >( id_, csts ) );
}

core::scoring::constraints::ConstraintCOPs const
ConstraintGenerator::lookup_stored_constraints( std::string const & id )
{
	if ( id == "" ) {
		utility_exit_with_message( "ID is not set! The constraint set returned will be empty. Something is being mis-used." );
	}
	// find the constraint set in map
	std::map< std::string, core::scoring::constraints::ConstraintCOPs >::const_iterator cst_it( cst_map_.find( id ) );
	if ( cst_it == cst_map_.end() ) {
		utility_exit_with_message( "Tried to remove constraints that aren't stored for ID=" + id );
	}
	return cst_it->second;
}

} //namespace moves
} //namespace protocols
