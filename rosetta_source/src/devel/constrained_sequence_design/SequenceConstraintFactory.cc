// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   
/// @brief
/// @author Javier Castellanos (javiercv@uw.edu)

// Unit Headers
#include <devel/constrained_sequence_design/SequenceConstraintFactory.hh>

// Package Headers
#include <devel/constrained_sequence_design/SequenceConstraintCreator.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.hh>

// Utility Headers
#include <basic/Tracer.hh>

namespace devel {
namespace constrained_sequence_design {

static basic::Tracer TR( "devel.constrained_sequence_design.SeuquenceConstraintFactory" );

SequenceConstraintFactory * SequenceConstraintFactory::instance_( 0 );
SequenceConstraintFactory::SequenceConstraintFactory() { } 
SequenceConstraintFactory::~SequenceConstraintFactory() { } 

SequenceConstraintFactory *
SequenceConstraintFactory::get_instance() {
	if ( ! instance_ )
		instance_ = new SequenceConstraintFactory;
	return instance_;
} // get_instance

void 
SequenceConstraintFactory::factory_register( SequenceConstraintCreatorOP c)
{
	if(creator_map_.count(c->keyname()) > 0)
		utility_exit_with_message( "Error: SequenceConstraint with name " + c->keyname() + " has already been registered in SequenceConstraintFactory");
	else
		creator_map_[ c->keyname() ] = c;
} // factory_register

/// brief return a new SequenceConstraint object of class constraint_type
SequenceConstraintOP 
SequenceConstraintFactory::newSequenceConstraint(const std::string & constraint_type)
{
	SequenceConstraintMap::iterator iter( creator_map_.find( constraint_type ) );
	if ( iter != creator_map_.end() ) {
		if ( ! iter->second )
			utility_exit_with_message( "Error: SequenceConstraintCreatorOP prototype for " + constraint_type + " is NULL!");
		return iter->second->create();
	} else {
		TR << "Available sequence constraints: ";
		for( SequenceConstraintMap::const_iterator it = creator_map_.begin(); it != creator_map_.end(); ++ it)
			TR << it->first << ", ";
		TR << std::endl;
		utility_exit_with_message( constraint_type + " is not known to the SequenceConstraintFactory. Was it registered via a SequenceConstraintRegistrator in one of the init.cc files (devel/init.cc or protocols/init.cc)?" );
		return NULL;
	}
} // newSequenceConstraint(const std::string & constriant_type_)

//SequenceConstraintOP 
//SequenceConstraintFactory::newSequenceConstraint(std::string type,  TagPtr const tag) {
//} //newSequenceConstraint(std::string type, TagPtr const  tag)

	
SequenceConstraintOP 
SequenceConstraintFactory::newSequenceConstraint(
	TagPtr const tag,
	protocols::moves::DataMap & data,
	Pose const & pose )
{
	SequenceConstraintOP cnstr( newSequenceConstraint( tag->getName() ) );
	runtime_assert( cnstr );
	cnstr->parse_my_tag( tag, data, pose );
	return cnstr;
} // newSequenceConstraint

} // namespace devel
} // namespace constrained_sequence_design 
