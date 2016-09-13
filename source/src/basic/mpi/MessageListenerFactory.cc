// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file MessageListenerFactory.cc
///
/// @brief A factory for lazily initializing message listeners. This should be used in conjunction the the MPIWorkPoolJobDistributor's message listening functionality
/// @author Tim Jacobs

#include <basic/mpi/MessageListener.fwd.hh>
#include <basic/mpi/DbMoverMessageListener.hh>
#include <basic/mpi/MessageListenerFactory.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/thread/threadsafe_creation.hh>

// Basic headers
#include <basic/Tracer.hh>

// C++ headers
#include <map>

namespace basic {
namespace mpi {

static THREAD_LOCAL basic::Tracer TR( "basic.mpi.MessageListenerFactory" );

MessageListenerFactory *
MessageListenerFactory::create_singleton_instance()
{
	return new MessageListenerFactory;
}

MessageListenerFactory::MessageListenerFactory()
{
	listeners_.clear();
}

MessageListenerOP
MessageListenerFactory::get_listener(
	listener_tags tag
){

	//if we already made this listener then return it, otherwise create a new one
	if ( listeners_.count( tag ) ) {
		TR.Debug << "Found existing listener for tag, returning it" << std::endl;
		return listeners_[tag];
	}

	MessageListenerOP listener;
	switch ( tag ) {
	case DATABASE_PROTOCOL_AND_BATCH_ID_TAG :
		TR.Debug << "Creating a new DbMoverMessageListener" << std::endl;
		listener = MessageListenerOP( new DbMoverMessageListener() );
		break;

	default :
		utility_exit_with_message("ERROR: you specified an invalid message listener");
		break;
	}
	listeners_[tag]=listener;
	return listener;
}

} //namespace
} //namespace

