// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file MessageListenerFactory.cc
///
/// @brief A factory for lazily initializing message listeners. This should be used in conjunction the the MPIWorkPoolJobDistributor's message listening functionality
/// @author Tim Jacobs

#include <basic/message_listening/MessageListener.fwd.hh>
#include <basic/message_listening/DbMoverMessageListener.hh>
#include <basic/message_listening/MessageListenerFactory.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/thread/threadsafe_creation.hh>

// Basic headers
#include <basic/Tracer.hh>

// C++ headers
#include <map>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

namespace basic{
namespace message_listening{

static thread_local basic::Tracer TR( "basic.message_listening.MessageListenerFactory" );

#if defined MULTI_THREADED && defined CXX11
std::atomic< MessageListenerFactory * > MessageListenerFactory::instance_( 0 );
#else
MessageListenerFactory * MessageListenerFactory::instance_( 0 );
#endif

#ifdef MULTI_THREADED
#ifdef CXX11

std::mutex MessageListenerFactory::singleton_mutex_;

std::mutex & MessageListenerFactory::singleton_mutex() { return singleton_mutex_; }

#endif
#endif

/// @brief static function to get the instance of ( pointer to) this singleton class
MessageListenerFactory * MessageListenerFactory::get_instance()
{
	boost::function< MessageListenerFactory * () > creator = boost::bind( &MessageListenerFactory::create_singleton_instance );
	utility::thread::safely_create_singleton( creator, instance_ );
	return instance_;
}

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
	if(listeners_.count( tag )){
		TR.Debug << "Found existing listener for tag, returning it" << std::endl;
		return listeners_[tag];
	}

	MessageListenerOP listener;
	switch ( tag ) {
		case DATABASE_PROTOCOL_AND_BATCH_ID_TAG:
			TR.Debug << "Creating a new DbMoverMessageListener" << std::endl;
			listener = new DbMoverMessageListener();
			break;

		default:
			utility_exit_with_message("ERROR: you specified an invalid message listener");
			break;
	}
	listeners_[tag]=listener;
	return listener;
}

} //namespace
} //namespace

