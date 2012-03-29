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
/// @author tim

#include <core/types.hh>

#include <protocols/jd2/message_listening/MessageListener.fwd.hh>
#include <protocols/jd2/message_listening/DbMoverMessageListener.hh>
#include <protocols/jd2/message_listening/MessageListenerFactory.hh>

#include <utility/exit.hh>

#include <basic/Tracer.hh>

#include <map>

namespace protocols{
namespace jd2{
namespace message_listening{

static basic::Tracer TR("protocols.jd2.message_listening.MessageListenerFactory");

MessageListenerFactory* MessageListenerFactory::instance_(0);

MessageListenerFactory* MessageListenerFactory::get_instance()
{
	if(instance_ == 0)
	{
		instance_ = new MessageListenerFactory();
	}
	return instance_;
}

MessageListenerFactory::MessageListenerFactory()
{
	listeners_.clear();
}

MessageListenerOP MessageListenerFactory::get_listener(listener_tags tag){

	//if we already made this listener then return it, otherwise create a new one
	if(listeners_.count( tag )){
		TR << "Found existing listener for tag, returning it" << std::endl;
		return listeners_[tag];
	}

	MessageListenerOP listener;
	switch ( tag ) {
		case DB_TAG:
			TR << "Creating a new DbMoverMessageListener" << std::endl;
			listener = new DbMoverMessageListener();
			break;

		default:
			utility_exit_with_message("ERROR: you specified an invalid message listener");
			break;
	}
	listeners_[tag]=listener;
	return listener;
}

} //namespace message_listening
} //namespace jd2
} //namespace protocols

