// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file MessageListenerFactor.hh
///
/// @brief
/// @author Tim Jacobs



#ifndef INCLUDED_basic_message_listening_MessageListenerFactory_HH
#define INCLUDED_basic_message_listening_MessageListenerFactory_HH

#include <basic/message_listening/MessageListener.fwd.hh>

// C++ headers
#include <map>

#ifdef MULTI_THREADED
#ifdef CXX11

// C++11 Headers
#include <atomic>
#include <mutex>

#endif
#endif

namespace basic {
namespace message_listening {

class MessageListenerFactory {

public:

	static MessageListenerFactory* get_instance();

	MessageListenerOP get_listener(listener_tags tag);

#ifdef MULTI_THREADED
#ifdef CXX11
public:

	/// @brief This public method is meant to be used only by the
	/// utility::thread::safely_create_singleton function and not meant
	/// for any other purpose.  Do not use.
	static std::mutex & singleton_mutex();

private:
	static std::mutex singleton_mutex_;
#endif
#endif

private:
	MessageListenerFactory();
	MessageListenerFactory(MessageListenerFactory const &);
	MessageListenerFactory const & operator = (MessageListenerFactory const &);
	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static MessageListenerFactory * create_singleton_instance();

private:
#if defined MULTI_THREADED && defined CXX11
	static std::atomic< MessageListenerFactory * > instance_;
#else
	static MessageListenerFactory * instance_;
#endif

	std::map<listener_tags, MessageListenerOP> listeners_;

};

} //namespace
} //namespace
#endif
