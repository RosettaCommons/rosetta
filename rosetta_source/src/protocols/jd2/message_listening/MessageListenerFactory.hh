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
/// @author tim



#ifndef INCLUDED_protocols_jd2_message_listening_MessageListenerFactory_HH
#define INCLUDED_protocols_jd2_message_listening_MessageListenerFactory_HH

#include <protocols/jd2/message_listening/MessageListener.fwd.hh>

#include <map>

namespace protocols {
namespace jd2 {
namespace message_listening {    

enum listener_tags {
	DB_TAG = 42
};

class MessageListenerFactory {

public:
    MessageListenerOP get_listener(listener_tags tag);
    
private:
    std::map<listener_tags, MessageListenerOP> listeners_;
        
};

} //namespace message_listening
} //namespace jd2
} //namespace protocols
#endif
