// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file DbMoverMessageListener.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_basic_message_listening_DbMoverMessageListener_hh
#define INCLUDED_basic_message_listening_DbMoverMessageListener_hh

#include <basic/message_listening/MessageListener.hh>
#include <basic/message_listening/DbMoverMessageListener.fwd.hh>

#include <numeric/types.hh>

#include <string>
#include <map>

namespace basic {
namespace message_listening {

class DbMoverMessageListener : public MessageListener{

public:
	DbMoverMessageListener();

	///@brief max_batch_id_ should only be set once.  Check ahead of time to avoid
	/// unnecessary SELECT statements.
	bool max_batch_id_set()
	{
		return initialized_;
	}

	///@brief set that max batch id.  This is only set if batch_ids_ is empty.
	///max_batch_id_ is used to ensure that a unique batch_id is chosen when
	///writing to a database that already has batch(es) stored in it.
	void set_max_batch_id(numeric::Size max_batch_id)
	{
		if(!initialized_)
		{
			max_batch_id_ = max_batch_id;
			initialized_ = true;
		}
	}

	///@brief receive the protocol id and batch id from the slave
	virtual
	void
	receive(
		std::string const & data);

	///@brief check to see if we have a protocol id and batch id. If we
	///have them then tell them to the slave. If we don't have them then
	///tell the slave to make them
	virtual
	bool
	request(
		std::string const & identifier,
		std::string & return_data);

	void
	deserialize_data(
		std::string const & data);

private:

	numeric::Size max_batch_id_;
	bool initialized_;
	numeric::Size protocol_id_;
	std::map<std::string, numeric::Size> batch_ids_;

};

} //namespace
} //namespace
#endif
