// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/wum2/EndPoint.hh
/// @brief   Non-MPI version of EndPoint
/// @details This class is required because SingleNode role needs to use an EndPoint that is not MPI dependent (ie just a wrapper for 2 queues)
/// @author  Ken Jung

#ifndef INCLUDED_protocols_wum2_EndPoint_hh
#define INCLUDED_protocols_wum2_EndPoint_hh


#include <protocols/wum2/EndPoint.fwd.hh>
#include <protocols/wum2/WorkUnit.fwd.hh>
#include <protocols/wum2/WUQueue.hh>
#include <boost/cstdint.hpp>
#include <boost/function.hpp>

namespace protocols {
namespace wum2 {

#ifndef __native_client__

struct StatusRequest;
struct StatusResponse;

class EndPoint {

public:
	EndPoint( boost::function < boost::uint64_t () > role_available_mem );
	virtual ~EndPoint(){}

	// real memory usage
	virtual boost::uint64_t current_mem() {
		return inq_.current_mem() +
			outq_.current_mem();
	}

	WUQueue & inq() { return inq_; }
	WUQueue & outq() { return outq_; }

	boost::uint64_t max_outgoing_wu_mem() {
		// for now, return whole queue size
		return outq_.current_mem();
	}

	// MPI derived class uses these
	// this is easier than doing static casts to use MPI_EndPoint fxns
	// non-void functions made pure virtual to avoid compiler warnings ~ Labonte
	virtual void check_and_act_status_request( boost::function< void ( StatusResponse & , int ) > /*functor*/ ) {}
	virtual void check_and_act_clearcommand() {}
	virtual void cleanup_reqs() {}
	virtual bool has_open_status( int /*rank*/ ) = 0;
	virtual void send_status_request( int /*rank*/ ){}
	virtual void listen_wu_sendrecv( StatusResponse & /*r*/, int /*requesting_node*/ ){}
	virtual bool initiate_wu_sendrecv( StatusResponse & /*r*/ ) = 0;
	virtual void act_on_status_response( boost::function<bool ( StatusResponse & r )> /*functor*/ ) {}

protected:

	WUQueue inq_;
	WUQueue outq_;

	boost::function< boost::uint64_t () > role_available_mem_;

};

struct StatusResponse {
	int rank; // rank of node responding to statusrequest
	boost::uint64_t incoming_allocated; // how much memory is allocated for incoming wu
	boost::uint64_t outq_current_mem; // how much memory outgoing queue is using
};

struct StatusRequest {
	int rank; // rank of node sending statusrequest
	boost::uint64_t max_outgoing_wu_mem; // the maximum amount of memory the outgoing wus will have
};

#endif

} // wum2
} // protocols

#endif
