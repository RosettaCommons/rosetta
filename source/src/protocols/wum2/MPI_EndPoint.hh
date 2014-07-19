// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/wum2/MPI_EndPoint.hh
/// @brief  Handles communication between different roles in mpi
// Every role in mpi needs one of these for each different role it communicates with
// for example, slave only communicates with master so slave needs one
// master communicates with slave and pool, so master needs two

// role_available_mem is a functor to the role's available memory function
//
// this allows a role with multiple endpoints have its endpoints be aware of
// other endpoint memory usage and act accordingly
/// @author Ken Jung

#ifndef INCLUDED_protocols_wum2_MPI_EndPoint_hh
#define INCLUDED_protocols_wum2_MPI_EndPoint_hh

#ifdef USEBOOSTMPI
#include <boost/mpi.hpp>
#include <boost/function.hpp>
#include <boost/tuple/tuple.hpp>

#include <protocols/wum2/MPI_EndPoint.fwd.hh>
#include <protocols/wum2/EndPoint.hh>
#include <protocols/wum2/WUQueueBuffer.hh>

#include <set>

namespace protocols {
namespace wum2 {

enum {
  CLEARCOMMAND = 100,
  STATUSREQUEST,
  STATUSRESPONSE,
  WORKUNITVEC,
};

class MPI_EndPoint : public EndPoint {

public:
	MPI_EndPoint( boost::mpi::communicator world, boost::function< uint64_t () > role_available_mem );
	~MPI_EndPoint(){}

	int mpi_rank() { return world_.rank(); }

    // real memory usage
	virtual boost::uint64_t current_mem() {
		return inq_.current_mem() +
				inbuf_.current_mem() +
				outq_.current_mem() +
				outbuf_.current_mem();
    }

	// checks and responds to a status request, then calls functor
	// int refers to rank of node requesting the StatusResponse
	virtual void check_and_act_status_request( boost::function< void ( StatusResponse & , int ) > functor );

	// default functor for check_and_act_status_request
	// opens irecv from asking node in anticipation of wu isend from asking node
	// -> master sending new wu to slave
	// opens isend to asking node in anticipation of wu irecv from asking node
	// -> slave sending completed wu back to master
	virtual void listen_wu_sendrecv( StatusResponse & r, int requesting_node );

	// isend to that node requesting a StatusResponse
	// fills up inbound_statusresponse_, acts_on_status_response required to clear it
	virtual void send_status_request( int rank );

	// for each inbound StatusResponse, calls functor on it and then deletes it if functor return true
	virtual void act_on_status_response( boost::function<bool ( StatusResponse & r )> functor );

	// default functor for acts_on_status_response
	// opens irecv from asking node in anticipation of wu isend from asking node
	// -> master sending new wu to slave
	// opens isend to asking node in anticipation of wu irecv from asking node
	// -> slave sending completed wu back to master
	virtual bool initiate_wu_sendrecv( StatusResponse & r );

	// deletes reqs and corresponding buffers that have been completed succesfully
	// also moves stuff from inbuf to inq
	virtual void cleanup_reqs();

	// setup an irecv from rank with WUs from rank's outbound totaling up to mem_size
	void receive_wus( int rank, boost::uint64_t mem_size );

	// setup an isend to rank with WUs from outbound totaling up to mem_size
	void send_wus( int rank, boost::uint64_t mem_size );

	// checks and act to a clear queue command
	virtual void check_and_act_clearcommand();

	virtual bool has_open_status( int rank ) {
		return open_status_.count( rank );
	}

private:

	boost::mpi::communicator world_;

	WUQueueBuffer inbuf_;
	WUQueueBuffer outbuf_;

	// these hold the buffers for statusresponse sending and receiving
	std::list< boost::tuple< boost::mpi::request, StatusResponse > > outbound_statusresponse_;
	std::list< boost::tuple< boost::mpi::request, StatusResponse > > inbound_statusresponse_;

	// these hold the buffers for statusrequest sending and receiving
	std::list< boost::tuple< boost::mpi::request, StatusRequest > > outbound_statusrequest_;

	std::set<int> open_status_; // holds ranks of nodes who have not responded to a statusrequest

	// keep a irecv open for status response and clearcommand always
	boost::tuple< boost::mpi::request, int > clearcommand_channel_;
	boost::tuple< boost::mpi::request, StatusRequest > statusrequest_channel_;

};


} // wum2
} // protocols

#endif
#endif
