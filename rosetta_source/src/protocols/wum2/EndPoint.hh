// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/wum2/EndPoint.fwd.hh
/// @brief  Handles communication between different roles in mpi
// Every role in mpi needs one of these for each different role it communicates with
// for example, slave only communicates with master so slave needs one
// master communicates with slave and pool, so master needs two

// role_available_mem is a functor to the role's available memory function
//
// this allows a role with multiple endpoints have its endpoints be aware of 
// other endpoint memory usage and act accordingly
/// @author Ken Jung

#ifndef INCLUDED_protocols_wum2_EndPoint_hh
#define INCLUDED_protocols_wum2_EndPoint_hh

#ifdef USEBOOSTMPI
// this is useless without mpi
#include <boost/mpi.hpp>
#include <boost/function.hpp>
#include <boost/tuple/tuple.hpp>

#include <protocols/wum2/EndPoint.fwd.hh>
#include <protocols/wum2/WorkUnit.fwd.hh>
#include <protocols/wum2/WUQueue.hh>
#include <protocols/wum2/WUQueueBuffer.hh>
#include <boost/cstdint.hpp>

#include <set>

namespace protocols {
namespace wum2 {

using namespace boost;

struct StatusResponse {
  int rank; // rank of node responding to statusrequest
  uint64_t incoming_allocated; // how much memory is allocated for incoming wu
  uint64_t outq_current_mem; // how much memory outgoing queue is using
#ifdef USEBOOSTSERIALIZE
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version) {
				ar & rank;
				ar & incoming_allocated;
				ar & outq_current_mem;
		}
#endif
};

struct StatusRequest {
  int rank; // rank of node sending statusrequest
  uint64_t max_outgoing_wu_mem; // the maximum amount of memory the outgoing wus will have
#ifdef USEBOOSTSERIALIZE
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version) {
				ar & rank;
				ar & max_outgoing_wu_mem;
		}
#endif
};

enum {
  CLEARCOMMAND = 100,
  STATUSREQUEST,
  STATUSRESPONSE,
  WORKUNITVEC,
};

class EndPoint {

  public:
    EndPoint( mpi::communicator world, uint64_t reserved_mem, uint64_t reserved_mem_multiplier, function< uint64_t () > role_available_mem );
    ~EndPoint(){}

    int mpi_rank() { return world_.rank(); }

    // real memory usage
    uint64_t current_mem() {
      return inq_.current_mem() +
        inbuf_.current_mem() +
        outq_.current_mem() +
        outbuf_.current_mem();
    }

    // need to buffer the reported current mem so that irecv doesn't preallocate too much mem
    uint64_t buffered_current_mem() {
      return  current_mem() + reserved_mem_ * reserved_mem_multiplier_;
    }

    // checks and responds to a status request, then calls functor
    // int refers to rank of node requesting the StatusResponse
    void check_and_act_status_request( function< void ( StatusResponse & , int ) > functor );

    // default functor for check_and_act_status_request
    // opens irecv from asking node in anticipation of wu isend from asking node
    // -> master sending new wu to slave
    // opens isend to asking node in anticipation of wu irecv from asking node
    // -> slave sending completed wu back to master
    void listen_wu_sendrecv( StatusResponse & r, int requesting_node );

    // isend to that node requesting a StatusResponse
    // fills up inbound_statusresponse_, acts_on_status_response required to clear it
    void send_status_request( int rank );

    // for each inbound StatusResponse, calls functor on it and then deletes it if functor return true
    void act_on_status_response( function<bool ( StatusResponse & r )> functor );

    // default functor for acts_on_status_response
    // opens irecv from asking node in anticipation of wu isend from asking node
    // -> master sending new wu to slave
    // opens isend to asking node in anticipation of wu irecv from asking node
    // -> slave sending completed wu back to master
    bool initiate_wu_sendrecv( StatusResponse & r );

    // deletes reqs and corresponding buffers that have been completed succesfully
    // also moves stuff from inbuf to inq
    void cleanup_reqs();
    
    // setup an irecv from rank with WUs from rank's outbound totaling up to mem_size
    void receive_wus( int rank, uint64_t mem_size );

    // setup an isend to rank with WUs from outbound totaling up to mem_size
    void send_wus( int rank, uint64_t mem_size );

    // checks and act to a clear queue command
    void check_and_act_clearcommand();



		bool has_open_status( int rank ) {
			return open_status_.count( rank );
		}


    // pops front of inq
    WorkUnitSP inq_popfront() {
      return inq_.pop_front();
    }
    void outq_pushback( WorkUnitSP wu) {
      outq_.push_back( wu );
    }

		uint64_t max_outgoing_wu_mem() {
			// for now, return whole queue size
			return outq_.current_mem();
		}


  private:

    mpi::communicator world_;

    // this is the amount of memory a single job may generate
    uint64_t reserved_mem_;

    // how many WUs the slave could possibly complete before master finishes sending more WU
    int reserved_mem_multiplier_;

    WUQueue inq_;
    WUQueueBuffer inbuf_;
    WUQueue outq_;
    WUQueueBuffer outbuf_;

    // these hold the buffers for statusresponse sending and receiving
    std::list< tuple< mpi::request, StatusResponse > > outbound_statusresponse_;
    std::list< tuple< mpi::request, StatusResponse > > inbound_statusresponse_;

    // these hold the buffers for statusrequest sending and receiving
    std::list< tuple< mpi::request, StatusRequest > > outbound_statusrequest_;

		std::set<int> open_status_; // holds ranks of nodes who have not responded to a statusrequest

    // keep a irecv open for status response and clearcommand always
    tuple< mpi::request, int > clearcommand_channel_;
    tuple< mpi::request, StatusRequest > statusrequest_channel_;

    function< uint64_t () > role_available_mem_;


};


} // wum2
} // protocols

#endif
#endif
