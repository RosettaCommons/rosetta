// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/wum2/MPI_EndPoint.cc
/// @brief  Handles communication between different roles in mpi
// Every role in mpi needs one of these for each different role it communicates with
// for example, slave only communicates with master so slave needs one
// master communicates with slave and pool, so master needs two

// role_available_mem is a functor to the role's available memory function
//
// this allows a role with multiple endpoints have its endpoints be aware of
// other endpoint memory usage and act accordingly
/// @author Ken Jung

#ifdef USEBOOSTMPI
// this is useless without mpi
#include <protocols/wum2/MPI_EndPoint.hh>
#include <protocols/wum2/WorkUnit.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/els.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

namespace protocols {
namespace wum2 {

using namespace boost;

static thread_local basic::Tracer TR( "protocols.wum2.MPI_EndPoint" );

MPI_EndPoint::MPI_EndPoint( mpi::communicator world, function< boost::uint64_t () > role_available_mem  )
  : world_( world ),
  EndPoint( role_available_mem ) {
    // set up initial channels
    clearcommand_channel_.get<0>() = world_.irecv( mpi::any_source, CLEARCOMMAND, clearcommand_channel_.get<1>() );
    statusrequest_channel_.get<0>() = world_.irecv( mpi::any_source, STATUSREQUEST, statusrequest_channel_.get<1>() );
}

void MPI_EndPoint::check_and_act_status_request( function< void ( StatusResponse & , int ) > functor ) {
  if ( statusrequest_channel_.get<0>().test().is_initialized() ) {
    mpi::request req;
    StatusResponse m;
    m.rank = mpi_rank();
		boost::uint64_t free_mem = role_available_mem_();
		boost::uint64_t available_work = statusrequest_channel_.get<1>().max_outgoing_wu_mem;
    m.incoming_allocated = free_mem < available_work ? free_mem : available_work;
    m.outq_current_mem = outq_.current_mem();

    tuple< mpi::request, StatusResponse > tmp = make_tuple ( req , m );
    outbound_statusresponse_.push_back( tmp );
		std::list< tuple< mpi::request, StatusResponse > >::reverse_iterator itr = outbound_statusresponse_.rbegin();
    itr->get<0>() = world_.isend( statusrequest_channel_.get<1>().rank , STATUSRESPONSE, itr->get<1>() );

    // call functor
    functor( itr->get<1>(), statusrequest_channel_.get<1>().rank );
		// refresh statusrequest channel
		statusrequest_channel_.get<0>() = world_.irecv( mpi::any_source, STATUSREQUEST, statusrequest_channel_.get<1>() );
  }
}

void MPI_EndPoint::listen_wu_sendrecv( StatusResponse & r, int requesting_node ) {
  // assumption here is that the asking node will do an isend/irecv to clear out this node's queues
  // opens irecv from asking node in anticipation of wu isend from asking node
  // -> master sending new wu to slave
  // opens isend to asking node in anticipation of wu irecv from asking node
  // -> slave sending completed wu back to master
  receive_wus( requesting_node, r.incoming_allocated );
  send_wus( requesting_node, r.outq_current_mem );
}

void MPI_EndPoint::send_status_request( int rank ) {
	mpi::request req;
	StatusRequest m;
	m.rank = mpi_rank();
	m.max_outgoing_wu_mem = max_outgoing_wu_mem();
  tuple< mpi::request, StatusRequest > tmpa = make_tuple(req, m);
  outbound_statusrequest_.push_back( tmpa );
	std::list< tuple< mpi::request, StatusRequest > >::reverse_iterator jtr = outbound_statusrequest_.rbegin();
  jtr->get<0>() = world_.isend( rank, STATUSREQUEST, jtr->get<1>() );

	mpi::request reqn;
	StatusResponse n;
  tuple< mpi::request, StatusResponse > tmp = make_tuple(reqn, n);
  inbound_statusresponse_.push_back( tmp );
	std::list< tuple< mpi::request, StatusResponse > >::reverse_iterator itr = inbound_statusresponse_.rbegin();
  itr->get<0>() = world_.irecv( rank, STATUSRESPONSE, itr->get<1>() );

	open_status_.insert(rank);
}


void MPI_EndPoint::act_on_status_response( function<bool ( StatusResponse & r )> f ){
  std::list< tuple< mpi::request, StatusResponse > >::iterator itr = inbound_statusresponse_.begin();
	while( itr != inbound_statusresponse_.end() ) {
    if ( itr->get<0>().test().is_initialized() ) {
      if( f( itr->get<1>() ) ) {
        // only replace if status response send/recv can actually be fufilled
        // may not be fufillfed due to memory limits and current memory use
				open_status_.erase( itr->get<1>().rank ) ;
        itr = inbound_statusresponse_.erase( itr );
			} else {
				++itr;
			}
    } else {
			++itr;
		}
  }
}

bool MPI_EndPoint::initiate_wu_sendrecv( StatusResponse & r ) {
  // assumption here is that only reason this node initiated a
  // StatusRequest->StatusResponse chain was to send and recv WU to the responding node
  //
  // we can't delete the status response until both send and recv are fulfilled
  bool can_send = false;
  bool can_recv = false;

  if( r.outq_current_mem == 0 ) {
    can_recv = true;
  } else if (r.outq_current_mem > role_available_mem_() ) {
    // once we get in here, may never actually be able to get out of this
    can_recv = false;
  } else {
    can_recv = true;
    receive_wus( r.rank, r.outq_current_mem );
  }

  if( r.incoming_allocated == 0 ) {
    can_send = true;
  /*  just sending empty message seems less complicated, more overhead though
   *
   } else if (outq.size() == 0 ) {
    can_send = false;
    */
  } else {
    can_send = true;
    send_wus( r.rank, r.incoming_allocated );
  }
  return can_send && can_recv;
}

void MPI_EndPoint::cleanup_reqs() {
  // buffer reqs
  outbuf_.cleanup_reqs();
  std::vector< WorkUnitSP > tmp = inbuf_.cleanup_reqs();
  for( std::vector< WorkUnitSP >::iterator itr = tmp.begin(), end = tmp.end(); itr != end; ++itr ) {
    if( (*itr)->prioritize() ) {
      inq_.push_front( *itr );
    } else {
      inq_.push_back( *itr );
    }
  }
  // outbound statusresponse reqs
 	std::list< tuple< mpi::request, StatusResponse > >::iterator itr = outbound_statusresponse_.begin();
	while( itr != outbound_statusresponse_.end() ){
    if ( itr->get<0>().test().is_initialized() ) {
      // outbound statusresponse was sent successfully
      itr = outbound_statusresponse_.erase( itr );
		} else {
			++itr;
		}
  }
	// outbound statusrequest reqs
  std::list< tuple< mpi::request, StatusRequest > >::iterator jtr = outbound_statusrequest_.begin();
	while( jtr != outbound_statusrequest_.end() ) {
    if ( jtr->get<0>().test().is_initialized() ) {
      // outbound statusresponse was sent successfully
      jtr = outbound_statusrequest_.erase( jtr );
		} else {
			++jtr;
		}
	}
}

void MPI_EndPoint::receive_wus( int rank, boost::uint64_t mem_size ) {
  if( mem_size != 0 ) {
    WUQueueBuffer::riterator itr = inbuf_.allocate_buffer( mem_size );
    itr->get<1>() = world_.irecv( rank, WORKUNITVEC, *(itr->get<2>()) );
  }
}

void MPI_EndPoint::send_wus( int rank, boost::uint64_t mem_size ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
  // this is kind of inefficient, copying the vector twice, but easiest way to track mem use
  if( mem_size != 0 ) {
    std::vector< WorkUnitSP > tmp;
    boost::uint64_t current_size = 0;
		int counter = 0;
		// slave should not have any outbound limit
		// master should have low outbound limit to promote distributing work among slaves
		// otherwise will send all work to one slave
		// use the fact that master rank is always lower than slave rank to figure out who we are
		// only implement limit for master send for now
		// putting this in here is bad, should be one level up in master/slave/baserole

		int num_masters =
		(option[OptionKeys::els::num_traj]() / option[OptionKeys::els::traj_per_master]() ) +
		(!( option[OptionKeys::els::num_traj]() % option[OptionKeys::els::traj_per_master]() == 0 )) ;
		int num_slaves = (world_.size() - num_masters)/num_masters;
    while( outq_.size_front() && current_size + outq_.size_front() <= mem_size ) {
			current_size += outq_.size_front();
      tmp.push_back( outq_.pop_front() );

			if( rank > mpi_rank() )
				counter++;
			if( outq_.size() < num_slaves || counter >= 2 ) break;

    }
    WUQueueBuffer::riterator itr = outbuf_.allocate_buffer( current_size );
    itr->get<2>()->insert( itr->get<2>()->end(), tmp.begin(), tmp.end() );
    itr->get<1>() = world_.isend( rank, WORKUNITVEC, *(itr->get<2>()) );
  }
}

void MPI_EndPoint::check_and_act_clearcommand() {
  // 2 potential problems here:
  // 1) can't clear inbuf because of possible dangling isend -> zombie WUs
  // mark and sweep?
  // 2) "security issue" anyone can send a clear command to a node -> who cares about security
  if ( clearcommand_channel_.get<0>().test().is_initialized() ) {
    inq_.clear();
		// refresh clearcommand channel
		clearcommand_channel_.get<0>() = world_.irecv( mpi::any_source, CLEARCOMMAND, clearcommand_channel_.get<1>() );
  }
}

} // wum2
} // protocols
#endif
