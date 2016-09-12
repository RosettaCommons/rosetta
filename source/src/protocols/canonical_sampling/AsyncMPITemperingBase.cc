// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/canonical_sampling/AsyncMPITemperingBaseMover.cc
/// @brief AsyncMPITemperingBase methods implemented
/// @author Oliver Lange ( oliver.lange@tum.de )


// Unit Headers
#include <protocols/canonical_sampling/AsyncMPITemperingBase.hh>

// protocols headers
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/canonical_sampling/ThermodynamicObserver.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>

#include <protocols/rosetta_scripts/util.hh>

//#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/Job.hh>

// core headers
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>

// numeric headers
#include <numeric/random/random.hh>

// utility headers
#include <utility/file/file_sys_util.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tag/Tag.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

// C++ Headers
#include <cmath>

using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer tr( "protocols.canonical_sampling.AsyncMPITemperingBase" );

namespace protocols {
namespace canonical_sampling {
using namespace core;

//int const mpi_LEVEL_INFORM = 1;
int const mpi_size_XCHANGE_READY = 2;
int const mpi_size_SWAP_INVITE = 1;
int const mpi_size_FINISHED = 1;
#ifdef USEMPI
int const mpi_magic_number = 100;
int const mpi_XCHANGE_READY = 1 + mpi_magic_number;
int const mpi_FINISHED = 3 + mpi_magic_number;
int const mpi_FINISHED_ACK = 4 + mpi_magic_number;
int const mpi_size_FINISHED_ACK = 1;
int const mpi_SWAP_INVITE = 2 + mpi_magic_number;
#endif

AsyncMPITemperingBase::AsyncMPITemperingBase() :
	rank_( -1 ),
	comm_size_( 0 ),
	ready_for_exchange_( false ),
	finished_( false ),
	recv_level_buffer_( nullptr ),
	send_level_buffer_( nullptr ),
	swap_invite_send_buffer_( nullptr ),
	swap_invite_recv_buffer_( nullptr ),
	finished_recv_buffer_( nullptr ),
	finished_send_buffer_( nullptr )
{
#ifdef USEMPI
	recv_level_requests_ = NULL;
	swap_invite_send_requests_ = NULL;
	finished_recv_requests_ = NULL;
	mpi_comm_ = MPI_COMM_NULL;
#else
	utility_exit_with_message( "AsyncMPITemperingBase requires MPI build" );
#endif
	allocate_buffers();
}

AsyncMPITemperingBase::AsyncMPITemperingBase( AsyncMPITemperingBase const& other ) :
	Parent( other ),
	rank_( other.rank_ ),
	comm_size_( other.comm_size_ ),
	ready_for_exchange_( other.ready_for_exchange_ ),
	finished_( other.finished_ ),
	recv_level_buffer_( nullptr ),
	send_level_buffer_( nullptr ),
	swap_invite_send_buffer_( nullptr ),
	swap_invite_recv_buffer_( nullptr ),
	finished_recv_buffer_( nullptr ),
	finished_send_buffer_( nullptr )
{
#ifdef USEMPI
	recv_level_requests_ = NULL;
	swap_invite_send_requests_ = NULL;
	finished_recv_requests_ = NULL;
	mpi_comm_ = MPI_COMM_NULL;
#else
	utility_exit_with_message( "AsyncMPITemperingBase requires MPI build" );
#endif
	allocate_buffers();
}

AsyncMPITemperingBase& AsyncMPITemperingBase::operator=( AsyncMPITemperingBase const& other ) {
	if ( &other == this ) return *this;
	Parent::operator=( other );
	rank_ = other.rank_;
	comm_size_ = other.comm_size_;
	ready_for_exchange_ = other.ready_for_exchange_;
	finished_ = other.finished_;
	recv_level_buffer_ = nullptr;
	send_level_buffer_ = nullptr;
	swap_invite_send_buffer_ = nullptr;
	swap_invite_recv_buffer_ = nullptr;
	finished_recv_buffer_ = nullptr;
	finished_send_buffer_ = nullptr;
#ifdef USEMPI
	set_mpi_comm( other.mpi_comm() );
#endif
	allocate_buffers();
	return *this;
}

AsyncMPITemperingBase::~AsyncMPITemperingBase() {
	deallocate_buffers();
}

void AsyncMPITemperingBase::allocate_buffers() {
	if ( !comm_size_ ) return;
	send_level_buffer_ = new int[mpi_size_XCHANGE_READY];
	swap_invite_recv_buffer_ = new int[ mpi_size_SWAP_INVITE ];
	finished_send_buffer_ = new int[ mpi_size_FINISHED ];
	if ( rank() == 0 ) {
		recv_level_buffer_ = new int[mpi_size_XCHANGE_READY*comm_size_];
		swap_invite_send_buffer_ = new int[ mpi_size_SWAP_INVITE*comm_size_ ];
		finished_recv_buffer_ = new int[mpi_size_FINISHED*comm_size_];
#ifdef USEMPI
		recv_level_requests_ = new MPI_Request[ comm_size_ ];
		swap_invite_send_requests_ = new MPI_Request[ comm_size_ ];
		finished_recv_requests_ = new MPI_Request[ comm_size_ ];
#endif
	};
}

void AsyncMPITemperingBase::deallocate_buffers() {
	if ( send_level_buffer_ ) delete [] send_level_buffer_;
	if ( swap_invite_recv_buffer_ ) delete [] swap_invite_recv_buffer_;
	if ( recv_level_buffer_ ) delete [] recv_level_buffer_;
	if ( swap_invite_send_buffer_ ) delete [] swap_invite_send_buffer_;
	if ( finished_send_buffer_ ) delete [] finished_send_buffer_;
	if ( finished_recv_buffer_ ) delete [] finished_recv_buffer_;
	recv_level_buffer_ = nullptr;
	send_level_buffer_ = nullptr;
	swap_invite_send_buffer_ = nullptr;
	swap_invite_recv_buffer_ = nullptr;
	finished_send_buffer_ = nullptr;
	finished_recv_buffer_ = nullptr;
#ifdef USEMPI
	if ( swap_invite_send_requests_ ) delete [] swap_invite_send_requests_;
	if ( recv_level_requests_ ) delete [] recv_level_requests_;
	if ( finished_recv_requests_ ) delete [] finished_recv_requests_;
	recv_level_requests_ = NULL;
	swap_invite_send_requests_ = NULL;
	finished_recv_requests_ = NULL;
#endif
}

/// @brief callback executed before any Monte Carlo trials
void
AsyncMPITemperingBase::initialize_simulation(
	core::pose::Pose & pose,
	protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover,
	core::Size cycle   //non-zero if trajectory is restarted
) {
	Parent::initialize_simulation( pose, metropolis_hastings_mover, cycle );
#ifdef USEMPI
	set_mpi_comm( jd2::current_mpi_comm() );
#endif
	allocate_buffers();
	start_listening();
}

void AsyncMPITemperingBase::initialize_simulation(
	core::pose::Pose & pose,
	protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover,
	core::Size level,
	core::Real temperature,
	core::Size cycle
) {
	Parent::initialize_simulation( pose, metropolis_hastings_mover, level, temperature, cycle );
#ifdef USEMPI
	set_mpi_comm( jd2::current_mpi_comm() );
#endif
	allocate_buffers();
	start_listening();
}

void AsyncMPITemperingBase::listen_for_rank( int MPI_ONLY( other_rank ) ) {
#ifdef USEMPI
	runtime_assert( rank() == 0 );
	int err = MPI_Irecv( &recv_level_buffer_[other_rank*mpi_size_XCHANGE_READY],
		mpi_size_XCHANGE_READY,
		MPI_INT,
		other_rank,
		mpi_XCHANGE_READY,
		mpi_comm(),
		&recv_level_requests_[ other_rank ]
	);
	if ( err != MPI_SUCCESS ) tr.Error << "error in listen_for_rank " << other_rank << std::endl;
#endif
}

void AsyncMPITemperingBase::start_listening() {
	//start listening for mpi_XCHANGE_READY messages
	if ( rank() == 0 ) {
		for ( int other_rank = 0; other_rank < comm_size_; other_rank++ ) {
			listen_for_rank( other_rank );
#ifdef USEMPI
			int err = MPI_Irecv( &finished_recv_buffer_[ other_rank*mpi_size_FINISHED ],
				mpi_size_FINISHED,
				MPI_INT,
				other_rank,
				mpi_FINISHED,
				mpi_comm(),
				&finished_recv_requests_[ other_rank ]
			);
			if ( err != MPI_SUCCESS ) tr.Error << "error in when listening for mpi_FINISHED of " << other_rank << std::endl;
#endif
		}
	}
}

void
AsyncMPITemperingBase::finalize_simulation(
	pose::Pose& pose,
	protocols::canonical_sampling::MetropolisHastingsMover const & mhm
) {
	Parent::finalize_simulation( pose, mhm );
	//if replica runs end at different times that leads to hanging in output for some reason
	// this seems to help.
#ifdef USEMPI
	MPI_Barrier( mpi_comm() );
#endif
	deallocate_buffers();
}

void AsyncMPITemperingBase::send_swap_invites( ExchangeRequest const& request1, ExchangeRequest const& request2 ) {
	runtime_assert( rank() == 0 );
	tr.Trace << "send swap invite to " << request1.rank_ << " to partner with " << request2.rank_ << std::endl;
	swap_invite_send_buffer_[ request1.rank_*mpi_size_SWAP_INVITE ] = request2.rank_;
	runtime_assert( mpi_size_SWAP_INVITE == 1 );
#ifdef USEMPI
	int err = MPI_Isend( &swap_invite_send_buffer_[ request1.rank_*mpi_size_SWAP_INVITE ],
		mpi_size_SWAP_INVITE,
		MPI_INT,
		request1.rank_,
		mpi_SWAP_INVITE,
		mpi_comm(),
		&swap_invite_send_requests_[ request1.rank_ ]
	);
	if ( err != MPI_SUCCESS ) tr.Error << "ERROR for send_swap_invites to " << request1.rank_ << std::endl;
#endif
}

void AsyncMPITemperingBase::ExchangeRequest::show( std::ostream& os ) const {
	os << "ExchangeRequest rank " << rank_ << " for (" << level_ << ":" << requested_level_ << ")";
}

std::ostream& operator<< (std::ostream& os, AsyncMPITemperingBase::ExchangeRequest const& er ) {
	er.show( os );
	return os;
}

void AsyncMPITemperingBase::remove_cycles_from_request_list() {
	typedef utility::vector1< std::pair< ExchangeRequest, ExchangeRequest > > MatchList;
	MatchList forced_requests;
	for ( ExchangeRequestMap::const_iterator it = exchange_requests_.begin(); it != exchange_requests_.end(); ++it ) {
		core::Size start_level( it->first );
		core::Size second_level_of_cycle( 0 );
		core::Size next_level = start_level;
		core::Size steps( 0 );
		while ( next_level ) {
			++steps;
			//visited[ next_level ] = true;
			ExchangeRequestMap::const_iterator match = exchange_requests_.find( next_level );
			if ( match != exchange_requests_.end() ) {
				next_level = match->second.requested_level_;
				if ( !second_level_of_cycle ) second_level_of_cycle = next_level;
				if ( next_level == start_level ) {
					runtime_assert( steps > 2 ); //we should have removed a pair in previous step
					tr.Warning << "cycle detected in request list starting at " << start_level << std::endl;
					tr.Warning << "destroy cycle by forcing ExchangeRequest of " << start_level << " with " << it->second.requested_level_ << std::endl;
					ExchangeRequest const& second_request = exchange_requests_[ second_level_of_cycle ];
					ExchangeRequest new_request( second_request.rank_, second_level_of_cycle, start_level );
					//overwrite the second request to have it point back to the first...
					exchange_requests_[ second_level_of_cycle ] = new_request;
					forced_requests.push_back( std::make_pair( it->second, new_request ) );
					next_level = 0;
				}
			} else {
				next_level = 0;
			}
		}
	}
	process_match_list( forced_requests );
}

void AsyncMPITemperingBase::process_match_list( MatchList const& matches ) {
	for ( auto const & matche : matches ) {
		exchange_requests_.erase( exchange_requests_.find( matche.first.level_ ) );
		exchange_requests_.erase( exchange_requests_.find( matche.second.level_ ) );
		send_swap_invites( matche.first, matche.second );
		send_swap_invites( matche.second, matche.first );
	}
}

void AsyncMPITemperingBase::process_matching_exchange_requests() {
	//any matching ExchangeRequests ?
	runtime_assert( rank() == 0 );

	MatchList matches;
	for ( ExchangeRequestMap::const_iterator it = exchange_requests_.begin(); it != exchange_requests_.end(); ++it ) {
		core::Size level( it->first );
		runtime_assert( it->first == it->second.level_ );
		core::Size req_level( it->second.requested_level_ );
		if ( level > req_level ) continue;
		ExchangeRequestMap::const_iterator match = exchange_requests_.find( req_level );
		if ( match != exchange_requests_.end() ) {
			if ( match->second.requested_level_ != level ) continue;
			tr.Debug << "found exchange match : " << match->second << " and " << it->second << std::endl;
			matches.push_back( std::make_pair( match->second, it->second ) );
		}
	}
	process_match_list( matches );
	remove_cycles_from_request_list();
}

bool AsyncMPITemperingBase::receive_exchange_requests() {
	//any new ExchangeRequests ?
	runtime_assert( rank() == 0 );
	int outcount( 1 );
#ifdef USEMPI
	runtime_assert( comm_size_ < 1000 );
	MPI_Status statuses[ 1000 ]; //avoid dynamic memoery allocation
	int finished_ranks[ 1000 ];
	int err = MPI_Testsome( comm_size_, recv_level_requests_, &outcount, finished_ranks, statuses );
	if ( err != MPI_SUCCESS ) tr.Error << "ERROR in receive_exchange_requests() " << std::endl;
	for ( int i = 0; i<outcount; ++i ) {
		int other_rank = finished_ranks[ i ];
		int other_level = recv_level_buffer_[ other_rank*mpi_size_XCHANGE_READY ];
		int requested_level = recv_level_buffer_[ other_rank*mpi_size_XCHANGE_READY + 1 ];
		exchange_requests_[ other_level ] = ExchangeRequest( other_rank, other_level, requested_level );
		tr.Trace << "NEW " << ExchangeRequest( other_rank, other_level, requested_level ) << std::endl;
		listen_for_rank( other_rank );
	}
	if (outcount) {
		tr.Trace << "-------------------------------------------" << std::endl;
		tr.Trace << "CURRENT Xchange List:... " << std::endl;
		for ( ExchangeRequestMap::const_iterator it = exchange_requests_.begin(); it != exchange_requests_.end(); ++it ) {
			tr.Trace << it->first << ": " << it->second << std::endl;
		}
		tr.Trace << "-------------------------------------------" << std::endl;
	}
#endif
	return outcount > 0;
}

void
AsyncMPITemperingBase::find_exchange_partner( int& partner, bool& is_master ) {
	partner = swap_invite_recv_buffer_[ 0 ];
	is_master = rank() < partner;
}

bool AsyncMPITemperingBase::finished_simulation( core::Size trials, core::Size ntrials ) {
#ifdef USEMPI
	if ( !finished_ && Parent::finished_simulation( trials, ntrials ) ) {
		finished_ = true;
		int const master_rank( 0 );
		MPI_Isend( finished_send_buffer_, mpi_size_FINISHED, MPI_INT, master_rank, mpi_FINISHED, mpi_comm(), &finished_send_request_ );
		MPI_Irecv( finished_send_buffer_, mpi_size_FINISHED, MPI_INT, master_rank, mpi_FINISHED_ACK, mpi_comm(), &finished_ack_recv_request_ );
	}
	if ( rank() == 0 ) {
		int flag;
		MPI_Status statuses[ 1000 ]; //avoid dynamic memoery allocation
		runtime_assert( comm_size_ < 1000 );
		MPI_Testall( comm_size_, finished_recv_requests_, &flag, statuses );
		if ( flag ) {
			int dummy( 1 );
			MPI_Request send_ack_request[ 1000 ];
			for ( int i = 1; i<comm_size_; ++i ) {
				MPI_Isend( &dummy, 1, MPI_INT, i, mpi_FINISHED_ACK, mpi_comm(), &send_ack_request[ i-1 ] );
			}
			MPI_Waitall( comm_size_-1, send_ack_request, statuses );
			return true;
		}
	}
	if ( finished_ && rank() > 0 ) {
		int flag( 1 );
		MPI_Status status;
		MPI_Test( &finished_ack_recv_request_, &flag, &status );
		return flag;
	}
	return false;
#else
	return Parent::finished_simulation( trials, ntrials );
#endif
}

bool AsyncMPITemperingBase::time_for_temp_move() {

	//are we initialized ??
	//check that buffers are allocated, just check one representative for all
	runtime_assert( send_level_buffer_ );

	if ( rank() == 0 ) {
		bool new_requests = receive_exchange_requests();
		if ( new_requests ) process_matching_exchange_requests();
	}

	if ( !ready_for_exchange_ ) {
		if ( Parent::time_for_temp_move() ) {
			ready_for_exchange_ = true;
			core::Size next_level_for_exchange( next_exchange_level() );
			runtime_assert( next_level_for_exchange > 0 && (int) next_level_for_exchange <= comm_size_ );
			tr.Trace << rank() << " at level " << temperature_level() << " is ready for exchange with " << next_level_for_exchange << std::endl;
			runtime_assert( mpi_size_XCHANGE_READY == 2 );
			send_level_buffer_[0]=temperature_level();
			send_level_buffer_[1]=next_level_for_exchange;
#ifdef USEMPI
			int const master_rank( 0 );
			int err1 = MPI_Isend(	send_level_buffer_,	mpi_size_XCHANGE_READY,	MPI_INT, master_rank,	mpi_XCHANGE_READY, mpi_comm(), &send_level_request_	);
			int err2 = MPI_Irecv(	swap_invite_recv_buffer_,	mpi_size_SWAP_INVITE,	MPI_INT, master_rank,	mpi_SWAP_INVITE, mpi_comm(), &swap_invite_recv_request_ );
			if ( err1 != MPI_SUCCESS ) tr.Error << "ERROR in MPI_Isend of ready_for_temp_move() " << std::endl;
			if ( err2 != MPI_SUCCESS ) tr.Error << "ERROR in MPI_Irecv of ready_for_temp_move() " << std::endl;
#endif
		}
	}
	if ( ready_for_exchange_ ) {
#ifdef USEMPI
		int flag( 1 );
		MPI_Status status;
		int err = MPI_Test( &swap_invite_recv_request_, &flag, &status );
		if ( err != MPI_SUCCESS ) tr.Error << "ERROR in 1st MPI_Test of ready_for_temp_move() " << std::endl;
		if ( flag ) {
			tr.Trace << "received command to attempt swap with " << swap_invite_recv_buffer_[ 0 ] << std::endl;
			//sanity check, the send_level_request should be completed a long time ago, since we already got the answer...
			int err = MPI_Test( &send_level_request_, &flag, &status );
			if ( err != MPI_SUCCESS ) tr.Error << "ERROR in 2nd MPI_Test of ready_for_temp_move() " << std::endl;
			runtime_assert( flag );
			//... sanity is served.

			ready_for_exchange_ = false;
			reset_temp_counter();
			return true; //exchange partner is now in swap_invite_recv_buffer_
		}
#endif
	}
	return false;
}

#ifdef USEMPI
void AsyncMPITemperingBase::set_mpi_comm( MPI_Comm const& mpi_comm ) {
	if ( mpi_comm != MPI_COMM_NULL ) {
		MPI_Comm_dup( mpi_comm, &mpi_comm_ );
		MPI_Comm_rank( mpi_comm_, &rank_ );
		int communicator_size;
		MPI_Comm_size( mpi_comm_, &communicator_size );
		comm_size_ = communicator_size;
		if ( communicator_size != (int) n_temp_levels() ) {
			std::ostringstream os;
			os << "For HamiltonianExchange the number of exchange cells " << n_temp_levels()
				 << "\n has to be consistent with the option -run:n_replica "
				 << communicator_size;
			utility_exit_with_message( os.str() );
		}
	} else {
		mpi_comm_ = MPI_COMM_NULL;
	}
}
#endif

} //moves
} //protocols

