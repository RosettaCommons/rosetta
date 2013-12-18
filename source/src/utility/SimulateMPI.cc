// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/SimulateMPI.cc
/// @brief
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#include <utility/SimulateMPI.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>

namespace utility {

std::string
msg_name( simulate_mpi_message_type msg_type ) {
	switch( msg_type ) {
	case smpi_char		: return "smpi_char";
	case smpi_integer : return "smpi_integer";
	case smpi_string	: return "smpi_string";
	case smpi_double	: return "smpi_double";
	case smpi_integers: return "smpi_integers";
	case smpi_doubles : return "smpi_doubles";
	}
	return "unknown message type";
}

SimulateMPIMessage::SimulateMPIMessage() :
	index_( 0 ),
	src_( 0 ),
	dst_( 0 ),
	processed_( false )
{}

void SimulateMPIMessage::src( platform::Size source ) { src_ = source; }
void SimulateMPIMessage::dst( platform::Size destination ) { dst_ = destination; }

void SimulateMPIMessage::mark_as_processed() { processed_ = true; }
bool SimulateMPIMessage::processed() const { return processed_; }

void SimulateMPIMessage::set_index( platform::Size setting ) { index_ = setting; }

void SimulateMPIMessage::set_char_msg( char setting ) {
	msg_type_ = smpi_char;
	char_msg_ = setting;
}

void SimulateMPIMessage::set_integer_msg( int setting ) {
	msg_type_ = smpi_integer;
	integer_msg_ = setting;
}

void SimulateMPIMessage::set_string_msg( std::string const & setting ) {
	msg_type_ = smpi_string;
	string_msg_ = setting;
}

void SimulateMPIMessage::set_double_msg( double setting ) {
	msg_type_ = smpi_double;
	double_msg_ = setting;
}

void SimulateMPIMessage::set_integers_msg( utility::vector1< int > const & setting ) {
	msg_type_ = smpi_integers;
	integers_msg_ = setting;
}

void SimulateMPIMessage::set_doubles_msg( utility::vector1< double > const & setting ) {
	msg_type_ = smpi_doubles;
	doubles_msg_ = setting;
}

// initialize private static data
SimulateMPIData * SimulateMPI::simulation_( 0 );
int SimulateMPI::rank_( 0 );

SimulateMPIData::SimulateMPIData( platform::Size nprocs ) :
	mpi_nprocs_( nprocs ),
	nmessages_( 0 ),
	messages_for_node_( nprocs ),
	messages_from_node_( nprocs ),
	messages_( nprocs )
{
	for ( platform::Size ii = 0; ii < nprocs; ++ii ) {
		messages_[ ii ].resize( nprocs );
	}
}

void SimulateMPIData::queue_message( SimulateMPIMessageOP msg )
{
	++nmessages_;
	msg->set_index( nmessages_ );
	all_messages_.push_back( msg );
	messages_for_node_[ msg->dst() ].push_back( msg );
	messages_from_node_[ msg->src() ].push_back( msg );
	messages_[ msg->dst() ][ msg->src() ].push_back( msg );
}

SimulateMPIMessageOP
SimulateMPIData::pop_next_message_for_node_of_type( platform::Size dst, simulate_mpi_message_type msg_type )
{
	SimulateMPIMessageOP most_recent_message( 0 );
	platform::Size most_recent_index( 0 );

	for ( platform::Size ii = 0; ii < mpi_nprocs_; ++ii ) {
		clear_processed_msgs( messages_[ dst ][ ii ] );
		if ( messages_[ dst ][ ii ].empty() ) continue;
		SimulateMPIMessageOP iimsg = messages_[ dst ][ ii ].front();
		if ( iimsg->msg_type() == msg_type ) {
			if ( most_recent_index == 0 || iimsg->index() < most_recent_index ) {
				most_recent_message = iimsg;
				most_recent_index = iimsg->index();
			}
		}
	}
	if ( most_recent_index == 0 ) {
		throw excn::EXCN_Msg_Exception( "SimulatedMPIData could not pop a message of type '" + msg_name( msg_type ) + "' for node " + to_string( dst ) );
	}
	most_recent_message->mark_as_processed();
	return most_recent_message;
}

SimulateMPIMessageOP
SimulateMPIData::pop_next_message_of_type(
	platform::Size dst,
	platform::Size src,
	simulate_mpi_message_type msg_type
)
{
	clear_processed_msgs( messages_[ dst ][ src ] );
	SimulateMPIMessageOP msg;
	if ( messages_[ dst ][ src ].empty() ) {
		throw excn::EXCN_Msg_Exception( "Could not retrieve a " + msg_name( msg_type ) + " message to " + to_string( dst ) + " from " + to_string( src ) + "; message queue is empty" );
	}

	msg = messages_[ dst ][ src ].front();
	if ( msg->msg_type() != msg_type ) {
		throw excn::EXCN_Msg_Exception( "Could not retrieve a " + msg_name( msg_type ) + " message to " + to_string( dst ) + " from " + to_string( src ) + "; next message is of type " + msg_name( msg->msg_type() ) );
	}

	msg->mark_as_processed();
	return msg;
}

/// @details Erase elements of the queue that have already been processed
/// until we arrive at an element that has not yet been processed and then
/// return.
void SimulateMPIData::clear_processed_msgs( MsgQueue & message_queue )
{
	for ( MsgQueue::iterator iter = message_queue.begin(),
			iter_end = message_queue.end();
			iter != iter_end; /*noinc*/ ) {
		if ( (*iter)->processed() ) {
			MsgQueue::iterator next_iter = iter;
			++next_iter;
			message_queue.erase( iter );
			iter = next_iter;
		} else {
			return;
		}
	}
}

void
SimulateMPI::initialize_simulation( int nprocs ) {
	simulation_ = new SimulateMPIData( nprocs );
	rank_ = 0;
}

bool
SimulateMPI::simulate_mpi() {
	return simulation_;
}


void
SimulateMPI::set_mpi_rank(int value) {
	assert(simulation_);
	assert(0 <= value);
	assert(value < mpi_nprocs());
	rank_ = value;
}


int
SimulateMPI::mpi_rank() {
	assert(simulation_);
	return rank_;
}

int
SimulateMPI::mpi_nprocs() {
	return simulation_->mpi_nprocs();
}

int
SimulateMPI::receive_integer_from_anyone()
{
	SimulateMPIMessageOP msg = simulation_->pop_next_message_for_node_of_type( rank_, smpi_integer );
	return msg->integer_msg();
}

std::string
SimulateMPI::receive_string_from_node(
	int source
) {
	assert(simulation_);
	assert(0 <= source);
	assert(source < mpi_nprocs());

	SimulateMPIMessageOP msg = simulation_->pop_next_message_of_type( rank_, source, smpi_string );
	return msg->string_msg();
}

void
SimulateMPI::send_string_to_node(
	int destination,
	std::string const & message
) {
	assert(simulation_);
	assert(0 <= destination);
	assert(destination < mpi_nprocs());

	SimulateMPIMessageOP msg = new SimulateMPIMessage;
	msg->src( rank_ );
	msg->dst( destination );
	msg->set_string_msg( message );
	simulation_->queue_message( msg );
}


char
SimulateMPI::receive_char_from_node(
	int source
) {
	assert(simulation_);
	assert(0 <= source);
	assert(source < mpi_nprocs());

	SimulateMPIMessageOP msg = simulation_->pop_next_message_of_type( rank_, source, smpi_char );
	return msg->char_msg();
}

void
SimulateMPI::send_char_to_node(
	int destination,
	char message
) {
	assert(simulation_);
	assert(0 <= destination);
	assert(destination < mpi_nprocs());

	SimulateMPIMessageOP msg = new SimulateMPIMessage;
	msg->src( rank_ );
	msg->dst( destination );
	msg->set_char_msg( message );
	simulation_->queue_message( msg );
}


int
SimulateMPI::receive_integer_from_node(
	int source
) {
	assert(simulation_);
	assert(0 <= source);
	assert(source < mpi_nprocs());

	SimulateMPIMessageOP msg = simulation_->pop_next_message_of_type( rank_, source, smpi_integer );
	return msg->integer_msg();
}

void
SimulateMPI::send_integer_to_node(
	int destination,
	int message
) {
	assert(simulation_);
	assert(0 <= destination);
	assert(destination < mpi_nprocs());

	SimulateMPIMessageOP msg = new SimulateMPIMessage;
	msg->src( rank_ );
	msg->dst( destination );
	msg->set_integer_msg( message );
	simulation_->queue_message( msg );
}

vector1< int >
SimulateMPI::receive_integers_from_node(
	int source
) {
	assert(simulation_);
	assert(0 <= source);
	assert(source < mpi_nprocs());

	SimulateMPIMessageOP msg = simulation_->pop_next_message_of_type( rank_, source, smpi_integers );
	return msg->integers_msg();
}

void
SimulateMPI::send_integers_to_node(
	int destination,
	vector1< int > const & message
) {
	assert(simulation_);
	assert(0 <= destination);
	assert(destination < mpi_nprocs());

	SimulateMPIMessageOP msg = new SimulateMPIMessage;
	msg->src( rank_ );
	msg->dst( destination );
	msg->set_integers_msg( message );
	simulation_->queue_message( msg );
}


double
SimulateMPI::receive_double_from_node(
	int source
) {
	assert(simulation_);
	assert(0 <= source);
	assert(source < mpi_nprocs());

	SimulateMPIMessageOP msg = simulation_->pop_next_message_of_type( rank_, source, smpi_double );
	return msg->double_msg();
}

void
SimulateMPI::send_double_to_node(
	int destination,
	double message
) {
	assert(simulation_);
	assert(0 <= destination);
	assert(destination < mpi_nprocs());

	SimulateMPIMessageOP msg = new SimulateMPIMessage;
	msg->src( rank_ );
	msg->dst( destination );
	msg->set_double_msg( message );
	simulation_->queue_message( msg );
}

vector1< double >
SimulateMPI::receive_doubles_from_node(
	int source
) {
	assert(simulation_);
	assert(0 <= source);
	assert(source < mpi_nprocs());

	SimulateMPIMessageOP msg = simulation_->pop_next_message_of_type( rank_, source, smpi_doubles );
	return msg->doubles_msg();
}

void
SimulateMPI::send_doubles_to_node(
	int destination,
	vector1< double > const & message
) {
	assert(simulation_);
	assert(0 <= destination);
	assert(destination < mpi_nprocs());

	SimulateMPIMessageOP msg = new SimulateMPIMessage;
	msg->src( rank_ );
	msg->dst( destination );
	msg->set_doubles_msg( message );
	simulation_->queue_message( msg );

}


} // utility
