// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/SimulateMPI.hh
/// @brief
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_utility_SimulateMPI_hh
#define INCLUDED_utility_SimulateMPI_hh


// MPI Headers have to be #included first
#ifdef USEMPI
#include <mpi.h>
#endif

// Unit headers
#include <utility/SimulateMPI.fwd.hh>

// Package headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <list>
#include <string>

namespace utility {

/// @detail This is for unit testing mpi code in a single processor. The
///idea is to buffer the messages in the SimulateMPIData stored in the
///SimulateMPI. To use this class, call initialize_simulation( nprocs ),
///then set the mpi rank can be set manually, and the functions in
///mpi_util are usable. By setting the mpi_rank to a different
///processor, other messages can be sent and received. See
///test/utility/simulate_mpi.cxxtest for examples.

enum simulate_mpi_message_type {
	smpi_char = 1,
  smpi_integer,
	smpi_string,
	smpi_double,
	smpi_integers,
	smpi_doubles
};

class SimulateMPIMessage : public pointer::ReferenceCount {
public:
	SimulateMPIMessage();
	void src( platform::Size source );
	void dst( platform::Size destination );
	platform::Size src() const { return src_; }
	platform::Size dst() const { return dst_; }
	void mark_as_processed();
	bool processed() const;

	/// @brief the SimulateMPIData class is responsible for setting the index of a message
	void set_index( platform::Size setting );

	void set_char_msg( char setting );
	void set_integer_msg( int setting );
	void set_string_msg( std::string const & setting );
	void set_double_msg( double setting );
	void set_integers_msg( utility::vector1< int > const & setting );
	void set_doubles_msg( utility::vector1< double > const & setting );

	platform::Size index() const { return index_; }
	simulate_mpi_message_type msg_type() const { return msg_type_; }
	char char_msg() const { return char_msg_; }
	int  integer_msg() const { return integer_msg_; }
	std::string const & string_msg() const { return string_msg_; }
	double double_msg() const { return double_msg_; }
	utility::vector1< int > const & integers_msg() const { return integers_msg_; }
	utility::vector1< double > const & doubles_msg() const { return doubles_msg_; }

private:
	platform::Size index_;
	platform::Size src_;
	platform::Size dst_;
	bool processed_;
	simulate_mpi_message_type msg_type_;
	char char_msg_;
	int  integer_msg_;
	std::string string_msg_;
	double double_msg_;
	utility::vector1< int > integers_msg_;
	utility::vector1< double > doubles_msg_;
};


class SimulateMPIData {
public:
	typedef std::list< SimulateMPIMessageOP > MsgQueue;
public:
	SimulateMPIData( platform::Size nprocs );

	int mpi_nprocs() const { return	mpi_nprocs_; }

	void queue_message( SimulateMPIMessageOP msg );

	SimulateMPIMessageOP pop_next_message_for_node_of_type( platform::Size dst, simulate_mpi_message_type msg_type );
	SimulateMPIMessageOP pop_next_message_of_type( platform::Size dst, platform::Size src, simulate_mpi_message_type msg_type );

private:
	void clear_processed_msgs( MsgQueue & );

private:

	platform::Size mpi_nprocs_;
	platform::Size nmessages_;
	MsgQueue all_messages_;
	vector0< MsgQueue > messages_for_node_;
	vector0< MsgQueue > messages_from_node_;
	vector0< vector0< MsgQueue > > messages_;

};

/// @brief singleton class storing simulated MPI state
class SimulateMPI {

private:
	// private and unimplemented
	SimulateMPI();
	SimulateMPI( SimulateMPI const & src );

public:

	static
	bool
	simulate_mpi();

	static
	void
	initialize_simulation( int nprocs );

	static
	void
	set_mpi_rank( int value );

	static
	int
	mpi_rank();

	static
	int
	mpi_nprocs();

	static
	int
	receive_integer_from_anyone();

	static
	std::string
	receive_string_from_node(int source);

	static
	void
	send_string_to_node(
		int destination,
		std::string const & message);

	static
	char
	receive_char_from_node(int source);

	static
	void
	send_char_to_node(
		int destination,
		char message);

	static
	int
	receive_integer_from_node(int source);

	static
	void
	send_integer_to_node(
		int destination,
		int message);

	static
	vector1< int >
	receive_integers_from_node(int source);

	static
	void
	send_integers_to_node(
		int destination,
		vector1< int > const & message);

	static
	double
	receive_double_from_node(int source);

	static
	void
	send_double_to_node(
		int destination,
		double message);

	static
	vector1< double >
	receive_doubles_from_node(int source);

	static
	void
	send_doubles_to_node(
		int destination,
		vector1< double > const & message);

private:

	static SimulateMPIData * simulation_;
	static int rank_;

};


} // utility

#endif
