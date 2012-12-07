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

#include <utility/SimulateMPI.fwd.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray2A.hh>
#include <string>

namespace utility {

///@detail This is for unit testing mpi code in a single processor. The
///idea is to buffer the messages in the SimulateMPIData stored in the
///SimulateMPI. To use this class, call initialize_simulation( nprocs ),
///then set the mpi rank can be set manually, and the functions in
///mpi_util are usable. By setting the mpi_rank to a different
///processor, other messages can be sent and received. See
///test/utility/simulate_mpi.cxxtest for examples.



struct SimulateMPIData {
	SimulateMPIData();
	SimulateMPIData(SimulateMPIData const &);

	int mpi_rank_;
	int mpi_nprocs_;

	ObjexxFCL::FArray2D< std::string > proc_string_buf_;
	ObjexxFCL::FArray2D< char > proc_char_buf_;
	ObjexxFCL::FArray2D< int > proc_integer_buf_;
	ObjexxFCL::FArray2D< vector1< int > > proc_integers_buf_;
	ObjexxFCL::FArray2D< double > proc_double_buf_;
	ObjexxFCL::FArray2D< vector1< double > > proc_doubles_buf_;

};

///@brief singleton class storing simulated MPI state
class SimulateMPI {

private:
	SimulateMPI();

	SimulateMPI( SimulateMPI const & src );

public:

	static
	bool
	simulate_mpi();

	static
	SimulateMPI *
	get_simulation();

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
	void
	set_mpi_nprocs( int value );

	static
	int
	mpi_nprocs();

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

};


} // utility

#endif
