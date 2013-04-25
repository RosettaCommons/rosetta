// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/mpi_util.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <utility/mpi_util.hh>
#include <utility/SimulateMPI.hh>

namespace utility {


#ifdef USEMPI

#define MPI_ONLY(x) x

int
mpi_rank() {
	int return_val( 0 );
	MPI_Comm_rank( MPI_COMM_WORLD, & return_val);/* get current process id */
	return return_val;
}


int
mpi_nprocs()
{
	int return_val( 0 );
	MPI_Comm_size( MPI_COMM_WORLD, & return_val);/* get number of processes */
	return return_val;
}

std::string
receive_string_from_node(
	int source
) {
	std::string return_val;
	int len( 0 );
	int tag( 1 );
	MPI_Status stat;
	MPI_Recv( &len, 1, MPI_INT, source, tag, MPI_COMM_WORLD, & stat );
	char * str = new char[ len + 1 ];
	str[ len ] = '\0'; // ? do I need null terminated strings?
	MPI_Recv( str, len, MPI_CHAR, source, tag, MPI_COMM_WORLD, & stat );
	return_val = std::string( str, len );
	delete [] str;
	return return_val;
}

void
send_string_to_node(
	int destination,
	std::string const & message
) {
	int tag( 1 );
	int len( message.size() );
	MPI_Send( &len, 1, MPI_INT, destination, tag, MPI_COMM_WORLD );
	MPI_Send( const_cast< char * > (message.c_str()), len, MPI_CHAR, destination, tag, MPI_COMM_WORLD );
}


char
receive_char_from_node(
	int source
) {
	char return_val = 0;
	int tag( 1 );
	MPI_Status stat;
	MPI_Recv( &return_val, 1, MPI_CHAR, source, tag, MPI_COMM_WORLD, & stat );
	return return_val;
}

void
send_char_to_node(
	int destination,
	char message
) {
	int tag( 1 );
	MPI_Send( &message, 1, MPI_CHAR, destination, tag, MPI_COMM_WORLD );
}

int
receive_integer_from_node(
	int source
) {
	int return_val(0);
	int tag( 1 );
	MPI_Status stat;
	MPI_Recv( &return_val, 1, MPI_INT, source, tag, MPI_COMM_WORLD, & stat );
	return return_val;
}


void
send_integer_to_node(
	int destination,
	int message
) {
	int tag( 1 );
	MPI_Send( &message, 1, MPI_INT, destination, tag, MPI_COMM_WORLD );
}


utility::vector1< int >
receive_integers_from_node(
	int source
) {
	utility::vector1< int > return_val;
	int len( 0 );
	int tag( 1 );
	MPI_Status stat;
	MPI_Recv( &len, 1, MPI_INT, source, tag, MPI_COMM_WORLD, & stat );
	if ( len != 0 ) {
		return_val.resize( len );
		int * intarray = new int[ len ];
		MPI_Recv( intarray, len, MPI_INT, source, tag, MPI_COMM_WORLD, & stat );
		for ( int ii = 0; ii < len; ++ii ) return_val[ ii + 1 ] = intarray[ ii ];
		delete [] intarray;
	}
	return return_val;
}


void
send_integers_to_node(
	int destination,
	utility::vector1< int > const & message
) {
	int tag( 1 );
	int len( message.size() );
	MPI_Send( &len, 1, MPI_INT, destination, tag, MPI_COMM_WORLD );
	if ( len != 0 ) {
		MPI_Send( const_cast< int * > (&message[1]), len, MPI_INT, destination, tag, MPI_COMM_WORLD );
	}
}

double
receive_double_from_node(
	int source
) {
	double return_val(0);
	int tag( 1 );
	MPI_Status stat;
	MPI_Recv( &return_val, 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, & stat );
	return return_val;
}


void
send_double_to_node(
	int destination,
	double message
) {
	int tag( 1 );
	MPI_Send( &message, 1, MPI_DOUBLE, destination, tag, MPI_COMM_WORLD );
}



utility::vector1< double >
receive_doubles_from_node(
	int source
) {
	utility::vector1< double > return_val;
	int len( 0 );
	int tag( 1 );
	MPI_Status stat;
	MPI_Recv( &len, 1, MPI_INT, source, tag, MPI_COMM_WORLD, & stat );
	if ( len != 0 ) {
		return_val.resize( len );
		double * doublearray = new double[ len ];
		MPI_Recv( doublearray, len, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, & stat );
		for ( int ii = 0; ii < len; ++ii ) return_val[ ii + 1 ] = doublearray[ ii ];
		delete [] doublearray;
	}
	return return_val;
}


void
send_doubles_to_node(
	int destination,
	utility::vector1< double > const & message
) {
	int tag( 1 );
	int len( message.size() );
	MPI_Send( &len, 1, MPI_INT, destination, tag, MPI_COMM_WORLD );
	if ( len != 0 ) {
		MPI_Send( const_cast< double * > (&message[1]), len, MPI_DOUBLE, destination, tag, MPI_COMM_WORLD );
	}
}

////////////////////////////
#else // USEMPI is not used
///////////////////////////

#define MPI_ONLY(x)

int
mpi_rank() {
	if(SimulateMPI::simulate_mpi()){
		return SimulateMPI::mpi_rank();
	} else {
		int return_val( 0 );
		return return_val;
	}
}


int
mpi_nprocs()
{
	if(SimulateMPI::simulate_mpi()){
		return SimulateMPI::mpi_nprocs();
	} else {
		int return_val( 0 );
		return return_val;
	}
}

std::string
receive_string_from_node(
	int source
) {
	if(SimulateMPI::simulate_mpi()){
		return SimulateMPI::receive_string_from_node(source);
	} else {
		std::string return_val;
		return return_val;
	}
}

void
send_string_to_node(
	int destination,
	std::string const & message)
{
	if(SimulateMPI::simulate_mpi()){
		SimulateMPI::send_string_to_node(destination, message);
	} else {
		return;
	}
}

char
receive_char_from_node(
	int source
) {
	if(SimulateMPI::simulate_mpi()){
		return SimulateMPI::receive_char_from_node(source);
	} else {
		char return_val = 0;
		return return_val;
	}
}

void
send_char_to_node(
	int destination,
	char message)
{
	if(SimulateMPI::simulate_mpi()){
		SimulateMPI::send_char_to_node(destination, message);
	} else {
		return;
	}
}

int
receive_integer_from_node(
	int source
) {
	if(SimulateMPI::simulate_mpi()){
		return SimulateMPI::receive_integer_from_node(source);
	} else {
		int return_val = 0;
		return return_val;
	}
}

void
send_integer_to_node(
	int destination,
	int message)
{
	if(SimulateMPI::simulate_mpi()){
		SimulateMPI::send_integer_to_node(destination, message);
	} else {
		return;
	}
}

utility::vector1< int >
receive_integers_from_node(
	int source
) {
	if(SimulateMPI::simulate_mpi()){
		return SimulateMPI::receive_integers_from_node(source);
	} else {
		utility::vector1< int > return_val;
		return return_val;
	}
}

void
send_integers_to_node(
	int destination,
	utility::vector1< int > const & message)
{
	if(SimulateMPI::simulate_mpi()){
		SimulateMPI::send_integers_to_node(destination, message);
	} else {
		return;
	}
}



double
receive_double_from_node(
	int source
) {
	if(SimulateMPI::simulate_mpi()){
		return SimulateMPI::receive_double_from_node(source);
	} else {
		double return_val = 0;
		return return_val;
	}
}

void
send_double_to_node(
	int destination,
	double message)
{
	if(SimulateMPI::simulate_mpi()){
		SimulateMPI::send_double_to_node(destination, message);
	} else {
		return;
	}
}

utility::vector1< double >
receive_doubles_from_node(
	int source
) {
	if(SimulateMPI::simulate_mpi()){
		return SimulateMPI::receive_doubles_from_node(source);
	} else {
		utility::vector1< double > return_val;
		return return_val;
	}
}

void
send_doubles_to_node(
	int destination,
	utility::vector1< double > const & message)
{
	if(SimulateMPI::simulate_mpi()){
		SimulateMPI::send_doubles_to_node(destination, message);
	} else {
		return;
	}
}


#endif // USEMPI



}

