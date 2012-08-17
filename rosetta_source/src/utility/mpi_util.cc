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

#ifdef USEMPI // MPI version
	#define MPI_ONLY(x) x
#else // Non MPI Versin
	#define MPI_ONLY(x)
#endif


namespace utility {

int
mpi_rank() {
	int return_val( 0 );
#ifdef USEMPI
	MPI_Comm_rank( MPI_COMM_WORLD, & return_val);/* get current process id */
#endif
	return return_val;
}


int
mpi_nprocs()
{
	int return_val( 0 );
#ifdef USEMPI
	MPI_Comm_size( MPI_COMM_WORLD, & return_val);/* get number of processes */
#endif
	return return_val;
}

std::string
receive_string_from_node( int MPI_ONLY( source ) )
{
	std::string return_val;
#ifdef USEMPI
	int len( 0 );
	int tag( 1 );
	MPI_Status stat;
	MPI_Recv( &len, 1, MPI_INT, source, tag, MPI_COMM_WORLD, & stat );
	char * str = new char[ len + 1 ];
	str[ len ] = '\0'; // ? do I need null terminated strings?
	MPI_Recv( str, len, MPI_CHAR, source, tag, MPI_COMM_WORLD, & stat );
	return_val = std::string( str, len );
	delete [] str;
#endif
	return return_val;
}

void
send_string_to_node( int MPI_ONLY( destination ), std::string const & MPI_ONLY( message ) )
{
#ifdef USEMPI
	int tag( 1 );
	int len( message.size() );
	MPI_Send( &len, 1, MPI_INT, destination, tag, MPI_COMM_WORLD );
	MPI_Send( const_cast< char * > (message.c_str()), len, MPI_CHAR, destination, tag, MPI_COMM_WORLD );
#endif
}


char
receive_char_from_node( int MPI_ONLY( source ) )
{
	char return_val;
#ifdef USEMPI
	int tag( 1 );
	MPI_Status stat;
	MPI_Recv( &return_val, 1, MPI_CHAR, source, tag, MPI_COMM_WORLD, & stat );
#endif
	return return_val;
}

void
send_char_to_node( int MPI_ONLY( destination ), char MPI_ONLY( message ) )
{
#ifdef USEMPI
	int tag( 1 );
	MPI_Send( &message, 1, MPI_CHAR, destination, tag, MPI_COMM_WORLD );
#endif
}

int
receive_integer_from_node( int MPI_ONLY( source ) )
{
	int return_val(0);
#ifdef USEMPI
	int tag( 1 );
	MPI_Status stat;
	MPI_Recv( &return_val, 1, MPI_INT, source, tag, MPI_COMM_WORLD, & stat );
#endif
	return return_val;
}


void
send_integer_to_node( int MPI_ONLY( destination ), int MPI_ONLY( message ) )
{
#ifdef USEMPI
	int tag( 1 );
	MPI_Send( &message, 1, MPI_INT, destination, tag, MPI_COMM_WORLD );
#endif
}


utility::vector1< int >
receive_integers_from_node( int MPI_ONLY( source ) )
{
	utility::vector1< int > return_val;
#ifdef USEMPI
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
#endif
	return return_val;
}


void
send_integers_to_node( int MPI_ONLY( destination ), utility::vector1< int > const & MPI_ONLY( message ) )
{
#ifdef USEMPI
	int tag( 1 );
	int len( message.size() );
	MPI_Send( &len, 1, MPI_INT, destination, tag, MPI_COMM_WORLD );
	if ( len != 0 ) {
		MPI_Send( const_cast< int * > (&message[1]), len, MPI_INT, destination, tag, MPI_COMM_WORLD );
	}
#endif
}

double
receive_double_from_node( int MPI_ONLY( source ) )
{
	double return_val(0);
#ifdef USEMPI
	int tag( 1 );
	MPI_Status stat;
	MPI_Recv( &return_val, 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, & stat );
#endif
	return return_val;
}


void
send_double_to_node( int MPI_ONLY( destination ), double MPI_ONLY( message ) )
{
#ifdef USEMPI
	int tag( 1 );
	MPI_Send( &message, 1, MPI_DOUBLE, destination, tag, MPI_COMM_WORLD );
#endif
}



utility::vector1< double >
receive_doubles_from_node( int MPI_ONLY( source ) )
{
	utility::vector1< double > return_val;
#ifdef USEMPI
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
#endif
	return return_val;
}


void
send_doubles_to_node( int MPI_ONLY( destination ), utility::vector1< double > const & MPI_ONLY( message ) )
{
#ifdef USEMPI
	int tag( 1 );
	int len( message.size() );
	MPI_Send( &len, 1, MPI_INT, destination, tag, MPI_COMM_WORLD );
	if ( len != 0 ) {
		MPI_Send( const_cast< double * > (&message[1]), len, MPI_DOUBLE, destination, tag, MPI_COMM_WORLD );
	}
#endif
}


}

