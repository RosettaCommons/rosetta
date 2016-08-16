// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   apps/pilot/andrew/apl_msd.cc
/// @brief  Multistate design executable as written by apl.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


/// Core headers
#include <devel/init.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/mpi_util.hh>
#include <utility/string_util.hh>

using basic::t_info;
static THREAD_LOCAL basic::Tracer TR( "app.andrew.mpi_test", t_info );


int main( int argc, char ** argv )
{
	try {

	using namespace utility;

	devel::init( argc, argv );

	TR << "Hello World from node " << mpi_rank() << " of " << mpi_nprocs() << std::endl;

	if ( mpi_rank() == 0 ) {
		vector1< int > onetwothree( 3 ); onetwothree[ 1 ] = 1; onetwothree[ 2 ] = 2; onetwothree[ 3 ] = 3;
		for ( int ii = 1; ii < mpi_nprocs(); ++ii ) {
			send_integers_to_node( ii, onetwothree );
		}
		std::string whatup = "What's up, Doc?";
		for (int ii = 1; ii < mpi_nprocs(); ++ii ) {
			send_string_to_node( ii, whatup );
		}
		vector1< double > onetwothree_d( 3 ); onetwothree_d[ 1 ] = 1.25; onetwothree_d[ 2 ] = 2.5; onetwothree_d[ 3 ] = 3.75;
		for ( int ii = 1; ii < mpi_nprocs(); ++ii ) {
			send_doubles_to_node( ii, onetwothree_d );
		}
		for ( int ii = 1; ii < mpi_nprocs(); ++ii ) {
			send_double_to_node( ii, 1.25 * ii );
		}
	} else {
		vector1< int > message = receive_integers_from_node( 0 );
		TR << "Node " << mpi_rank() << " received integer message";
		for ( int ii = 1; ii <= message.size(); ++ii ) TR << " " << message[ ii ];
		TR << std::endl;
		std::string string_message = receive_string_from_node( 0 );
		TR << "Node " << mpi_rank() << " received string message " << string_message << std::endl;
		vector1< double > double_message = receive_doubles_from_node( 0 );
		TR << "Node " << mpi_rank() << " received double message";
		for ( int ii = 1; ii <= double_message.size(); ++ii ) TR << " " << double_message[ ii ];
		TR << std::endl;
		TR << "Node " << mpi_rank() << " recieved a single double " << receive_double_from_node( 0 ) << std::endl;
	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
