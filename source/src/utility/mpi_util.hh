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

#ifndef INCLUDED_utility_mpi_util_hh
#define INCLUDED_utility_mpi_util_hh


// MPI Headers have to be #included first
#ifdef USEMPI
#include <mpi.h>
#endif

// Utility headers
#include <utility/vector1.hh>

/// STL Headers
// Have to go full string because of Windows PyRosetta
#include <string>

namespace utility {

int
mpi_rank();

int
mpi_nprocs();

/// @brief Use MPI to wait until some node sends an integer -- usually its own mpi_rank
/// so that it can send further messages.
int
receive_integer_from_anyone();

/// @brief Use MPI to receive a string from a particular node.
std::string
receive_string_from_node( int source );

void
send_string_to_node( int source, std::string const & message );

/// @brief Use MPI to receive a single char from a particular node.
char
receive_char_from_node( int source );

void
send_char_to_node( int destination, char message );

/// @brief Use MPI to receive a single integer from a particular node.
int
receive_integer_from_node( int source );

void
send_integer_to_node( int destination, int message );

/// @brief Use MPI to receive a vector of integers from a particular node.
utility::vector1< int >
receive_integers_from_node( int source );

void
send_integers_to_node( int destination, utility::vector1< int > const & message );

/// @brief Use MPI to receive a single double from a particular node.
double
receive_double_from_node( int source );

void
send_double_to_node( int destination, double message );

/// @brief Use MPI to receive a vector of doubles from a particular node.
utility::vector1< double >
receive_doubles_from_node( int source );

void
send_doubles_to_node( int destination, utility::vector1< double > const & message  );

}

#endif
