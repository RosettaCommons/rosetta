// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/mpi/mpi_enums.hh
/// @author Tim Jacobs
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_basic_mpi_mpi_enum_hh
#define INCLUDED_basic_mpi_mpi_enum_hh

/// @brief Tags used to tag messeges sent by MPI functions used to decide whether a slave is requesting a new job id or
///flagging as job as being a bad input
namespace basic {
namespace mpi {


enum mpi_tags {
	NEW_JOB_ID_TAG = 10,
	BAD_INPUT_TAG = 20,
	JOB_SUCCESS_TAG = 30,
	JOB_FAILURE_TAG = 40,
	REQUEST_MESSAGE_TAG = 50,
	RECEIVE_MESSAGE_TAG = 60,
	JOB_GO_TAG = 70
};



} //mpi
} //basic


#endif
