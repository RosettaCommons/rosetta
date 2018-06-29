// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/moves/MoverStatus.cc
/// @brief
/// @author Nobuyasu Koga (nobuyasu@uw.edu)

#include <protocols/moves/MoverStatus.hh>
#include <utility/exit.hh>


namespace protocols {
namespace moves {

MoverStatus mstype_from_name( std::string const & name )
{
	MoverStatus ms;
	if ( name == "MS_SUCCESS" ) {
		ms = MS_SUCCESS;
	} else if ( name == "FAIL_RETRY" ) {
		ms = FAIL_RETRY;
	} else if ( name == "FAIL_DO_NOT_RETRY" ) {
		ms = FAIL_DO_NOT_RETRY;
	} else if ( name == "FAIL_BAD_INPUT" ) {
		ms = FAIL_BAD_INPUT;
	} else {
		//ms = FAIL_BAD_INPUT;
		utility_exit_with_message( "Invalid name of mover status !" );
	}
	return ms;
}

}//moves
}//protocols
