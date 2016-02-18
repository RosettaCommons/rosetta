// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/MoverStatus.hh
/// @brief  return status enum for Movers
/// @author Steven Lewis smlewi@gmail.com
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_protocols_moves_MoverStatus_hh
#define INCLUDED_protocols_moves_MoverStatus_hh

#include <string>

namespace protocols {
namespace moves {

/// @brief return status for movers - mover was successful, failed but can be retried, etc; used mostly by job dist.
// why not naming all thes flags MS_XXXX ???
// probably because the MS part is encoded in the type! a related question: why
// MS_SUCCESS and not just SUCCESS?
enum MoverStatus {
	MS_SUCCESS = 0,
	FAIL_RETRY,
	FAIL_DO_NOT_RETRY,
	FAIL_BAD_INPUT,
	FAIL,


	//Alternative names of enums.
	MS_FAIL_RETRY = FAIL_RETRY,
	MS_FAIL_DO_NOT_RETRY = FAIL_DO_NOT_RETRY,
	MS_FAIL_BAD_INPUT = FAIL_BAD_INPUT,
	MS_FAIL = FAIL
};

MoverStatus mstype_from_name( std::string const & name );

}//moves
}//protocols

#endif //INCLUDED_protocols_moves_MoverStatus_HH
