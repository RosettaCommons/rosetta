// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/match_enzdes_util/EnzdesLoopsFile.fwd.hh
///
/// @brief
/// @author Florian Richter, floric@u.washington.edu, april 2009


#ifndef INCLUDED_protocols_toolbox_match_enzdes_util_EnzdesLoopsFile_fwd_hh
#define INCLUDED_protocols_toolbox_match_enzdes_util_EnzdesLoopsFile_fwd_hh


#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace toolbox {
namespace match_enzdes_util {

class EnzdesLoopsFile;
typedef utility::pointer::shared_ptr< EnzdesLoopsFile > EnzdesLoopsFileOP;
typedef utility::pointer::shared_ptr< EnzdesLoopsFile const > EnzdesLoopsFileCOP;

class EnzdesLoopInfo;
typedef utility::pointer::shared_ptr< EnzdesLoopInfo > EnzdesLoopInfoOP;
typedef utility::pointer::shared_ptr< EnzdesLoopInfo const > EnzdesLoopInfoCOP;

} //match_enzdes_util
} //toolbox
} //protocols

#endif //
