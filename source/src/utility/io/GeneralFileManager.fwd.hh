// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file utility/io/GeneralFileManager.fwd.hh
/// @brief A singleton class for managing arbitrary files to ensure that they are loaded once and only once from disk.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_utility_io_GeneralFileManager_fwd_hh
#define INCLUDED_utility_io_GeneralFileManager_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace utility {
namespace io {

class GeneralFileManager;
class GeneralFileContents;

typedef utility::pointer::shared_ptr< GeneralFileManager > GeneralFileManagerOP;
typedef utility::pointer::shared_ptr< GeneralFileManager const > GeneralFileManagerCOP;
typedef utility::pointer::shared_ptr< GeneralFileContents > GeneralFileContentsOP;
typedef utility::pointer::shared_ptr< GeneralFileContents const > GeneralFileContentsCOP;

} //utility
} //io

#endif //INCLUDED_utility_io_GeneralFileManager_fwd_hh
