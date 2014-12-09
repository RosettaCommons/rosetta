// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/file_data_options.fwd.hh
///
/// @brief
/// @author Brian D. Weitzner brian.weitzner@gmail.com

#ifndef INCLUDED_core_io_pdb_file_data_options_FWD_HH
#define INCLUDED_core_io_pdb_file_data_options_FWD_HH


#include <utility/pointer/owning_ptr.hh>

// C++ headers

namespace core {
namespace io {
namespace pdb {

class FileDataOptions;

typedef utility::pointer::shared_ptr< FileDataOptions > FileDataOptionsOP;
typedef utility::pointer::shared_ptr< FileDataOptions const > FileDataOptionsCOP;

} // namespace pdb
} // namespace io
} // namespace core


#endif // INCLUDED_core_io_pdb_file_data_options_FWD_HH
