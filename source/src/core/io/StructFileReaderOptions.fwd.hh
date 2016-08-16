// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/StructFileReaderOptions.fwd.hh
///
/// @brief
/// @author Brian D. Weitzner brian.weitzner@gmail.com

#ifndef INCLUDED_core_io_StructFileReaderOptions_FWD_HH
#define INCLUDED_core_io_StructFileReaderOptions_FWD_HH


#include <utility/pointer/owning_ptr.hh>

// C++ headers

namespace core {
namespace io {

class StructFileReaderOptions;

typedef utility::pointer::shared_ptr< StructFileReaderOptions > StructFileReaderOptionsOP;
typedef utility::pointer::shared_ptr< StructFileReaderOptions const > StructFileReaderOptionsCOP;

} // namespace io
} // namespace core


#endif // INCLUDED_core_io_StructFileReaderOptions_FWD_HH
