// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/chromophore/ChromophoreDataReader.fwd.hh
/// @brief ChromophoreDataReader class forward declarations header
/// @author Nina Bozhanova (nbozhanova@gmail.com)

#ifndef INCLUDED_protocols_chromophore_ChromophoreDataReader_fwd_hh
#define INCLUDED_protocols_chromophore_ChromophoreDataReader_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace chromophore {

class ChromophoreDataReader;

using ChromophoreDataReaderOP = utility::pointer::shared_ptr< ChromophoreDataReader >;
using ChromophoreDataReaderCOP = utility::pointer::shared_ptr< ChromophoreDataReader const >;

} // chromophore
} // protocols

#endif //INCLUDED_protocols_chromophore_ChromophoreDataReader_fwd_hh
