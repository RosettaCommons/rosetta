// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/LoopsFileLoader.fwd.hh
/// @brief
/// @author

#ifndef INCLUDED_protocols_loops_LoopsFileLoader_FWD_HH
#define INCLUDED_protocols_loops_LoopsFileLoader_FWD_HH

//utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace loops {

class LoopsFileLoader;
typedef utility::pointer::shared_ptr< LoopsFileLoader > LoopsFileLoaderOP;
typedef utility::pointer::shared_ptr< LoopsFileLoader const > LoopsFileLoaderCOP;


} // namespace loops
} // namespace protocols

#endif //INCLUDED_protocols_resource_manager_ResourceLoader_FWD_HH
