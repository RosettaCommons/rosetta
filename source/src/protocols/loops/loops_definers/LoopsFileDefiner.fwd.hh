// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/loops_definers/LoopsFileDefiner.fwd.hh
/// @brief  LoopsFileDefiner class forward declarations header
/// @author Matthew O'Meara (mattjomeara@gmail.com)


#ifndef INCLUDED_protocols_loops_loops_definers_LoopsFileDefiner_FWD_HH
#define INCLUDED_protocols_loops_loops_definers_LoopsFileDefiner_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace loops {
namespace loops_definers {

// Forward
class LoopsFileDefiner;

typedef utility::pointer::shared_ptr< LoopsFileDefiner > LoopsFileDefinerOP;
typedef utility::pointer::shared_ptr< LoopsFileDefiner const > LoopsFileDefinerCOP;

} //namespace
} //namespace
} //namespace

#endif //include guard
