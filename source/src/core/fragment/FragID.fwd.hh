// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragments/FragID.fwd.hh
/// @brief  set of fragments for a certain alignment frame
/// @author Oliver Lange (olange@u.washington.edu)
/// @author James Thompson
/// @date   Wed Oct 20 12:08:31 2007
///
#ifndef INCLUDED_core_fragment_FragID_FWD_HH
#define INCLUDED_core_fragment_FragID_FWD_HH

#include <utility/vector1.fwd.hh>

namespace core {
namespace fragment {

class FragID;
typedef utility::vector1< FragID > FragID_List;

} //fragment
} //core

#endif
