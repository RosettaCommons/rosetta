// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/hotspot_hashing/SearchPatternRotSetOp.hh
/// @brief  Creates rigid body variants from search pattern during repacking.
/// @author Alex Ford (fordas@uw.edu)

#ifndef INCLUDED_protocols_hotspot_hashing_SearchPatternRotSetOp_fwd_hh
#define INCLUDED_protocols_hotspot_hashing_SearchPatternRotSetOp_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace hotspot_hashing {

class SearchPatternRotSetOp;
typedef utility::pointer::shared_ptr< SearchPatternRotSetOp > SearchPatternRotSetOpOP;
typedef utility::pointer::shared_ptr< SearchPatternRotSetOp const > SearchPatternRotSetOpCOP;

class AddSearchPatternRotSetOp;
typedef utility::pointer::shared_ptr< AddSearchPatternRotSetOp > AddSearchPatternRotSetOpOP;
typedef utility::pointer::shared_ptr< AddSearchPatternRotSetOp const > AddSearchPatternRotSetOpCOP;

}
}

#endif
