// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
/// @file   FlexPepDockingFlags.fwd.hh
///
/// @brief flags structure for FlexPepDocking protocols - forward declarations
/// @date January 1, 2009
/// @author Barak Raveh

#ifndef INCLUDED_protocols_flexpep_docking_FlexPepDockingFlags_fwd_hh
#define INCLUDED_protocols_flexpep_docking_FlexPepDockingFlags_fwd_hh

#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace flexpep_docking {

class FlexPepDockingFlags;
typedef utility::pointer::weak_ptr< FlexPepDockingFlags > FlexPepDockingFlagsAP;
typedef utility::pointer::weak_ptr< FlexPepDockingFlags const > FlexPepDockingFlagsCAP;
typedef utility::pointer::shared_ptr< FlexPepDockingFlags > FlexPepDockingFlagsOP;
typedef utility::pointer::shared_ptr< FlexPepDockingFlags const > FlexPepDockingFlagsCOP;

} // namespace flexpep_docking
} // namespace protocols

#endif

