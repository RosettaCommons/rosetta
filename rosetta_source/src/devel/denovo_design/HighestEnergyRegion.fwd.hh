// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/flxbb/filters/HighestEnergyRegion.fwd.hh
/// @brief Restrict design to residues in regions with poor total energy
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_flxbb_filters_highestenergyregion_fwd_hh
#define INCLUDED_protocols_flxbb_filters_highestenergyregion_fwd_hh

// utilty headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace flxbb {
namespace filters {

  // Forward declaration
  class HighestEnergyRegionOperation;
	class DesignByFragmentQualityOperation;

  // Types
  typedef utility::pointer::owning_ptr< HighestEnergyRegionOperation > HighestEnergyRegionOperationOP;
  typedef utility::pointer::owning_ptr< HighestEnergyRegionOperation const > HighestEnergyRegionOperationCOP;
  typedef utility::pointer::owning_ptr< DesignByFragmentQualityOperation > DesignByFragmentQualityOperationOP;
  typedef utility::pointer::owning_ptr< DesignByFragmentQualityOperation const > DesignByFragmentQualityOperationCOP;

}
}
}
#endif


