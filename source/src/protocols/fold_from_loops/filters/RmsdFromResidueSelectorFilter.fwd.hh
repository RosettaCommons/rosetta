// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   RmsdFromResidueSelectorFilter.fwd.hh
/// @brief  Evaluate RMSD between two poses allowing to select the regions to compare in each pose through ResidueSelector
/// @author Jaume Bonet (jaume.bonet@gmail.com)

#ifndef INCLUDED_protocols_fold_from_loops_filters_RmsdFromResidueSelectorFilter_fwd_hh
#define INCLUDED_protocols_fold_from_loops_filters_RmsdFromResidueSelectorFilter_fwd_hh
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace fold_from_loops {
namespace filters {

class RmsdFromResidueSelectorFilter;

typedef utility::pointer::shared_ptr< RmsdFromResidueSelectorFilter > RmsdFromResidueSelectorFilterOP;
typedef utility::pointer::shared_ptr< RmsdFromResidueSelectorFilter const > RmsdFromResidueSelectorFilterCOP;

}
}
}

#endif
