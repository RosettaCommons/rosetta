// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/simple_filter/SidechainDepthFilter.hh
/// @brief Filter for measuring minimum distance from any sidechain atom in a reside
//  to the closest water molecule
/// @author Hahnbeom Park

#ifndef INCLUDED_protocols_simple_filters_SidechainDepthFilter_fwd_hh
#define INCLUDED_protocols_simple_filters_SidechainDepthFilter_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace simple_filters {

class ResidueDepthFilter;
struct ResidueDepthData;

typedef utility::pointer::shared_ptr< ResidueDepthData > ResidueDepthDataOP;
typedef utility::pointer::shared_ptr< ResidueDepthData const > ResidueDepthDataCOP;


}
}

#endif
