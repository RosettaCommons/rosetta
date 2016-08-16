// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_basic_datacache_HierarchicalDataMap_fwd_hh
#define INCLUDED_basic_datacache_HierarchicalDataMap_fwd_hh

#include <utility/pointer/access_ptr.hh>

namespace basic {
namespace datacache {

class HierarchicalDataMap;

typedef utility::pointer::shared_ptr< HierarchicalDataMap > HierarchicalDataMapOP;
typedef utility::pointer::shared_ptr< HierarchicalDataMap const > HierarchicalDataMapCOP;
typedef utility::pointer::weak_ptr< HierarchicalDataMap > HierarchicalDataMapAP;
typedef utility::pointer::weak_ptr< HierarchicalDataMap const > HierarchicalDataMapCAP;

}
}

#endif
