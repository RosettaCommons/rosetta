// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <basic/datacache/HierarchicalDataMap.hh>

namespace basic {
namespace datacache {

HierarchicalDataMap::HierarchicalDataMap() {}

HierarchicalDataMap::~HierarchicalDataMap() = default;

void HierarchicalDataMap::set_parent(HierarchicalDataMapCAP parent) {
	parent_ = parent;
}

void HierarchicalDataMap::unset_parent() {
	parent_ = HierarchicalDataMapCAP();
}

}
}
