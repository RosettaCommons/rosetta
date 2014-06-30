// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ./src/protocols/indexed_structure_store/filters/FragmentLookupFilter.hh
/// @brief Header for FragmentLookupFilter class.
/// @detailed
/// @author Alex Ford (fordas@uw.edu)


#ifndef INCLUDED_protocols_indexed_structure_store_filters_FragmentLookupFilter_fwd_hh
#define INCLUDED_protocols_indexed_structure_store_filters_FragmentLookupFilter_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace protocols {
namespace indexed_structure_store {
namespace filters {

class FragmentLookupFilter;

typedef utility::pointer::owning_ptr<FragmentLookupFilter> FragmentLookupFilterOP;
typedef utility::pointer::owning_ptr<FragmentLookupFilter const> FragmentLookupFilterCOP;
typedef utility::pointer::access_ptr<FragmentLookupFilter> FragmentLookupFilterAP;
typedef utility::pointer::access_ptr<FragmentLookupFilter const> FragmentLookupFilterCAP;

}
}
}

#endif
