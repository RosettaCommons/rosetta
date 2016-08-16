// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/indexed_structure_store/filters/FragmentLookupFilterCreator.hh
/// @brief  FilterCreators for the FragmentLookupFilter
/// @author Alex Ford (fordas@uw.edu)

#ifndef INCLUDED_protocols_indexed_structure_store_filters_FragmentLookupFilterCreator_hh
#define INCLUDED_protocols_indexed_structure_store_filters_FragmentLookupFilterCreator_hh

// Package Headers
#include <protocols/filters/FilterCreator.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// c++ headers
#include <string>

namespace protocols {
namespace indexed_structure_store {
namespace filters {

class FragmentLookupFilterCreator : public protocols::filters::FilterCreator
{
public:
	virtual protocols::filters::FilterOP create_filter() const;
	virtual std::string keyname() const;
};


} //namespace filters
} //namespace indexed_structure_store
} //namespace protocols

#endif
