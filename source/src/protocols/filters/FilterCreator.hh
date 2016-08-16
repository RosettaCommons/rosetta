// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/filters/FilterCreator.hh
/// @brief  Base class for FilterCreators for the Filter load-time factory registration scheme
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_filters_FilterCreator_hh
#define INCLUDED_protocols_filters_FilterCreator_hh

// Unit Headers
#include <protocols/filters/Filter.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// c++ headers
#include <string>

namespace protocols {
namespace filters {

/// @brief Abstract base class for a Filter factory; the Creator class is responsible for
/// creating a particular filter class.
class FilterCreator : public utility::pointer::ReferenceCount
{
public:
	FilterCreator();
	virtual ~FilterCreator();

	virtual FilterOP create_filter() const = 0;
	virtual std::string keyname() const = 0;
};

typedef utility::pointer::shared_ptr< FilterCreator > FilterCreatorOP;
typedef utility::pointer::shared_ptr< FilterCreator const > FilterCreatorCOP;

} //namespace filters
} //namespace protocols

#endif
