// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/filters/FilterFactory.fwd.hh
/// @brief  Forward declaration of Factory class for Filters
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_filters_FilterFactory_fwd_hh
#define INCLUDED_protocols_filters_FilterFactory_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace filters {

class FilterFactory;

typedef utility::pointer::shared_ptr< FilterFactory > FilterFactoryOP;
typedef utility::pointer::shared_ptr< FilterFactory const > FilterFactoryCOP;

} //namespace filters
} //namespace protocols

#endif
