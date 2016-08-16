// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author John Karanicolas


#ifndef INCLUDED_basic_MetricValue_fwd_hh
#define INCLUDED_basic_MetricValue_fwd_hh

#include <utility/pointer/owning_ptr.hh>

#include <utility/down_cast.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <cassert>
#include <iosfwd>


namespace basic {

class MetricValueBase;

typedef utility::pointer::shared_ptr< MetricValueBase > MetricValueBaseOP;
typedef utility::pointer::shared_ptr< MetricValueBase const > MetricValueBaseCOP;

} // namespace basic


#endif
