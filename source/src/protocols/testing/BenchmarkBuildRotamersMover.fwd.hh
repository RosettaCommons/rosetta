// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/testing/BenchmarkBuildRotamersMover.fwd.hh
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_protocols_testing_BenchmarkBuildRotamersMover_fwd_hh
#define INCLUDED_protocols_testing_BenchmarkBuildRotamersMover_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace testing {

class BenchmarkBuildRotamersMover;
using BenchmarkBuildRotamersMoverOP = utility::pointer::shared_ptr< BenchmarkBuildRotamersMover >;
using BenchmarkBuildRotamersMoverCOP = utility::pointer::shared_ptr< BenchmarkBuildRotamersMover const >;

} // moves
} // protocols

#endif //INCLUDED_protocols_testing_BenchmarkBuildRotamersMover_fwd
