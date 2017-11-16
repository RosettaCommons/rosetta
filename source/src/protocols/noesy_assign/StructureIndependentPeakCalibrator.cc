// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file FragmentSampler.cc
/// @brief ab-initio fragment assembly protocol for proteins
/// @details
///   Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange

// Unit Headers
#include <protocols/noesy_assign/StructureIndependentPeakCalibrator.hh>

// Package Headers

// Project Headers

// Utility headers
#include <basic/Tracer.hh>

//// C++ headers

#include <utility/vector1.hh>


static basic::Tracer tr( "protocols.noesy_assign.crosspeaks" );

using core::Real;
using namespace core;
using namespace basic;
//using namespace basic::options;
//using namespace basic::options::OptionKeys;

namespace protocols {
namespace noesy_assign {

void StructureIndependentPeakCalibrator::collect_upperbound_statistics( core::Size peak, TypeCumulator const& types ) {
	collect_target_statistics( peaks()[ peak ]->distance_bound(), types);
}

}
}
