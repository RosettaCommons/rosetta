// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file FragmentSampler.cc
/// @brief ab-initio fragment assembly protocol for proteins
/// @detailed
///	  Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange

// Unit Headers
#include <protocols/noesy_assign/StructureIndependentPeakCalibrator.hh>

// Package Headers
// AUTO-REMOVED #include <protocols/noesy_assign/CrossPeakList.hh>

// Project Headers

// Utility headers
#include <basic/Tracer.hh>

//// C++ headers
// AUTO-REMOVED #include <cmath>

#include <utility/vector1.hh>


static thread_local basic::Tracer tr( "protocols.noesy_assign.crosspeaks" );

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
