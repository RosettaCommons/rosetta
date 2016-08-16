// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/PoseCommentsFeatures.cc
/// @brief  report feature data to database
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/Report.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <protocols/features/FeaturesReporter.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace features {

using utility::pointer::ReferenceCount;

Report::Report() :
	reporters_()
{}

Report::Report(Report const &) :
	ReferenceCount(),
	reporters_()
{}

Report::~Report() {}


} //namesapce
} //namespace
