// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/simple_metrics/SimpleMetric.cc
/// @brief The base class for Metrics in the Metric/Filter/Reporter system
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


// Unit header or inline function header
#include <core/simple_metrics/SimpleMetric.hh>

// NOTE: This file should have NO dependencies other than its header.


namespace core {
namespace simple_metrics {



SimpleMetric::SimpleMetric( std::string const & type_name ):
	utility::pointer::ReferenceCount(),
	type_( type_name )

{

}

SimpleMetric::~SimpleMetric(){}

SimpleMetric::SimpleMetric( SimpleMetric const & src ):
	utility::pointer::ReferenceCount(),
	type_(src.type_)
{

}



} //core
} //metrics


