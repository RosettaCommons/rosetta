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
/// @author Colin A. Smith


#ifndef INCLUDED_basic_MetricValueIO_hh
#define INCLUDED_basic_MetricValueIO_hh

#include <basic/MetricValue.fwd.hh>

#include <iosfwd>


namespace basic {

/// @brief check whether a MetricValue can be written to or read from a stream
bool
handles_metric_value(
	MetricValueBase const & metric_value
);

/// @brief write a MetricValue to a stream, returns true if successful
/// @details
/// This function currently supports some of the most often used MetricValues. The
/// type is always encoded in the output so that an object of the same class can be
/// generated during reading. The supported types and their formats is as follows:
///
/// double                    double <num>
/// int                           Int <num>
/// std::size_t                    std::size_t <num>
/// bool                          Bool <0 or 1>
/// utility::vector1<double>  Real[ <num1> <num2> ... ]
/// utility::vector1<int>         Int[ <num1> <num2> ... ]
/// utility::vector1<size_t>  size_t[ <num1> <num2> ... ]
/// utility::vector1<bool>        Bool[ <bool1> <bool2> ... ]
bool
write_metric_value(
	std::ostream & os,
	MetricValueBase const & metric_value
);

/// @brief read a MetricValue from a stream, returns true if successful
bool
read_metric_value(
	std::istream & is,
	MetricValueBase & metric_value
);

/// @brief read a MetricValue from a stream, returns NULL if unsuccessful
MetricValueBaseOP
read_metric_value(
	std::istream & is
);


} // namespace basic

#endif
