// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/output/MultipleOutputter.cc
/// @brief  Definition of the %MultipleOutputter class's methods
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/jd3/output/MultipleOutputter.hh>

// Package headers
#include <protocols/jd3/output/MultipleOutputSpecification.hh>

namespace protocols {
namespace jd3 {
namespace output {

MultipleOutputter::MultipleOutputter() {}
MultipleOutputter::~MultipleOutputter() {}

/// @brief Invoke write_output on all of the ResultOutputters this %MultipleOutputter contains.
/// This class expects the OutputSpecification to be of type MultipleOutputSpecification and will
/// hand the ResultSpecifications that this class contains to the corresponding ResultOutputter.
void MultipleOutputter::write_output(
	OutputSpecification const & spec,
	JobResult const & result
)
{
	typedef MultipleOutputSpecification MOS;
	debug_assert( dynamic_cast< MOS const * > (&spec) );
	MOS const & mos_spec( static_cast< MOS const & > (spec) );
	debug_assert( mos_spec.output_specifications().size() == outputters_.size() );
	for ( core::Size ii = 1; ii <= mos_spec.output_specifications().size(); ++ii ) {
		outputters_[ ii ]->write_output( *mos_spec.output_specifications()[ ii ], result );
	}
}

/// @brief Invoke flush on all of the ResultOutputters this %MultipleOutputter contains
void MultipleOutputter::flush()
{
	for ( auto outputter : outputters_ ) {
		outputter->flush();
	}
}

void MultipleOutputter::append_outputter( ResultOutputterOP outputter )
{
	outputters_.push_back( outputter );
}

utility::vector1< ResultOutputterOP > const & MultipleOutputter::outputters() const
{
	return outputters_;
}


} // namespace output
} // namespace jd3
} // namespace protocols
