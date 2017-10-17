// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/output/ResultOutputter.hh
/// @brief  Definition of the %ResultOutputter class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_output_ResultOutputter_HH
#define INCLUDED_protocols_jd3_output_ResultOutputter_HH

// Unit headers
#include <protocols/jd3/output/ResultOutputter.fwd.hh>

// Package headers
#include <protocols/jd3/output/OutputSpecification.fwd.hh>
#include <protocols/jd3/JobResult.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace jd3 {
namespace output {

/// @brief The %ResultOutputter
class ResultOutputter : public utility::pointer::ReferenceCount
{
public:

	ResultOutputter();
	virtual ~ResultOutputter();

	virtual
	void write_output(
		OutputSpecification const & specification,
		JobResult const & result
	) = 0;

	/// @brief Output from an outputter may be held back and only flushed when requested
	/// by the JobQueen; I/O can be expensive, so it's a good idea to gather up the
	/// results of many outputs before flushing them to disk.
	virtual
	void flush() = 0;

};

} // namespace output
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_output_ResultOutputter_HH
