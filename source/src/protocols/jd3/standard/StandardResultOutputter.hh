// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/standard/StandardResultOutputter.hh
/// @brief  Definition of the %StandardResultOutputter class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_standard_StandardResultOutputter_HH
#define INCLUDED_protocols_jd3_standard_StandardResultOutputter_HH

// Unit headers
#include <protocols/jd3/standard/StandardResultOutputter.fwd.hh>

// Package headers
#include <protocols/jd3/output/ResultOutputter.hh>
#include <protocols/jd3/output/OutputSpecification.fwd.hh>
#include <protocols/jd3/JobResult.fwd.hh>

namespace protocols {
namespace jd3 {
namespace output {

/// @brief The %StandardResultOutputter
class StandardResultOutputter : output::ResultOutputter
{
public:

	StandardResultOutputter();
	virtual ~StandardResultOutputter();

	void write_output(
		OutputSpecification const & specification,
		JobResult const & result
	) override;

	void flush() override;

	void set_primary_outputter( pose_outputter::PoseOutputterOP outputter );
	void append_secondary_outputter( pose_outputter::SecondaryPoseOutputterOP outputter );

private:
	pose_outputter::PoseOutputterOP primary_outputter_;
	utility::vector1< pose_outputter::SecondaryPoseOutputterOP > secondary_outputters_;

};

} // namespace output
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_output_StandardResultOutputter_HH
