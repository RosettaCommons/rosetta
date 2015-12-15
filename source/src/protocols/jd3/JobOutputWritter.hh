// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/JobOutputWriter.hh
/// @brief  The declaration for class protocols::jd3::JobOutputWriter
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_JobOutputWriter_hh
#define INCLUDED_protocols_jd3_JobOutputWriter_hh

// Unit headers
#include <protocols/jd3/JobOutputWriter.fwd.hh>

// Package headers
#include <protocols/jd3/InnerJobOutputWriter.fwd.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <core/types.hh>

//C++ headers
#include <string>
#include <list>
#include <map>
#include <set>

#include <utility/vector1.hh>


namespace protocols {
namespace jd3 {

/// @brief %JobOutputWriter class is responsible for taking data out of the JobOutput
/// object and writing it to some sort of output destination (a file, or a database).
class JobOutputWriter : public utility::pointer::ReferenceCount
{
public:

	JobOutputWriter();
	virtual ~JobOutputWriter();

	virtual
	void write_output_for_job( LarvalJobOP, JobOutputOP ) = 0;

}; // JobOutputWriter

} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_JobOutputWriter_HH
