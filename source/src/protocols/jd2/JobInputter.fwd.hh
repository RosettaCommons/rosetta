// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/JobInputter.fwd.hh
/// @brief  header file for JobInputter class, part of August 2008 job distributor as planned at RosettaCon08
/// @author Steven Lewis smlewi@gmail.com

#ifndef INCLUDED_protocols_jd2_JobInputter_fwd_hh
#define INCLUDED_protocols_jd2_JobInputter_fwd_hh

#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace jd2 {

class JobInputter;
typedef utility::pointer::access_ptr< JobInputter > JobInputterAP;
typedef utility::pointer::access_ptr< JobInputter const > JobInputterCAP;
typedef utility::pointer::owning_ptr< JobInputter > JobInputterOP;
typedef utility::pointer::owning_ptr< JobInputter const > JobInputterCOP;

// wrap enum in a descriptive namespace to reduce conflicts and keep down
// namespace pollution
namespace JobInputterInputSource {

/// @enum describes the type of input source that the JobInputter is using
/// Also please add the lookup in job_inputter_input_source_to_string()
enum Enum {
	NONE,
	UNKNOWN,
	POSE,
	SILENT_FILE,
	PDB_FILE,
	ATOM_TREE_FILE,
	DATABASE,
	MAKE_ROT_LIB,
	RESOURCE_MANAGED_JOB,
	SCREENING_FILE
};

}

}//jd2
}//protocols

#endif //INCLUDED_protocols_jd2_JobInputter_FWD_HH
