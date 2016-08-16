// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/JobDistributorFactory
/// @brief  JobDistributorFactory class, part of August 2008 job distributor as planned at RosettaCon08.
/// @author Andrew Leaver-Fay
/// @author Steven Lewis smlewi@gmail.com

#ifndef INCLUDED_protocols_jd2_JobDistributorFactory_hh
#define INCLUDED_protocols_jd2_JobDistributorFactory_hh

// Unit headers
#include <protocols/jd2/JobDistributorFactory.fwd.hh>

// Package headers
#include <protocols/jd2/JobDistributor.fwd.hh>
#include <protocols/jd2/JobInputter.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Parser.fwd.hh>

namespace protocols {
namespace jd2 {

class JobDistributorFactory {
public:
	static
	JobDistributor *
	create_job_distributor();

	static
	JobInputterOP
	create_job_inputter();

	static
	JobOutputterOP
	create_job_outputter();

	static
	JobOutputterOP
	create_job_outputter( JobOutputterOP default_jobout );

	static
	ParserOP
	create_parser();

};

} // namespace jd2
} // namespace protocols

#endif //INCLUDED_protocols_jd2_JobDistributorFactory_HH
