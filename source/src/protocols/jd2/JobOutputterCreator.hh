// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/JobOutputterCreator.hh
/// @brief  Base class for JobOutputterCreators for the JobOutputter load-time factory registration scheme
/// @author Steven Lewis smlewi@gmail.com, Brian Weitzner brian.weitzner@gmail.com

#ifndef INCLUDED_protocols_jd2_JobOutputterCreator_hh
#define INCLUDED_protocols_jd2_JobOutputterCreator_hh

// Unit Headers
#include <protocols/jd2/JobOutputter.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// c++ headers
#include <string>

namespace protocols {
namespace jd2 {

/// @brief Abstract base class for a JobOutputter factory; the Creator class is responsible for
/// creating a particular mover class.
class JobOutputterCreator : public utility::pointer::ReferenceCount
{
public:
	JobOutputterCreator();
	virtual ~JobOutputterCreator();

	virtual JobOutputterOP create_JobOutputter() const = 0;
	virtual std::string keyname() const = 0;
};

typedef utility::pointer::shared_ptr< JobOutputterCreator > JobOutputterCreatorOP;
typedef utility::pointer::shared_ptr< JobOutputterCreator const > JobOutputterCreatorCOP;

} //namespace jd2
} //namespace protocols

#endif
