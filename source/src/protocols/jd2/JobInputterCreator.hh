// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/JobInputterCreator.hh
/// @brief  Base class for JobInputterCreators for the JobInputter load-time factory registration scheme
/// @author Steven Lewis smlewi@gmail.com, Brian Weitzner brian.weitzner@gmail.com

#ifndef INCLUDED_protocols_jd2_JobInputterCreator_hh
#define INCLUDED_protocols_jd2_JobInputterCreator_hh

// Unit Headers
#include <protocols/jd2/JobInputter.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// c++ headers
#include <string>

namespace protocols {
namespace jd2 {

/// @brief Abstract base class for a JobInputter factory; the Creator class is responsible for
/// creating a particular mover class.
class JobInputterCreator : public utility::pointer::ReferenceCount
{
public:
	JobInputterCreator();
	~JobInputterCreator() override;

	virtual JobInputterOP create_JobInputter() const = 0;
	virtual std::string keyname() const = 0;
};

typedef utility::pointer::shared_ptr< JobInputterCreator > JobInputterCreatorOP;
typedef utility::pointer::shared_ptr< JobInputterCreator const > JobInputterCreatorCOP;

} //namespace jd2
} //namespace protocols

#endif
