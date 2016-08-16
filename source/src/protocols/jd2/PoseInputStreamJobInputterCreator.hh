// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/PoseInputStreamJobInputterCreator.hh
/// @brief
/// @author Steven Lewis smlewi@gmail.com, Brian Weitzner brian.weitzner@gmail.com

#ifndef INCLUDED_protocols_jd2_PoseInputStreamJobInputterCreator_hh
#define INCLUDED_protocols_jd2_PoseInputStreamJobInputterCreator_hh

#include <protocols/jd2/JobInputterCreator.hh>

namespace protocols {
namespace jd2 {

class PoseInputStreamJobInputterCreator : public protocols::jd2::JobInputterCreator {
public:
	virtual JobInputterOP create_JobInputter() const;
	virtual std::string keyname() const;

};

} //jd2
} //protocols

#endif //INCLUDED_protocols_jd2_PoseInputStreamJobInputterCreator_hh

