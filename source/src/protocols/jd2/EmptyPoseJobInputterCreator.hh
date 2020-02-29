// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/EmptyPoseJobInputterCreator.hh
/// @brief  creator for EmptyPoseJobInputter
/// @author Danny Farrell

#ifndef INCLUDED_protocols_jd2_EmptyPoseJobInputterCreator_hh
#define INCLUDED_protocols_jd2_EmptyPoseJobInputterCreator_hh

#include <protocols/jd2/JobInputterCreator.hh>

namespace protocols {
namespace jd2 {

class EmptyPoseJobInputterCreator : public protocols::jd2::JobInputterCreator {
public:
	protocols::jd2::JobInputterOP create_JobInputter() const override;
	std::string keyname() const override;

};

} // jd2
} // protocols

#endif //INCLUDED_protocols_jd2_EmptyPoseJobInputterCreator_hh

