// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/make_rot_lib/MakeRotLibJobInputterCreator.hh
/// @brief  MakeRotLibJobInputterCreator function declarations.
/// @author P. Douglas Renfrew ( renfrew@nyu.edu )

#ifndef INCLUDED_protocols_make_rot_lib_MakeRotLibJobInputterCreator_hh
#define INCLUDED_protocols_make_rot_lib_MakeRotLibJobInputterCreator_hh

#include <protocols/jd2/JobInputterCreator.hh>

namespace protocols {
namespace make_rot_lib {

class MakeRotLibJobInputterCreator : public protocols::jd2::JobInputterCreator {
public:
	virtual jd2::JobInputterOP create_JobInputter() const;
	virtual std::string keyname() const;

};

} //make_rot_lib
} //protocols

#endif //INCLUDED_protocols_make_rot_lib_MakeRotLibJobInputterCreator_hh
