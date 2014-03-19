// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

///@file   protocols/pose_selectors/BasicPoseSelectorCreators.hh
///@brief  Collection of PoseSelector creators for basic selectors
///@author Luki Goldschmidt <luki@uw.edu>

#ifndef INCLUDED_protocols_pose_selectors_BasicPoseSelectorCreators_hh
#define INCLUDED_protocols_pose_selectors_BasicPoseSelectorCreators_hh

#include <protocols/rosetta_scripts/PoseSelectorCreator.hh>
#include <string>

namespace protocols {
namespace pose_selectors {

class AndSelectorCreator : public protocols::rosetta_scripts::PoseSelectorCreator {
public:
	virtual protocols::rosetta_scripts::PoseSelectorOP create_selector() const;
	virtual std::string keyname() const { return "AndSelector"; }
};

class OrSelectorCreator : public protocols::rosetta_scripts::PoseSelectorCreator {
public:
	virtual protocols::rosetta_scripts::PoseSelectorOP create_selector() const;
	virtual std::string keyname() const { return "OrSelector"; }
};

class TopNByPropertyCreator : public protocols::rosetta_scripts::PoseSelectorCreator {
public:
	virtual protocols::rosetta_scripts::PoseSelectorOP create_selector() const;
	virtual std::string keyname() const { return "TopNByProperty"; }
};

}
}

#endif

