// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/remodel/RemodelEnzdesCstModule.hh
///
/// @brief this file handles merging constraint defined by enzdes type cstfile
/// @brief and blueprint definition of positions and add them to the pose
/// @author Possu Huang, possu@u.washington.edu, Jan 2010


#ifndef INCLUDED_protocols_forge_remodel_RemodelEnzdesCstModule_hh
#define INCLUDED_protocols_forge_remodel_RemodelEnzdesCstModule_hh

//#include <protocols/enzdes/enzdes/EnzdesRemodelProtocol.hh>
//#include <protocols/enzdes/EnzdesFlexBBProtocol.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/forge/remodel/RemodelData.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>


namespace protocols{
namespace forge {
namespace remodel {

class RemodelEnzdesCstModule : public protocols::toolbox::match_enzdes_util::EnzConstraintIO {

public:

	RemodelEnzdesCstModule(protocols::forge::remodel::RemodelData external_data);

	~RemodelEnzdesCstModule();

	void apply(core::pose::Pose & pose);
	void blueprint_cst_definition(core::pose::Pose & pose);
  void enable_constraint_scoreterms(core::scoring::ScoreFunctionOP scorefxn);
	void use_backbone_only_blocks();
	void use_all_blocks();
private: // data

	protocols::forge::remodel::RemodelData remodel_data_;
	core::Size cstblocksize_;
	core::scoring::ScoreFunctionOP scorefxn_;
	bool backbone_only_;

};

} // namespace remodel
} // namespace forge
} // namespace protocols

#endif
