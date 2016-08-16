// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/protein/checker/ProteinAtrRepChecker.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_protein_checker_ProteinAtrRepChecker_HH
#define INCLUDED_protocols_stepwise_protein_checker_ProteinAtrRepChecker_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/stepwise/modeler/protein/checker/ProteinAtrRepChecker.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.hh>

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace protein {
namespace checker {

class ProteinAtrRepChecker: public utility::pointer::ReferenceCount {

public:

	//Constructor
	ProteinAtrRepChecker( pose::Pose const & pose,
		utility::vector1< Size > const & moving_res_list );

	~ProteinAtrRepChecker();

	Real delta_atr_score() const{ return delta_atr_score_; }
	Real delta_rep_score() const{ return delta_rep_score_; }
	Real base_atr_score() const{ return base_atr_score_; }
	Real base_rep_score() const{ return base_rep_score_; }

public:

	bool
	check_screen( pose::Pose & current_pose_screen );

private:

	void
	get_base_atr_rep_score( core::pose::Pose const & pose );

	void
	initialize_scorefxn();

private:

	utility::vector1< Size > const moving_res_list_;

	Real rep_cutoff_, atr_cutoff_;
	Real base_atr_score_;
	Real base_rep_score_;
	Real delta_atr_score_;
	Real delta_rep_score_;

	core::scoring::ScoreFunctionOP atr_rep_screening_scorefxn_;

};

} //checker
} //protein
} //modeler
} //stepwise
} //protocols

#endif
