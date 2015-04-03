// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/flxbb/FlxbbDesign_main.cc
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit Headers
#include <protocols/flxbb/FlxbbDesign_main.hh>
#include <protocols/flxbb/FlxbbDesign.hh>

// Project headers
#include <protocols/moves/Mover.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/constraints/util.hh>
#include <basic/options/option.hh>
#include <protocols/jd2/JobDistributor.hh>
// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


using namespace basic::options;
using namespace basic::options::OptionKeys;

using namespace core;

namespace protocols {
namespace flxbb{

	void FlxbbDesign_main(){

		using namespace protocols::jobdist;
		using namespace protocols::moves;
		using protocols::jd2::JobDistributor;

		scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
		scoring::ScoreFunctionOP scorefxn_design = scorefxn;
		scoring::ScoreFunctionOP scorefxn_relax = scorefxn;

		if( option[ in::file::fullatom ]() ) {
			scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn( *scorefxn_relax );
		} else{
			scoring::constraints::add_constraints_from_cmdline_to_scorefxn( *scorefxn_relax );
		}

		MoverOP protocol;
		protocol = MoverOP( new FlxbbDesign( scorefxn_design, scorefxn_relax ) );
		JobDistributor::get_instance()->go( protocol );

	}

} // namespace flxbb
} // namespace protocols
