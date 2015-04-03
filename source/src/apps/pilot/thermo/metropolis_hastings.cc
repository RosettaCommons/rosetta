// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file metropolis_hastings.cc
/// @brief Metropolis-Hastings Monte Carlo app
/// @author
/// @details
// devel headers
#include <devel/init.hh>
// protocols headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>
#include <protocols/moves/MonteCarlo.hh>
// core headers
#include <basic/options/option_macros.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>
//#include <basic/Tracer.hh>
// numeric headers
//#include <numeric/random/random.hh>
// utility headers
#include <utility/file/file_sys_util.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
//basic::Tracer TR("apps.metropolis_hastings");
OPT_1GRP_KEY(Integer, mh, ntrials)
OPT_1GRP_KEY(Real, mh, kt)
OPT_1GRP_KEY(Real, mh, br_weight)
OPT_1GRP_KEY(Real, mh, sm_weight)
OPT_1GRP_KEY(Real, mh, sh_weight)
OPT_1GRP_KEY(Real, mh, sc_weight)
OPT_1GRP_KEY(Real, mh, sc_mc_weight)
OPT_1GRP_KEY(Integer, mh, sc_mc_ntrials)
OPT_1GRP_KEY(Real, mh, sc_prob_uniform)
OPT_1GRP_KEY(Real, mh, sc_prob_withinrot)
int
main( int argc, char * argv [] )
{
try {
	NEW_OPT(mh::ntrials, "number of total Monte Carlo trials", 1000);
	NEW_OPT(mh::kt, "value of kT for Monte Carlo", 0.6);
	NEW_OPT(mh::br_weight, "proportion of backrub moves", 0.4);
	NEW_OPT(mh::sm_weight, "proportion of small moves", 0.1);
	NEW_OPT(mh::sh_weight, "proportion of shear moves", 0.1);
	NEW_OPT(mh::sc_weight, "proportion of side chain moves", 0.4);
	NEW_OPT(mh::sc_mc_weight, "proportion of side chain Monte Carlo moves", 0.0);
	NEW_OPT(mh::sc_mc_ntrials, "number of side chain Monte Carlo trials", 100);
	NEW_OPT(mh::sc_prob_uniform, "probability of uniformly sampling chi angles", 0.1);
	NEW_OPT(mh::sc_prob_withinrot, "probability of sampling within the current rotamer", 0.2);
	devel::init(argc, argv);
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	protocols::canonical_sampling::MetropolisHastingsMoverOP metropolis_hastings_mover(new protocols::canonical_sampling::MetropolisHastingsMover);
	core::scoring::ScoreFunctionOP score_fxn(core::scoring::get_score_function());
	protocols::moves::MonteCarloOP monte_carlo(new protocols::moves::MonteCarlo(*score_fxn, option[ mh::kt ]));
	metropolis_hastings_mover->set_monte_carlo(monte_carlo);
	metropolis_hastings_mover->set_ntrials(option[ mh::ntrials ]);
	bool preserve_cbeta = option[ mh::br_weight ] ? true : false;
	metropolis_hastings_mover->add_backrub_mover(option[ mh::br_weight ]);
	metropolis_hastings_mover->add_small_mover(option[ mh::sm_weight ]);
	metropolis_hastings_mover->add_shear_mover(option[ mh::sh_weight ]);
	metropolis_hastings_mover->add_sidechain_mover(option[ mh::sc_weight ], option[ mh::sc_prob_uniform ], option[ mh::sc_prob_withinrot ], preserve_cbeta);
	metropolis_hastings_mover->add_sidechain_mc_mover(option[ mh::sc_mc_weight ], option[ mh::sc_prob_uniform ], option[ mh::sc_prob_withinrot ], preserve_cbeta, option[ mh::sc_mc_ntrials ]);
	protocols::jd2::JobDistributor::get_instance()->go(metropolis_hastings_mover);
} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl; }
		return -1;
}
