// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   rosetta/benchmark/ScoreEach.bench.cc
///
/// @brief  Scoring Analytic Etable benchmark
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#include <apps/benchmark/performance/ScoreEach.bench.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

class ScoreAnalyticEtableBenchmark : public ScoreEachBenchmark
{
public:

	ScoreAnalyticEtableBenchmark(
		std::string name,
		core::scoring::ScoreType score_type,
		core::Size base_scale_factor
	) :
		ScoreEachBenchmark(name, score_type, base_scale_factor)
	{}

	virtual void setUp() {
		using namespace core::scoring;
		using namespace core::scoring::methods;

		ScoreFunctionOP scorefxn( new ScoreFunction() );
		EnergyMethodOptions energy_method_options;
		energy_method_options.analytic_etable_evaluation(true);
		scorefxn->set_energy_method_options(energy_method_options);
		set_scorefxn(scorefxn);

		ScoreEachBenchmark::setUp();
	}
};

ScoreAnalyticEtableBenchmark Score_analytic_etable_fa_atr_("core.scoring.Score_analytic_etable_100x_fa_atr",fa_atr,100);
ScoreAnalyticEtableBenchmark Score_analytic_etable_fa_rep_("core.scoring.Score_analytic_etable_100x_fa_rep",fa_rep,100);
ScoreAnalyticEtableBenchmark Score_analytic_etable_fa_sol_("core.scoring.Score_analytic_etable_100x_fa_sol",fa_sol,100);
ScoreAnalyticEtableBenchmark Score_analytic_etable_fa_intra_atr_("core.scoring.Score_analytic_etable_100x_fa_intra_atr",fa_intra_atr,100);
ScoreAnalyticEtableBenchmark Score_analytic_etable_fa_intra_rep_("core.scoring.Score_analytic_etable_100x_fa_intra_rep",fa_intra_rep,100);
ScoreAnalyticEtableBenchmark Score_analytic_etable_fa_intra_sol_("core.scoring.Score_analytic_etable_100x_fa_intra_sol",fa_intra_sol,100);
