// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   rosetta/benchmark/LigandDock.bench.hh
///
/// @brief Dock the ligand in the 7cpa complex.
/// Use all options (flexible ligand, flexible backbone)
/// @author Gordon Lemmon


#include <apps/benchmark/performance/performance_benchmark.hh>
#include <basic/options/option.hh>

//Auto Headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <utility/vector1.hh>

class LigandDockScriptBenchmark : public PerformanceBenchmark
{
public:
	LigandDockScriptBenchmark(std::string name) : PerformanceBenchmark(name) {};

	virtual void setUp() {
		basic::options::option.load_options_from_file("ligand_dock/ligand_dock_script_flags.txt");
	};

	virtual void run(core::Real scaleFactor) {
		protocols::moves::MoverOP mover;

		for(int i=0; i<scaleFactor; i++) {
			protocols::jd2::JobDistributor::get_instance()->go(mover);
		}
	};

	virtual void tearDown() {};
};
