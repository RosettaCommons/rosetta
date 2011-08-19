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


#include "benchmark.hh"
#include "init_util.hh"
#include "core/io/pdb/pose_io.hh"
#include "core/pose/Pose.hh"
#include "protocols/ligand_docking/LigandDockProtocol.hh"
#include "basic/options/option.hh"
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
//#include <platform/types.hh>
//#include <core/types.hh>
//#include <core/chemical/AA.hh>
//#include <core/chemical/ResidueType.fwd.hh>
//#include <core/chemical/ResidueTypeSet.fwd.hh>
//#include <core/conformation/Conformation.fwd.hh>
//#include <core/conformation/Residue.fwd.hh>
//#include <core/conformation/signals/XYZEvent.fwd.hh>
//#include <core/grid/CartGrid.fwd.hh>
//#include <core/io/atom_tree_diffs/atom_tree_diff.hh>
//#include <core/import_pose/file_data.fwd.hh>
//#include <core/io/pdb/file_data.hh>
//#include <core/scoring/Energies.fwd.hh>
//#include <core/scoring/EnergyMap.fwd.hh>
//#include <core/scoring/ScoreFunction.fwd.hh>
//#include <core/scoring/ScoreFunctionInfo.fwd.hh>
//#include <core/scoring/constraints/Constraint.fwd.hh>
//#include <core/scoring/constraints/ConstraintSet.fwd.hh>
//#include <core/scoring/constraints/Constraints.fwd.hh>
//#include <core/scoring/constraints/Constraints.hh>
//#include <core/scoring/constraints/Func.fwd.hh>
//#include <core/scoring/constraints/XYZ_Func.fwd.hh>
//#include <core/id/SequenceMapping.fwd.hh>
//#include <basic/MetricValue.fwd.hh>
//// AUTO-REMOVED #include <basic/OStream.fwd.hh>
//#include <utility/stream_util.hh>
//#include <basic/Tracer.fwd.hh>
//#include <basic/Tracer.hh>
//#include <basic/datacache/BasicDataCache.fwd.hh>
//#include <protocols/filters/Filter.fwd.hh>
//#include <protocols/match/downstream/ExternalGeomSampler.fwd.hh>
//#include <protocols/moves/DataMap.fwd.hh>
//#include <protocols/moves/MonteCarlo.fwd.hh>
//#include <protocols/moves/Mover.fwd.hh>
//#include <protocols/moves/Mover.hh>
//#include <protocols/moves/MoverStatistics.hh>
//#include <protocols/moves/MoverStatus.hh>
//#include <utility/Bound.fwd.hh>
//#include <utility/Bound.hh>
//#include <utility/down_cast.hh>
//#include <utility/exit.hh>

#include <protocols/jd2/JobDistributor.hh>

class LigandDockScriptBenchmark : public Benchmark
{
public:
	LigandDockScriptBenchmark(std::string name) : Benchmark(name) {};

	virtual void setUp() {
		basic::options::option.load_options_from_file("ligand_dock/ligand_dock_script_flags.txt");
	};

	virtual void run(int scaleFactor) {
		protocols::moves::MoverOP mover;

		for(int i=0; i<scaleFactor; i++) {
			protocols::jd2::JobDistributor::get_instance()->go(mover);
		}
	};

	virtual void tearDown() {};
};
