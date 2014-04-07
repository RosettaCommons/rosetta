// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/cmiles/helix_rotate.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Utility headers
#include <devel/init.hh>

// Project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/moves/HelixRotate.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/rational_mc/RationalMonteCarlo.hh>
#include <protocols/viewer/viewers.hh>

void* viewer_main(void* ) {
  using core::pose::Pose;
  using core::scoring::ScoreFunctionOP;
  using core::scoring::ScoreFunctionFactory;
  using protocols::loops::Loop;
  using protocols::moves::HelixRotate;
  using protocols::moves::MoverOP;
  using protocols::simple_moves::rational_mc::RationalMonteCarlo;

  Pose input  = *core::import_pose::pose_from_pdb("/work/tex/casp9_benchmark/meval/fast_cm/T0538/2kruA_1.pdb_full_length.pdb");
  Pose output = *core::import_pose::pose_from_pdb("/work/tex/casp9_benchmark/meval/fast_cm/T0538/2kruA_1.pdb_full_length.pdb");

  core::util::switch_to_residue_type_set(input, core::chemical::CENTROID);
  core::util::switch_to_residue_type_set(output, core::chemical::CENTROID);

  // Translate the specified sheet
  Loop helix(41, 52);
  double dist = 10.0;

  MoverOP base_mover = new HelixRotate(helix, dist);
  ScoreFunctionOP score = ScoreFunctionFactory::create_score_function("score0");

  RationalMonteCarlo mc(base_mover, score, 10000, 2.0, false);
  mc.apply(output);
}

int main(int argc, char* argv[]) {
	try {

		devel::init(argc, argv);
		protocols::viewer::viewer_main(viewer_main);
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
