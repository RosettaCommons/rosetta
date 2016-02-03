// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/cmiles/sheet_translate.cc
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
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/rational_mc/RationalMonteCarlo.hh>
#include <protocols/moves/SheetTranslate.hh>
#include <protocols/viewer/viewers.hh>

void run(protocols::moves::MoverOP base_mover, core::pose::Pose* pose) {
  using core::scoring::ScoreFunctionFactory;
  using protocols::simple_moves::rational_mc::RationalMonteCarlo;
  assert(pose);

  RationalMonteCarlo mc(
      base_mover,
      ScoreFunctionFactory::create_score_function("score0"),
      800,
      10.0,
      false);

  mc.apply(*pose);
}

void* viewer_main(void* ) {
  using core::pose::Pose;
  using protocols::loops::Loop;
  using protocols::moves::MoverOP;
  using protocols::moves::SheetTranslate;

  Pose input  = *core::import_pose::pose_from_file("/work/tex/casp9_benchmark/meval/fast_cm/T0552/2oxgZ_1.pdb_full_length.pdb", core::import_pose::PDB_file);
  Pose output = *core::import_pose::pose_from_file("/work/tex/casp9_benchmark/meval/fast_cm/T0552/2oxgZ_1.pdb_full_length.pdb", core::import_pose::PDB_file);

  core::util::switch_to_residue_type_set(input, core::chemical::CENTROID);
  core::util::switch_to_residue_type_set(output, core::chemical::CENTROID);

  // Translate the specified sheet
  Loop sheet(51, 53);
  double dist = 0.01;

  MoverOP fw_mover = new SheetTranslate(sheet, +dist);
  MoverOP bw_mover = new SheetTranslate(sheet, -dist);

  for (unsigned i = 1; i <= 10; ++i) {
    run(fw_mover, &output);
    run(bw_mover, &output);
  }
}

int main(int argc, char* argv[]) {
    try{
  devel::init(argc, argv);
  protocols::viewer::viewer_main(viewer_main);
    } catch ( utility::excn::EXCN_Base const & e ) {
                             std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                }
       return 0;
}
