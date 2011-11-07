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

// C/C++ headers
#include <iostream>

// Utility headers
#include <devel/init.hh>
#include <numeric/xyzVector.hh>

// Project headers
#include <core/id/NamedAtomID.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/moves/SheetTranslate.hh>

int main(int argc, char* argv[]) {
  using core::id::NamedAtomID;
  using core::pose::Pose;
  using numeric::xyzVector;
  using protocols::loops::Loop;
  using protocols::moves::SheetTranslate;

  devel::init(argc, argv);
  const Pose input = *core::import_pose::pose_from_pdb("/work/tex/casp9_benchmark/meval/fast_cm/T0552/2oxgZ_1.pdb_full_length.pdb");
  Pose output = *core::import_pose::pose_from_pdb("/work/tex/casp9_benchmark/meval/fast_cm/T0552/2oxgZ_1.pdb_full_length.pdb");

  // Translate the specified sheet
  Loop sheet(51, 53);
  double dist = 4.3;

  SheetTranslate mover(sheet, dist);
  mover.apply(output);

  for (unsigned i = sheet.start(); i <= sheet.stop(); ++i) {
    xyzVector<double> xyz_input = input.xyz(NamedAtomID("CA", i));
    xyzVector<double> xyz_output = output.xyz(NamedAtomID("CA", i));
    std::cout << "Distance between residue " << i << " in input and output => "
              << xyz_input.distance(xyz_output) << std::endl;
  }

  // Write result to disk
  core::io::pdb::dump_pdb(output, "output.pdb");
}
