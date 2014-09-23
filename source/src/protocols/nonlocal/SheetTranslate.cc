// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/moves/SheetTranslate.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <protocols/nonlocal/SheetTranslate.hh>

// C/C++ headers
#include <iostream>
#include <string>

// Utility header
#include <basic/Tracer.hh>
#include <numeric/xyzVector.hh>

// Project headers
#include <core/conformation/Conformation.hh>
#include <core/id/NamedAtomID.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/nonlocal/StarTreeBuilder.hh>

// Package headers
#include <protocols/moves/Mover.hh>

//Auto Headers
#include <utility/vector1.hh>
namespace protocols {
namespace nonlocal {

static thread_local basic::Tracer TR( "protocols.nonlocal.SheetTranslate" );

SheetTranslate::SheetTranslate() {
  initialize(protocols::loops::Loop(), 0.0);
}

SheetTranslate::SheetTranslate(const protocols::loops::Loop& sheet, double distance_ang) {
  initialize(sheet, distance_ang);
}

void SheetTranslate::initialize(const protocols::loops::Loop& sheet, double distance) {
  sheet_ = sheet;
  distance_ = distance;
}

void SheetTranslate::apply(core::pose::Pose& pose) {
  using core::id::NamedAtomID;
  using core::kinematics::FoldTree;
  using core::kinematics::Jump;
  using core::kinematics::Stub;
  using numeric::xyzVector;
  using protocols::loops::Loops;
  using protocols::nonlocal::StarTreeBuilder;

  if (!is_valid()) {
    TR.Warning << "SheetTranslate::apply() invoked with invalid or incomplete information." << std::endl;
    TR.Warning << "  sheet_ => " << get_sheet() << std::endl;
    TR.Warning << "  distance_ => " << get_distance() << std::endl;
    return;
  }

  // Retain a copy of the input fold tree, since we're responsible for restoring it
  FoldTree input_tree = pose.fold_tree();

  // Configure new kinematics
  Loops chunks;
  decompose_structure(pose.total_residue(), &chunks);

  StarTreeBuilder builder;
  builder.set_up(chunks, &pose);
  TR.Debug << pose.fold_tree() << std::endl;

  // Define the axis of translation as the vector between the first and last residues of the sheet
  xyzVector<double> axis = pose.xyz(NamedAtomID("CA", sheet_.stop())) - pose.xyz(NamedAtomID("CA", sheet_.start()));

  // Translation along the axis
  unsigned jump_num = jump_containing_sheet(chunks);
  Jump jump = pose.jump(jump_num);
  jump.translation_along_axis(pose.conformation().upstream_jump_stub(jump_num), axis, get_distance());
  pose.set_jump(jump_num, jump);

  // Restore input fold tree
  builder.tear_down(&pose);
  pose.fold_tree(input_tree);
}

unsigned SheetTranslate::jump_containing_sheet(const protocols::loops::Loops& chunks) const {
  for (unsigned i = 1; i <= chunks.num_loop(); ++i) {
    if (chunks[i].start() == sheet_.start())
      return i;
  }
  return 0;  // invalid
}

void SheetTranslate::decompose_structure(unsigned num_residues, protocols::loops::Loops* chunks) const {
  using protocols::loops::Loop;
  assert(chunks);
  assert(num_residues > 0);

  const unsigned start = get_sheet().start();
  const unsigned stop = get_sheet().stop();

  // Residues 1 to (sheet - 1)
  if (start > 1) {
    chunks->add_loop(Loop(1, start - 1));
  }

  // Sheet
  chunks->add_loop(Loop(start, stop));

  // Residues (sheet + 1) to end
  if (stop < num_residues) {
    chunks->add_loop(Loop(stop + 1, num_residues));
  }

  chunks->sequential_order();
}

bool SheetTranslate::is_valid() const {
  return sheet_.start() > 0 && sheet_.start() < sheet_.stop();
}

const protocols::loops::Loop& SheetTranslate::get_sheet() const {
  return sheet_;
}

void SheetTranslate::set_sheet(const protocols::loops::Loop& sheet) {
  sheet_ = sheet;
}

double SheetTranslate::get_distance() const {
  return distance_;
}

void SheetTranslate::set_distance(double distance_ang) {
  distance_ = distance_ang;
}

std::string SheetTranslate::get_name() const {
  return "SheetTranslate";
}

moves::MoverOP SheetTranslate::fresh_instance() const {
  return moves::MoverOP( new SheetTranslate() );
}

moves::MoverOP SheetTranslate::clone() const {
  return moves::MoverOP( new SheetTranslate(*this) );
}

}  // namespace moves
}  // namespace protocols
