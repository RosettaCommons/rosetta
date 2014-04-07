// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/cmiles/close_loops_ccd.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <devel/init.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <protocols/comparative_modeling/LoopRelaxMover.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

int main(int argc, char* argv[]) {
	try {

  using core::fragment::FragSetOP;
  using core::fragment::FragSetOP;
  using core::pose::PoseOP;
  using protocols::comparative_modeling::LoopRelaxMover;
  using protocols::loops::Loops;
  using protocols::loops::LoopsOP;
  using protocols::simple_moves::SwitchResidueTypeSetMover;
  using utility::vector1;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  devel::init(argc, argv);

  FragSetOP frag_sm = core::fragment::FragmentIO().read_data(option[in::file::frag3]());
  FragSetOP frag_lg = core::fragment::FragmentIO().read_data(option[in::file::frag9]());

  PoseOP pose = core::import_pose::pose_from_pdb(option[OptionKeys::in::file::s]()[1]);
  SwitchResidueTypeSetMover to_centroid(core::chemical::CENTROID);
  to_centroid.apply(*pose);
  pose->dump_pdb("starting_ccd.pdb");

  // Choose chainbreaks automatically
  LoopsOP empty = new protocols::loops::Loops();
  protocols::comparative_modeling::LoopRelaxMover closure;
  closure.remodel("quick_ccd");
  closure.intermedrelax("no");
  closure.refine("no");
  closure.relax("no");
  closure.loops(empty);

  utility::vector1<core::fragment::FragSetOP> fragments;
  fragments.push_back(frag_lg);
  fragments.push_back(frag_sm);
  closure.frag_libs(fragments);

  // Use atom pair constraints when available
  closure.cmd_line_csts(option[constraints::cst_fa_file].user());

  // Simple kinematics
  closure.apply(*pose);
  pose->dump_pdb("ending_ccd.pdb");

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
