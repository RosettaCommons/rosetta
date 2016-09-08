// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/cmiles/fragment_rmsd.cc
/// @author Christopher Miles (cmiles@uw.edu)

// C/C++ headers
#include <iostream>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <devel/init.hh>

// Project headers
#include <core/types.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragmentRmsd.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>

using core::fragment::FragSetCOP;
using core::pose::PoseCOP;

void check_lengths(FragSetCOP fragments, PoseCOP pose) {
	if ( fragments->max_pos() != pose->size() ) {
		std::cerr << "Residues in pose: " << pose->size() << std::endl;
		std::cerr << "Residues in fragments: " << fragments->max_pos() << std::endl;
		std::cerr << "Error: different number of residues in pose and fragments" << std::endl;
		exit(1);
	}
}

int main(int argc, char* argv[]) {
	try {

		using core::Size;
		using core::fragment::FragmentIO;
		using core::fragment::FragmentRmsd;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		devel::init(argc, argv);

		FragSetCOP fragments = core::fragment::FragmentIO().read_data(option[in::file::frag3]());
		PoseCOP pose = core::import_pose::pose_from_file(option[OptionKeys::in::file::s]()[1], core::import_pose::PDB_file);
		check_lengths(fragments, pose);

		const Size last = fragments->max_pos() - fragments->max_frag_length() + 1;

		std::cout << "position frag_idx rmsd" << std::endl;
		FragmentRmsd calc(fragments);
		for ( Size pos = 1; pos <= last; ++pos ) {
			for ( Size k = 1; k <= 200; ++k ) {
				std::cout << pos << " " << k << " " << calc.rmsd(pos, k, *pose) << std::endl;
			}
		}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
