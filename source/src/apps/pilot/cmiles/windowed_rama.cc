// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/cmiles/windowed_rama.cc
/// @author Christopher Miles (cmiles@uw.edu)

// C/C++ headers
#include <iostream>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <devel/init.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Ramachandran.hh>

using core::Real;
using core::Size;
using core::fragment::FragSetCOP;
using core::pose::PoseCOP;
using utility::vector1;

void compute_windowed_rama(PoseCOP pose, FragSetCOP fragments, Size window, vector1<Real>* ramas) {
  using core::chemical::AA;
  using core::scoring::Ramachandran;
  assert(pose);
  assert(ramas);
  assert(fragments);

  const Size k = (window - 1) / 2;
  const Size n = pose->total_residue();
  ramas->resize(n);

  Ramachandran scorer;

  for (Size center = (1 + k); center <= (n - k); ++center) {
    Real score = 0;

    for (Size residue = (center - k); residue <= (center + k); ++residue) {
      AA amino = pose->aa(residue);
      Real phi = pose->phi(residue);
      Real psi = pose->psi(residue);

      score += scorer.eval_rama_score_residue(amino, phi, psi);
    }

    (*ramas)[center] = score;
  }
}

int main(int argc, char* argv[]) {
    try {
  using core::fragment::FragmentIO;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  devel::init(argc, argv);

  FragSetCOP fragments = core::fragment::FragmentIO().read_data(option[in::file::frag3]());
  PoseCOP pose = core::import_pose::pose_from_pdb(option[OptionKeys::in::file::s]()[1]);

  const Size window = fragments->max_frag_length();

  vector1<Real> ramas;
  compute_windowed_rama(pose, fragments, window, &ramas);

  for (Size i = 1; i <= ramas.size(); ++i) {
    std::cout << "rama " << i << " " << ramas[i] << std::endl;
  }
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }
        return 0;
}
