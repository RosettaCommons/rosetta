// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/cmiles/windowed_rmsd.cc
/// @author Christopher Miles (cmiles@uw.edu)

// C/C++ headers
#include <iostream>
#include <string>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <devel/init.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>

using core::Real;
using core::Size;
using core::pose::Pose;
using core::pose::PoseOP;
using core::pose::PoseCOP;
using std::string;
using utility::vector1;

/// @detail Ensure that the lengths of the models to evaluate match the reference
void check_lengths(PoseCOP reference, const vector1<PoseCOP>& models) {
  const Size n = reference->total_residue();

  for (Size i = 1; i <= models.size(); ++i) {
    const Size m = models[i]->total_residue();
    if (m != n) {
      std::cerr << "Residues in reference: " << n << std::endl;
      std::cerr << "Residues in model[" << i << "]: " << m << std::endl;
      utility_exit_with_message("Reference structure and model have unequal number of residues");
    }
  }
}

void compute_windowed_rmsd(const Pose& reference, const Pose& model, Size window, vector1<Real>* rmsds) {
  assert(rmsds);
  const Size k = (window - 1) / 2;
  const Size n = reference.total_residue();

  rmsds->resize(n);
  for (Size residue = 1 + k; residue <= (n - k); ++residue) {
    (*rmsds)[residue] = core::scoring::CA_rmsd(reference, model, residue - k, residue + k);
  }
}

void show(const string& filename, const vector1<Real>& rmsds) {
  for (Size i = 1; i <= rmsds.size(); ++i) {
    std::cout << filename << " " << i << " " << rmsds[i] << std::endl;
  }
}

int main(int argc, char* argv[]) {
  try {
    using namespace basic::options;
  using namespace basic::options::OptionKeys;
  devel::init(argc, argv);

  PoseOP reference = core::import_pose::pose_from_pdb(option[OptionKeys::in::file::native]());
  vector1<PoseOP> models = core::import_pose::poseOPs_from_pdbs(option[OptionKeys::in::file::s]());
  check_lengths(reference, models);

  const Size window = option[OptionKeys::evaluation::window_size]();

  std::cout << "filename resi rmsd" << std::endl;
  for (Size i = 1; i <= models.size(); ++i) {
    const string& filename = option[OptionKeys::in::file::s]()[i];
    const Pose& model = *(models[i]);

    vector1<Real> rmsds;
    compute_windowed_rmsd(*reference, model, window, &rmsds);
    show(filename, rmsds);
  }
  } catch ( utility::excn::EXCN_Base const & e ) {
                            std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                }
      return 0;
}
