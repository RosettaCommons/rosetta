// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief definition of generic functions returning information on a constellation
/// @author Andrea Bazzoli

#include <devel/constel/cnl_info.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>

using core::conformation::Residue;


namespace devel {
namespace constel {

///
/// @brief returns the center of mass of a constellation
///
/// @param[in] cnl pose residue indexes of the residues forming the
///   constellation.
/// @param[in] ps pose to which all residues in the constellation belong.
///
/// @remarks 1. It is assumed that the residues forming the constellation have
///   non-zero occupancy only for the atoms that belong to the constellation.
///   This is guaranteed if the residues forming the constellation were
///   previously passed as arguments to function
/// 	"SingResCnlCrea::zero_occ_for_deleted_atoms()".
///
///
xyzVector<Real> cnl_com(vector1<Size> const &cnl, Pose const &ps) {

  xyzVector<Real> com(0);
  Size nats = 0;
  for(Size i=1; i<=cnl.size(); ++i) {
    Residue const &rsd(ps.residue(cnl[i]));
    for (Size j=1; j<=rsd.natoms(); ++j)
      if (ps.pdb_info()->occupancy(cnl[i], j)) {
        ++nats;
        com += rsd.xyz(j);
      }
  }
  com /= nats;

  return com;
}

} // constel
} // devel 
