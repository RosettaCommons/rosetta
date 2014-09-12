// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @brief Definition of a filter to extract constellations containing aromatic rings
/// @author Andrea Bazzoli

#include "AromaticFilter.hh"
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>

namespace devel {
namespace constel {

using core::conformation::Residue;
using core::chemical::aa_his;

///
/// @brief: returns true if the given constellation contains at least one aromatic
/// 	ring; returns false otherwise.
///
/// @param[in] ps pose containing the constellation
/// @param[in] cnl pose indexes of the residues forming the constellation
///
/// @details Histidine is considered to be aromatic.
///
/// @remarks It is assumed that:
/// 	1. The residues forming the constellation have non-zero occupancy only for
///   	the atoms that belong to the constellation. This is guaranteed if the
/// 		residues forming the constellation had their indexes previously passed
/// 		as arguments to function "SingResCnlCrea::zero_occ_for_deleted_atoms()".
///
bool has_aromatic(Pose const& ps, utility::vector1<Size> const& cnl) {

	Size const N = cnl.size();
	for(Size i=1; i<=N; ++i) {

		Size const ri = cnl[i];
		Residue const& res = ps.residue(ri);

		if(res.aa() == aa_his)
			return true;

		if(res.is_aromatic()) {
			Size const cgi = res.atom_index("CG");
			if(ps.pdb_info()->occupancy( ri, cgi))
				return true;
		}
	}

	return false;
}

} // constel
} // devel
