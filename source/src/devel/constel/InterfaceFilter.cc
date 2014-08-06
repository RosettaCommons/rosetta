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

/// @brief Definition of a filter to extract constellations shared by multiple chains
/// @author Andrea Bazzoli

#include <devel/constel/InterfaceFilter.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>

namespace devel {
namespace constel {


///
/// @brief Returns true if a constellation is shared by multiple chains;
/// 	returns false otherwise.
///
/// @param[in] ps pose containing the constellation
/// @param[in] cnl pose indexes of the residues forming the constellation
///
bool at_interface(Pose const& ps, utility::vector1<Size> const& cnl) {

	char const CID1 = ps.pdb_info()->chain(cnl[1]);
	Size N = cnl.size();

	for(Size i=2; i<=N; ++i) {
		if(CID1 != ps.pdb_info()->chain(cnl[i]))
			return true;
	}

	return false;
}

} // constel
} // devel

