// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constel/InterfaceFilter.cc
/// @brief Definition of a filter to extract constellations shared by multiple chains
/// @author Andrea Bazzoli

#include <protocols/constel/InterfaceFilter.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>

namespace protocols {
namespace constel {

using core::pose::Pose;
using core::Size;

/// @brief Returns true if a constellation is shared by multiple chains;
///  returns false otherwise.
///
/// @param[in] ps pose containing the constellation
/// @param[in] cnl pose indexes of the residues forming the constellation
///
bool at_interface(Pose const& ps, utility::vector1<core::Size> const& cnl) {

	char const CID1 = ps.pdb_info()->chain(cnl[1]);
	core::Size N = cnl.size();

	for ( core::Size i=2; i<=N; ++i ) {
		if ( CID1 != ps.pdb_info()->chain(cnl[i]) ) {
			return true;
		}
	}

	return false;
}

} // constel
} // protocols

