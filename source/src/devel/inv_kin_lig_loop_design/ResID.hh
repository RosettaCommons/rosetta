// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/InvKinLigLoopDesign/ResID.hh
///
/// @brief
/// @author

#ifndef DEVEL_INVKINLIGLOOPDESIGN_RESID_HH
#define DEVEL_INVKINLIGLOOPDESIGN_RESID_HH

#include <iosfwd>
#include <map>
#include <string>

#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>


namespace devel {

namespace inv_kin_lig_loop_design {

using namespace std;

// ===============================================
// ==================== ResID ====================
// ===============================================

struct ResID;
typedef map<ResID,core::conformation::Residue*> resids_type;

struct ResID {
	std::string chain_id;
	int res_num;

	ResID();
	ResID( std::string const & chain_id, int res_num );
	ResID( core::conformation::Residue const& r );

	bool operator==(ResID const& other ) const;
	bool operator<(ResID const& other ) const;

	static resids_type get_resids( core::pose::Pose& pose );
};

istream& operator>>(istream& in, ResID& res_id );
ostream& operator<<(ostream& out, ResID const& res_id );


}

}

#endif // DEVEL_LOOPDESIGN_RESID_HH
