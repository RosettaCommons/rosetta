// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//////////////////////////////////////////////////////////////////////
///
/// @brief
/// Compare the rotamer recovery between a native protein and a list of other proteins
///
/// @details
/// This is an implementation taken from James Thompson. I am not even sure he knows I stole it
/// from him. The main function that is called is the get_rotamer_recovery() function. You can
/// pass this function a native pdb and a list of altered pdbs, or just 1 native and 1
/// alterd pdb. The rotamer recovery will be output to the screen. Output looks like:
/// # total = 1
///    resi_idx  nat_bb_bin      pct_bb    nat_rot1    pct_rot1    nat_rot2    pct_rot2    nat_rot3    pct_rot3    nat_rot4    pct_rot4
///           1           E      1.0000           1      1.0000           2      1.0000           1      1.0000         999      0.0000
///           2           B      1.0000           2      1.0000           1      1.0000         999      0.0000         999      0.0000
/// Where the # total is how many proteins compared.
/// resi_idx = residue index
/// nat_bb_bin = dssp naming for bb
/// pct_bb = how many match the bb bins?
/// nat_rot1 = chi 1
/// pct_rot1 = how many are correct
/// If 999 appears, that means that the amino acid does not have that chi angle
///
///
///
/// @author
/// @author James Thompson (original author)
/// @author Steven Combs (moved it to protocols for general use)
///
/////////////////////////////////////////////////////////////////////////


/// @file
/// @brief

#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_RotamerRecovery_hh
#define INCLUDED_protocols_toolbox_pose_metric_calculators_RotamerRecovery_hh

//#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
#include <utility/vector1.hh>

#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace pose_metric_calculators {

class RotamerRecovery{
public:
	RotamerRecovery(){

	}

	char torsion2big_bin(
		float const phi,
		float const psi,
		float const omega
	);

	utility::vector1< char > get_ss( core::pose::Pose & pose );

	void print_rot_vec(
		core::pack::dunbrack::RotVector rot_vec,
		std::ostream & out
	);

	utility::vector1< char > bb_bins_from_pose(
		core::pose::Pose const & pose
	);

	utility::vector1< core::pack::dunbrack::RotVector >
	rots_from_pose(
		core::pose::Pose const & pose
	);


	utility::vector1< utility::vector1< core::Real > >
	chis_from_pose(
		core::pose::Pose const & pose
	);


	void get_rotamer_recovery(core::pose::Pose & native, utility::vector1<core::pose::Pose> & compared_poses);

	void get_rotamer_recovery(core::pose::Pose & native, core::pose::Pose & compared_pose);

private:


};


}
}

#endif /* ROTAMERRECOVERY_HH_ */
