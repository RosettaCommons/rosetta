// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author

#include <core/scoring/cryst/util.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <basic/options/option.hh>
#include <basic/Tracer.hh>

#include <iostream>

namespace core {
namespace scoring {
namespace cryst {

static basic::Tracer TS( "core.scoring.cryst.util" );


///////
/////// fix hydrogen b factors to be 1.2 times that of the attached heavyatom
void fix_bfactorsH( core::pose::Pose & pose ) {
	for (Size resid = 1; resid <= pose.total_residue(); ++resid) {
		core::conformation::Residue const &rsd_i = pose.residue(resid);
		Size nheavyatoms = rsd_i.nheavyatoms();
		for (Size atmid = 1; atmid <= nheavyatoms; ++atmid) {
			Real hB = 1.2*pose.pdb_info()->temperature( resid, atmid );
			for (Size hid = rsd_i.attached_H_begin(atmid), hid_end = rsd_i.attached_H_end(atmid);
			          hid <= hid_end; ++hid) {
				pose.pdb_info()->temperature( resid, hid, hB );
			}
		}
	}
	pose.pdb_info()->obsolete( false );
}

///////
/////// fix b factors from missing sidechains
///////    if a sidechain b factor is == 0 and the CA one is not, set it to 1.5 times the CA b factor
void fix_bfactorsMissing( core::pose::Pose & pose ) {
	for (Size resid = 1; resid <= pose.total_residue(); ++resid) {
		if (!pose.residue_type(resid).is_protein()) continue;
		core::conformation::Residue const &rsd_i = pose.residue(resid);
		Real bCA = 1.5 * pose.pdb_info()->temperature( resid, 2 );
		for (Size atmid = rsd_i.first_sidechain_atom(); atmid <= rsd_i.nheavyatoms(); ++atmid) {
			Real bSC = pose.pdb_info()->temperature( resid, atmid );
			if (bSC == 0)
				pose.pdb_info()->temperature( resid, atmid, bCA );
		}
	}
	pose.pdb_info()->obsolete( false );

	// now fix attached H
	fix_bfactorsH( pose );
}


} // namespace cryst
} // namespace scoring
} // namespace core
