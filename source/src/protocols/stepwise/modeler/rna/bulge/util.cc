// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/rna/bulge/util.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/rna/bulge/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "protocols.stepwise.modeler.rna.bulge.util" );
using namespace core;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace bulge {

//////////////////////////////////////////////////////////////////////
utility::vector1< bool >
detect_base_contacts( core::pose::Pose const & pose ) {

	static Distance const CONTACT_DIST_CUTOFF( 4.0 );
	static Distance const NBR_DIST_CUTOFF( 12.0 );
	static Size const MIN_ATOM_CONTACTS( 5 );

	utility::vector1< bool > base_makes_contact( pose.total_residue(), false );

	// could be dramatically accelerated by using nbr list, etc.
	for ( Size i = 1; i <= pose.total_residue(); i++ ) {
		if ( !pose.residue_type( i ).is_RNA() ) continue;

		bool found_contact( false );
		Size num_contacts( 0 );
		for ( Size ii = pose.residue_type( i ).first_sidechain_atom()+1 /*just nucleobase, no O2'*/;
				ii < pose.residue_type( i ).nheavyatoms(); ii++ ) {
			if ( pose.residue_type( i ).is_virtual( ii ) ) continue;

			bool found_atom_contact( false );
			for ( Size j = 1; j <= pose.total_residue(); j++ ) {
				if ( i == j ) continue;
				if ( ( pose.residue( i ).nbr_atom_xyz() - pose.residue( j ).nbr_atom_xyz() ).length() > NBR_DIST_CUTOFF ) continue;

				for ( Size jj = 1;
						jj < pose.residue_type( i ).nheavyatoms(); jj++ ) {
					if ( pose.residue_type( j ).is_virtual( jj ) ) continue;
					if ( ( pose.residue( i ).xyz( ii ) - pose.residue( j ).xyz( jj ) ).length() < CONTACT_DIST_CUTOFF ) {
						//      TR << "FOUND CONTACT " << pose.pdb_info()->chain(i) << ":" << pose.pdb_info()->number( i ) << " " << pose.residue(i).atom_name(ii)
						//  << " with " << pose.pdb_info()->chain(j) << ":" << pose.pdb_info()->number( j ) << " " << pose.residue(j).atom_name(jj) << std::endl;
						found_atom_contact = true;
						break;
					}
				} // jj
				if ( found_atom_contact ) break;
			} // j
			if ( found_atom_contact ) num_contacts++;
			if ( num_contacts >= MIN_ATOM_CONTACTS ) {
				found_contact = true;
				break;
			}
		}
		if ( found_contact ) {
			base_makes_contact[ i ] = true;
		} else {
			TR.Debug << "BULGE RESIDUE at " << pose.pdb_info()->chain(i) << ":" << pose.residue_type(i).name1() << pose.pdb_info()->number( i ) << std::endl;
		}
	}
	return base_makes_contact;
}

} //bulge
} //rna
} //modeler
} //stepwise
} //protocols
