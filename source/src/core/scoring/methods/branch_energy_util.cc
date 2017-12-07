// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/branch_energyutil.cc
/// @brief  Utility functions for scoring branches.
/// @author Andrew Watkins

#include <core/scoring/methods/branch_energy_util.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <utility/vector1.hh>

#include <core/types.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/conformation/Residue.fwd.hh>


namespace core {
namespace scoring {
namespace methods {

void
find_relevant_connections_onersd( pose::Pose const & pose, Size const seqpos, ResidueAtomOverlaps & branch_connection ) {

	for ( Size ii = 1; ii <= pose.size(); ++ii ) {

		if ( ii == seqpos ) {
			// 'lower trial'
			auto const & ii_rsd = pose.residue( ii );


			if ( ii_rsd.type().has_variant_type( core::chemical::FIVEPRIME_CAP ) ) {
				// OK, this residue has FIVEPRIME_CAP. That means that it has SOME connection
				// whose connection atom is ZO3'
				for ( Size k = 1; k <= ii_rsd.connect_map_size(); k++ ) {
					if ( ii_rsd.atom_name( ii_rsd.residue_connect_atom_index( k ) ) != "ZO3'" ) continue;
					Size other( ii_rsd.connected_residue_at_resconn( k ) );
					if ( other == 0 ) continue;
					auto const & other_rsd = pose.residue( other );

					Size const m = ii_rsd.residue_connection_conn_id( k );
					if ( other_rsd.atom_name( other_rsd.residue_connect_atom_index( m ) ) != " P  " ) continue;

					branch_connection.res1 = ii;
					branch_connection.res2 = other;
					branch_connection.res1_ovl1_overlaps = "P";
					branch_connection.res1_ovl2_overlaps = "O5'";
					branch_connection.res2_ovu1_overlaps = "ZO3'";
					return;
				}
			} else if ( ii_rsd.type().has_variant_type( core::chemical::C2_BRANCH_POINT ) ) {
				// OK, this residue has C2_BRANCH_POINT. That means that it has SOME connection
				// whose connection atom is O2'
				for ( Size k = 1; k <= ii_rsd.connect_map_size(); k++ ) {
					if ( ii_rsd.atom_name( ii_rsd.residue_connect_atom_index( k ) ) != " O2'" ) continue;
					Size other( ii_rsd.connected_residue_at_resconn( k ) );
					if ( other == 0 ) continue;
					auto const & other_rsd = pose.residue( other );

					Size const m = ii_rsd.residue_connection_conn_id( k );
					if ( other_rsd.atom_name( other_rsd.residue_connect_atom_index( m ) ) != " P  " ) continue;

					branch_connection.res1 = ii;
					branch_connection.res2 = other;
					branch_connection.res1_ovl1_overlaps = "P";
					branch_connection.res1_ovl2_overlaps = "O5'";
					branch_connection.res2_ovu1_overlaps = "O2'";
					return;
				}
			}
		} else {
			// 'upper trial'
			for ( Size jj = 1; jj <= pose.size(); ++jj ) {
				if ( pose.residue_type( jj ).has_variant_type( core::chemical::FIVEPRIME_CAP ) ) {

					auto const & jj_rsd = pose.residue( jj );

					// OK, this residue has FIVEPRIME_CAP. That means that it has SOME connection
					// whose connection atom is ZO3'
					for ( Size k = 1; k <= jj_rsd.connect_map_size(); k++ ) {
						if ( jj_rsd.atom_name( jj_rsd.residue_connect_atom_index( k ) ) != "ZO3'" ) continue;

						if ( jj_rsd.connected_residue_at_resconn( k ) != ii ) break;
						auto const & other_rsd = pose.residue( jj );

						Size const m = jj_rsd.residue_connection_conn_id( k );
						if ( other_rsd.atom_name( other_rsd.residue_connect_atom_index( m ) ) != " P  " ) continue;

						branch_connection.res1 = jj;
						branch_connection.res2 = ii;
						branch_connection.res1_ovl1_overlaps = "P";
						branch_connection.res1_ovl2_overlaps = "O5'";
						branch_connection.res2_ovu1_overlaps = "ZO3'";
						return;
					}
				} else if ( pose.residue_type( jj ).has_variant_type( core::chemical::C2_BRANCH_POINT ) ) {

					auto const & jj_rsd = pose.residue( jj );

					// OK, this residue has C2_BRANCH_POINT. That means that it has SOME connection
					// whose connection atom is O2'
					for ( Size k = 1; k <= jj_rsd.connect_map_size(); k++ ) {
						if ( jj_rsd.atom_name( jj_rsd.residue_connect_atom_index( k ) ) != " O2'" ) continue;

						if ( jj_rsd.connected_residue_at_resconn( k ) != ii ) break;
						auto const & other_rsd = pose.residue( jj );

						Size const m = jj_rsd.residue_connection_conn_id( k );
						if ( other_rsd.atom_name( other_rsd.residue_connect_atom_index( m ) ) != " P  " ) continue;

						branch_connection.res1 = jj;
						branch_connection.res2 = ii;
						branch_connection.res1_ovl1_overlaps = "P";
						branch_connection.res1_ovl2_overlaps = "O5'";
						branch_connection.res2_ovu1_overlaps = "O2'";
						return;
					}
				}
			}
		}
	}
}

void
find_relevant_connections( pose::Pose const & pose, utility::vector1< ResidueAtomOverlaps > & branch_connections ) {
	//auto const & fold_tree = pose.fold_tree();
	// Later -- check to make sure FT doesn't run through any of these.

	// For now -- explicitly look for residues with KEY VARIANT TYPES
	// then check what they're bonded to by KEY CONNECTIONS.
	//
	// This architecture assumes that we loop through all rsds and only care about
	// THIS variant. Another option is to check for the variants of each residue,
	// THEN loop through that set, THEN deploy action for each variant that 'is relevant'.
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		if ( pose.residue_type( ii ).has_variant_type( core::chemical::FIVEPRIME_CAP ) ) {

			auto const & ii_rsd = pose.residue( ii );

			// OK, this residue has FIVEPRIME_CAP. That means that it has SOME connection
			// whose connection atom is ZO3'
			for ( Size k = 1; k <= ii_rsd.connect_map_size(); k++ ) {
				if ( ii_rsd.atom_name( ii_rsd.residue_connect_atom_index( k ) ) != "ZO3'" ) continue;
				Size other( ii_rsd.connected_residue_at_resconn( k ) );
				if ( other == 0 ) continue;

				auto const & other_rsd = pose.residue( other );
				Size const m = ii_rsd.residue_connection_conn_id( k );

				// SPACES?!?!?!
				if ( other_rsd.atom_name( other_rsd.residue_connect_atom_index( m ) ) != " P  " ) continue;

				ResidueAtomOverlaps foo;
				foo.res1 = ii;
				foo.res2 = other;
				foo.res1_ovl1_overlaps = "P";
				foo.res1_ovl2_overlaps = "O5'";
				foo.res2_ovu1_overlaps = "ZO3'";
				branch_connections.push_back( foo );
				break;
			}
		} else if ( pose.residue_type( ii ).has_variant_type( core::chemical::C2_BRANCH_POINT ) ) {

			auto const & ii_rsd = pose.residue( ii );

			// OK, this residue has C2_BRANCH_POINT. That means that it has SOME connection
			// whose connection atom is O2'
			for ( Size k = 1; k <= ii_rsd.connect_map_size(); k++ ) {
				if ( ii_rsd.atom_name( ii_rsd.residue_connect_atom_index( k ) ) != " O2'" ) continue;
				Size other( ii_rsd.connected_residue_at_resconn( k ) );
				if ( other == 0 ) continue;

				auto const & other_rsd = pose.residue( other );
				Size const m = ii_rsd.residue_connection_conn_id( k );

				// SPACES?!?!?!
				if ( other_rsd.atom_name( other_rsd.residue_connect_atom_index( m ) ) != " P  " ) continue;

				ResidueAtomOverlaps foo;
				foo.res1 = ii;
				foo.res2 = other;
				foo.res1_ovl1_overlaps = "P";
				foo.res1_ovl2_overlaps = "O5'";
				foo.res2_ovu1_overlaps = "O2'";
				branch_connections.push_back( foo );
				break;
			}
		}
	}
}

} // namespace methods
} // namespace scoring
} // namespace core

