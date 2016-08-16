// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/upstream/ProteinUpstreamBuilder.cc
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

// Unit headers
#include <protocols/match/upstream/UpstreamResTypeGeometry.hh>

// Package headers

// Project headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/Atom.hh>

// Utility headers
#include <utility/exit.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace upstream {

/// @details Auto-generated virtual destructor
UpstreamResTypeGeometry::~UpstreamResTypeGeometry() {}


UpstreamResTypeGeometry::UpstreamResTypeGeometry() :
	restype_name_( "UNINITIALIZED" ),
	N_atom_id_( 0 ),
	CA_atom_id_( 0 ),
	C_atom_id_( 0 ),
	O_atom_id_( 0 ),
	CB_atom_id_( 0 ),
	H_atom_id_( 0 ),
	HA_atom_id_( 0 )
{}

UpstreamResTypeGeometry::UpstreamResTypeGeometry( core::chemical::ResidueType const & res ) :
	restype_name_( res.name() ),
	N_atom_id_( 0 ),
	CA_atom_id_( 0 ),
	C_atom_id_( 0 ),
	O_atom_id_( 0 ),
	CB_atom_id_( 0 ),
	H_atom_id_( 0 ),
	HA_atom_id_( 0 )
{
	initialize_from_residue_type( res );
}

void
UpstreamResTypeGeometry::initialize_from_residue_type(
	core::chemical::ResidueType const & res
)
{
	if ( restype_name_ != res.name() ) {
		restype_name_ = res.name();
	}

	/// 1. Resize arrays that depend on the number of atoms
	Size const n_atoms = res.natoms();

	controlling_chi_for_atom_ =  res.last_controlling_chi();
	which_point_for_atom_.resize( n_atoms );
	std::fill( which_point_for_atom_.begin(), which_point_for_atom_.end(), 0 );
	restype_atom_id_2_nonchi_atom_id_.resize( n_atoms );
	std::fill( restype_atom_id_2_nonchi_atom_id_.begin(), restype_atom_id_2_nonchi_atom_id_.end(), 0 );
	nonchi_atom_id_2_restype_atom_id_.resize( 0 );
	nonchi_atoms_in_ideal_frame_.resize( 0 );

	/// 2. Resize arrays that depend on the number of chi
	Size const n_chi = res.nchi();

	chitip_atoms_.resize( n_chi );
	std::fill( chitip_atoms_.begin(), chitip_atoms_.end(), 0 );

	pre_chitip_transforms_.resize( n_chi );
	for ( Size ii = 1; ii <= n_chi; ++ii ) pre_chitip_transforms_[ ii ].set_identity();

	ht_for_chitip_atoms_.resize( n_chi );
	for ( Size ii = 1; ii <= n_chi; ++ii ) ht_for_chitip_atoms_[ ii ].set_identity();

	nonchitip_atoms_.resize( res.nchi() );
	for ( Size ii = 1; ii <= n_chi; ++ii ) nonchitip_atoms_[ ii ].clear();

	points_for_nonchitip_atoms_.resize( res.nchi() );
	for ( Size ii = 1; ii <= n_chi; ++ii ) points_for_nonchitip_atoms_[ ii ].clear();

	/// Quick atom indexing

	N_atom_id_  = ( res.has( "N"  ) ? res.atom_index("N")  : 0 );
	CA_atom_id_ = ( res.has( "CA" ) ? res.atom_index("CA") : 0 );
	C_atom_id_  = ( res.has( "C"  ) ? res.atom_index("C")  : 0 );
	O_atom_id_  = ( res.has( "O"  ) ? res.atom_index("O")  : 0 );
	CB_atom_id_ = ( res.has( "CB" ) ? res.atom_index("CB") : 0 );
	H_atom_id_  = ( res.has( "H"  ) ? res.atom_index("H")  : 0 );
	HA_atom_id_ = ( res.has( "HA" ) ? res.atom_index("HA") : 0 );

	if ( HA_atom_id_ == 0 && res.aa() == core::chemical::aa_gly ) {
		HA_atom_id_ = ( res.has( "2HA" ) ? res.atom_index( "2HA" ) : 0 );
	}

	if ( N_atom_id_ != 0 && CA_atom_id_ != 0 && C_atom_id_ != 0 ) {

		Vector halfpoint = 0.5 * ( res.atom( N_atom_id_ ).ideal_xyz() + res.atom( C_atom_id_ ).ideal_xyz() );
		HTReal ideal_frame( res.atom( N_atom_id_ ).ideal_xyz(), halfpoint, res.atom( CA_atom_id_ ).ideal_xyz() );

		/// backbone atoms, besides the cannonical 3
		for ( Size ii = 1; ii <= n_atoms; ++ii ) {
			if ( ii == N_atom_id_ || ii == CA_atom_id_ || ii == C_atom_id_ ) continue;
			if ( res.last_controlling_chi( ii ) != 0 ) continue;
			nonchi_atoms_in_ideal_frame_.push_back( ideal_frame.to_local_coordinate( res.atom( ii ).ideal_xyz() ));
			nonchi_atom_id_2_restype_atom_id_.push_back( ii );
			restype_atom_id_2_nonchi_atom_id_[ ii ] = nonchi_atoms_in_ideal_frame_.size();
		}
	}


	if ( nchi() == 0 ) return; // match from gly? can't see why you'd want to!


	for ( Size ii = 1; ii <= n_chi; ++ii ) {
		assert( res.chi_atoms( ii ).size() == 4 );

		Size const
			chiat1( res.chi_atoms( ii )[ 1 ] ),
			chiat2( res.chi_atoms( ii )[ 2 ] ),
			chiat3( res.chi_atoms( ii )[ 3 ] ),
			chiat4( res.chi_atoms( ii )[ 4 ] );

		if ( ii > 1 && chiat3 != chitip_atoms_[ ii-1 ] ) {
			// This is a chi dihedral defined by atoms that are rigidly controlled by
			// an intervening stretch of other atoms.  What is needed is the coordinate
			// frame that transforms from a frame defined at the chitip(ii-1) atoms from
			// atoms 2, 3, and 4 of chi(ii-1) to the frame at the third atom of chi(ii).

			Size const prevchi_at2( res.chi_atoms(ii-1)[2] ),
				prevchi_at3( res.chi_atoms(ii-1)[3] ),
				prevchi_at4( res.chi_atoms(ii-1)[4] );
			HTReal launch_frame( res.atom(prevchi_at2).ideal_xyz(),
				res.atom(prevchi_at3).ideal_xyz(),
				res.atom(prevchi_at4).ideal_xyz() );
			HTReal landing_frame( res.atom(chiat1).ideal_xyz(),
				res.atom(chiat2).ideal_xyz(),
				res.atom(chiat3).ideal_xyz());

			pre_chitip_transforms_[ ii ] = launch_frame.inverse() * landing_frame;

		}


		chitip_atoms_[ ii ] = chiat4;


		ht_for_chitip_atoms_[ ii ].set_xaxis_rotation_rad( -1 * res.icoor( chiat4 ).theta() );
		ht_for_chitip_atoms_[ ii ].walk_along_z( res.icoor( chiat4 ).d() );

		HTReal chi_tip_frame(
			res.atom( chiat2 ).ideal_xyz(),
			res.atom( chiat3 ).ideal_xyz(),
			res.atom( chiat4 ).ideal_xyz() );

		Size const n_nontip_ats_for_chi = res.atoms_last_controlled_by_chi( ii ).size() - 1;

		nonchitip_atoms_[ ii ].reserve( n_nontip_ats_for_chi );
		points_for_nonchitip_atoms_[ ii ].reserve( n_nontip_ats_for_chi );

		for ( Size jj = 1; jj <= res.atoms_last_controlled_by_chi( ii ).size(); ++jj ) {
			Size const jjatom = res.atoms_last_controlled_by_chi( ii )[ jj ];
			if ( jjatom == chiat4 ) continue;

			Vector jjloc_in_chitip_frame = chi_tip_frame.to_local_coordinate( res.atom( jjatom ).ideal_xyz() );
			nonchitip_atoms_[ ii ].push_back( jjatom );
			points_for_nonchitip_atoms_[ ii ].push_back( jjloc_in_chitip_frame );
			which_point_for_atom_[ jjatom ] = points_for_nonchitip_atoms_[ ii ].size();
		}
	}

}

bool
UpstreamResTypeGeometry::atom_has_nonchi_coordinate( Size restype_atomid ) const
{
	return restype_atom_id_2_nonchi_atom_id_[ restype_atomid ] != 0;
}

UpstreamResTypeGeometry::Vector const &
UpstreamResTypeGeometry::coordinate_for_nonchi_atom_in_ideal_frame( Size restype_atomid ) const
{
	runtime_assert( restype_atom_id_2_nonchi_atom_id_[ restype_atomid ] != 0 );
	return nonchi_atoms_in_ideal_frame_[ restype_atom_id_2_nonchi_atom_id_[ restype_atomid ] ];
}


}
}
}
