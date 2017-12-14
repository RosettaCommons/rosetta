// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/match_enzdes_util/LigandConformer.cc
/// @brief  Implementation for class to hold a ligand conformation
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

/// Unit headers
#include <protocols/toolbox/match_enzdes_util/LigandConformer.hh>

// Project headers
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Residue.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>

// Utility headers
#include <utility/vector1.functions.hh>

#include <core/id/AtomID.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace toolbox {
namespace match_enzdes_util {


LigandConformer::LigandConformer() :
	parent(),
	d12_( 1.0 ),
	d23_( 1.0 ),
	ang123_( 109.5 ),
	ignore_h_collisions_( false )
{}

LigandConformer::LigandConformer( LigandConformer const & /*other*/ ) = default;

LigandConformer::~LigandConformer() = default;

LigandConformer::Real
LigandConformer::atom1_atom2_distance() const
{
	debug_assert( ligand_restype_ );
	return d12_;
}


LigandConformer::Real
LigandConformer::atom2_atom3_distance() const
{
	debug_assert( ligand_restype_ );
	return d23_;
}

/// @brief Returns an angle in degrees between the three downstream atoms.

LigandConformer::Real
LigandConformer::atom1_atom2_atom3_angle() const
{
	debug_assert( ligand_restype_ );
	return ang123_;
}


LigandConformer::Real
LigandConformer::oatom1_oatom2_distance() const
{
	debug_assert( ligand_restype_ );
	return points_in_D3_frame_[ orientation_atoms_[ 1 ] ].distance(
		points_in_D3_frame_[ orientation_atoms_[ 2 ] ] );
}

LigandConformer::Real
LigandConformer::oatom2_oatom3_distance() const
{
	debug_assert( ligand_restype_ );
	return points_in_D3_frame_[ orientation_atoms_[ 2 ] ].distance(
		points_in_D3_frame_[ orientation_atoms_[ 3 ] ] );
}

LigandConformer::Real
LigandConformer::oatom1_oatom2_oatom3_angle() const
{
	debug_assert( ligand_restype_ );
	return numeric::constants::d::radians_to_degrees *
		numeric::angle_radians(
		points_in_D3_frame_[ orientation_atoms_[ 1 ] ],
		points_in_D3_frame_[ orientation_atoms_[ 2 ] ],
		points_in_D3_frame_[ orientation_atoms_[ 3 ] ] );
}


void
LigandConformer::coordinates_from_orientation(
	Real6 const & orientation,
	utility::vector1< core::id::AtomID > const & atom_indices,
	utility::vector1< Vector > & atom_coords
) const
{
	//for ( Size jj=1; jj != points_in_global_orintation_frame_.size(); ++jj ) {
	// std::cout << "atomno " << jj << " " << points_in_global_orintation_frame_[jj][0] << "  " << points_in_global_orintation_frame_[jj][1] << "  " << points_in_global_orintation_frame_[jj][2] << std::endl;
	//}

	//std::cout << "LigandConformer::coordinates from hit begin" << std::endl;
	//std::cout << "hit: " << hit.first()[1] << " "<< hit.first()[2] << " "<< hit.first()[3] << " "<< hit.first()[4] << " " << std::endl;
	//std::cout << "hit: " << hit.second()[1] << " "<< hit.second()[2] << " "<< hit.second()[3] << " "<< hit.second()[4]  << " " << hit.second()[5]  << " "<< hit.second()[6] << " " << std::endl;
	//std::cout << "atom_indices.size() " << atom_indices.size() << std::endl;
	//std::cout << "points_in_global_orintation_frame_.size() " << points_in_global_orintation_frame_.size() << std::endl;

	HTReal global_frame = frame_from_global_orientation( orientation );
	for ( Size ii = 1; ii <= atom_indices.size(); ++ii ) {
		//std::cout << "atom_indices[ ii ].atomno() " << ii << " " << atom_indices[ ii ].atomno() << std::endl;
		debug_assert( atom_indices[ ii ].rsd() == 1 );
		atom_coords[ ii ] = global_frame * points_in_global_orintation_frame_[ atom_indices[ ii ].atomno() ];
		//std::cout << "x: "  << atom_coords[ ii ][0] << " y: " << atom_coords[ ii ][1] << " z: "<< atom_coords[ ii ][2] << std::endl;
	}

	//std::cout << "LigandConformer::coordinates from hit begin" << std::endl;
}

void
LigandConformer::initialize_from_residue(
	Size D1,
	Size D2,
	Size D3,
	Size orientation_atom1,
	Size orientation_atom2,
	Size orientation_atom3,
	core::conformation::Residue const & residue
)
{
	//std::cout << "APL DEBUG LigandConformer::initialize_from_residue begin " << std::endl;

	atoms_123_[ 1 ] = D1;
	atoms_123_[ 2 ] = D2;
	atoms_123_[ 3 ] = D3;

	/// Measure the two distances and one angle needed to describe the geometry of this ligand
	/// to the ClassicMatchAlgorithm.
	d12_ = residue.xyz( D1 ).distance( residue.xyz( D2 ) );
	d23_ = residue.xyz( D2 ).distance( residue.xyz( D3 ) );
	ang123_ = numeric::constants::d::radians_to_degrees * numeric::angle_radians(
		residue.xyz( D1 ), residue.xyz( D2 ), residue.xyz( D3 ) );

	orientation_atoms_[ 1 ] = orientation_atom1;
	orientation_atoms_[ 2 ] = orientation_atom2;
	orientation_atoms_[ 3 ] = orientation_atom3;

	Size const natoms = residue.natoms();
	if ( natoms < 3 ) {
		utility_exit_with_message( "ERROR in LigandConformer: cannot build a residue with fewer than three atoms" );
	}

	ligand_restype_ = residue.type_ptr();
	HTReal D3frame( residue.xyz( D1 ), residue.xyz( D2 ), residue.xyz( D3 ) );
	HTReal oframe( residue.xyz( orientation_atom1 ), residue.xyz( orientation_atom2 ), residue.xyz( orientation_atom3 ) );

	oframe_in_D3frame_ = D3frame.inverse() * oframe;

	//std::cout << "APL DEBUG LigandConformer::initialize_from_residue natoms " << natoms << std::endl;
	points_in_global_orintation_frame_.resize( natoms );
	points_in_D3_frame_.resize( natoms );
	for ( Size ii = 1; ii <= natoms; ++ii ) {
		points_in_global_orintation_frame_[ ii ] = oframe.to_local_coordinate( residue.xyz( ii ) );
		points_in_D3_frame_[ ii ] = D3frame.to_local_coordinate( residue.xyz( ii ) );
	}

	bool const atom1_virtual = residue.atom_type( D1 ).name() == "X";
	bool const atom2_virtual = residue.atom_type( D2 ).name() == "X";
	bool const atom3_virtual = residue.atom_type( D3 ).name() == "X";

	bool atom1_heavy = residue.atom_type( D1 ).element() != "H" && ! atom1_virtual;
	bool atom2_heavy = residue.atom_type( D2 ).element() != "H" && ! atom2_virtual;
	bool atom3_heavy = residue.atom_type( D3 ).element() != "H" && ! atom3_virtual;

	Size nvirt( 0 );
	for ( Size ii = 1; ii <= natoms; ++ii ) if ( residue.atom_type( ii ).element() == "X" ) ++nvirt;

	Size const n_atoms_for_collision_check( ignore_h_collisions_ ?
		residue.nheavyatoms() - ( atom1_heavy ? 1 : 0 ) - ( atom2_heavy ? 1 : 0 ) - ( atom3_heavy ? 1 : 0 ) :
		natoms - nvirt - ( atom1_virtual ? 0 : 1 ) - ( atom2_virtual ? 0 : 1 ) - ( atom3_virtual ? 0 : 1 ));

	collision_check_id_2_restype_id_.resize( n_atoms_for_collision_check );
	restype_id_2_collision_check_id_.resize( natoms, 0 );

	utility::vector1< bool >    selected( natoms, false );
	selected[ D1 ] = selected[ D2 ] = selected[ D3 ] = true;
	create_collcheck_ordering( selected, 1 );

	oats_in_D3_frame_[ 1 ] = D3frame.to_local_coordinate( residue.xyz( orientation_atoms_[ 1 ] ));
	oats_in_D3_frame_[ 2 ] = D3frame.to_local_coordinate( residue.xyz( orientation_atoms_[ 2 ] ));
	oats_in_D3_frame_[ 3 ] = D3frame.to_local_coordinate( residue.xyz( orientation_atoms_[ 3 ] ));

	//std::cout << "APL DEBUG LigandConformer::initialize_from_residue end " << std::endl;

}


LigandConformer::Real6
LigandConformer::global_orientation_from_frame3(
	HTReal const & frame3
) const
{
	HTReal global_frame = frame3 * oframe_in_D3frame_;
	//std::cout.precision( 12 );
	//std::cout << "Global frame" << std::endl;
	//std::cout << "  " << global_frame.xx() << " " << global_frame.yx() << " " << global_frame.zx() << " " << global_frame.px() << std::endl;
	//std::cout << "  " << global_frame.xy() << " " << global_frame.yy() << " " << global_frame.zy() << " " << global_frame.py() << std::endl;
	//std::cout << "  " << global_frame.xz() << " " << global_frame.yz() << " " << global_frame.zz() << " " << global_frame.pz() << std::endl;
	//std::cout.precision( 6 );

	Vector euler_angles = global_frame.euler_angles_deg();
	for ( Size ii = 1; ii <= 3; ++ii ) if ( euler_angles( ii ) < 0 ) euler_angles( ii ) += 360.0;
	Vector oat3_coords = global_frame.point(); //frame3 * points_in_at3_frame_[ restype_id_2_at3_frame_id_[ orientation_atoms_[ 3 ]]];

	Real6 global_coords;
	global_coords[ 1 ] = oat3_coords( 1 );
	global_coords[ 2 ] = oat3_coords( 2 );
	global_coords[ 3 ] = oat3_coords( 3 );
	global_coords[ 4 ] = euler_angles( 1 );
	global_coords[ 5 ] = euler_angles( 2 );
	global_coords[ 6 ] = euler_angles( 3 );


	/*if ( false ) {

	Vector oat1_coord = frame3 * oats_in_at3_frame_[ 1 ];
	Vector oat2_coord = frame3 * oats_in_at3_frame_[ 2 ];
	Vector oat3_coord = frame3 * oats_in_at3_frame_[ 3 ];
	HTReal global_frame2( oat1_coord, oat2_coord, oat3_coord );
	Vector euler_angles2 = global_frame2.euler_angles_deg();
	std::cout << "Euler angle comparison.";
	std::cout << " 1: " << euler_angles( 1 ) << " vs " << euler_angles2( 1 ) << " ";
	std::cout << " 2: " << euler_angles( 2 ) << " vs " << euler_angles2( 2 ) << " ";
	std::cout << " 3: " << euler_angles( 3 ) << " vs " << euler_angles2( 3 ) << std::endl;

	for ( Size ii = 1; ii <= points_in_D3_frame_.size(); ++ii ) {
	Vector ii3loc = frame3 * points_in_D3_frame_[ ii ];
	Vector iigloc = global_frame * points_in_global_orintation_frame_[ at3_frame_id_2_restype_id_[ ii ]];
	for ( Size jj = 1; jj <= 3; ++jj ) {
	std::cout << ii << " " << jj << " " << ii3loc( jj ) << " " << iigloc( jj ) << std::endl;
	}
	}

	HTReal global_frame3;
	global_frame3.from_euler_angles_deg( euler_angles );
	global_frame3.set_point( oat3_coord );
	//std::cout.precision( 12 );
	std::cout << "Global frame3" << std::endl;
	std::cout << "  " << global_frame3.xx() << " " << global_frame3.yx() << " " << global_frame3.zx() << " " << global_frame3.px() << std::endl;
	std::cout << "  " << global_frame3.xy() << " " << global_frame3.yy() << " " << global_frame3.zy() << " " << global_frame3.py() << std::endl;
	std::cout << "  " << global_frame3.xz() << " " << global_frame3.yz() << " " << global_frame3.zz() << " " << global_frame3.pz() << std::endl;
	//std::cout.precision( 6 );
	}*/

	return global_coords;
}

LigandConformer::HTReal
LigandConformer::frame_from_global_orientation(
	Real6 orientation
) const
{
	Vector oat3_coord(orientation[ 1 ], orientation[ 2 ], orientation[ 3 ] );
	Vector euler_deg( orientation[ 4 ], orientation[ 5 ], orientation[ 6 ] );

	//std::cout << "oat3_coord: " << oat3_coord(1) << " " << oat3_coord(2) << " " << oat3_coord(3) << std::endl;
	//std::cout << "euler_deg: " << euler_deg(1) << " " << euler_deg(2) << " " << euler_deg(3) << std::endl;

	HTReal oframe;

	oframe.from_euler_angles_deg( euler_deg );
	oframe.set_point( oat3_coord );

	//std::cout.precision( 12 );
	//std::cout << "Global frame" << std::endl;
	//std::cout << "  " << oframe.xx() << " " << oframe.yx() << " " << oframe.zx() << " " << oframe.px() << std::endl;
	//std::cout << "  " << oframe.xy() << " " << oframe.yy() << " " << oframe.zy() << " " << oframe.py() << std::endl;
	//std::cout << "  " << oframe.xz() << " " << oframe.yz() << " " << oframe.zz() << " " << oframe.pz() << std::endl;
	//std::cout.precision( 6 );

	Vector euler_deg2 = oframe.euler_angles_deg();
	//std::cout << "reverse euler angles: " << euler_deg2(1) << " " << euler_deg2(2) << " " << euler_deg2(3) << std::endl;

	return oframe;
}

void
LigandConformer::ignore_h_collisions( bool setting )
{
	if ( ligand_restype_ != nullptr ) {
		utility_exit_with_message( "ERROR: ignore_h_collisions_ must be set before the downstream restype is initialized" );
	} else {
		ignore_h_collisions_ = setting;
	}
}

void
LigandConformer::move_atoms_to_collcheck_begin(
	utility::vector1< Size > const & restype_atnos_to_move_early
)
{
	Size start_from( 1 );
	utility::vector1< bool > selected( points_in_D3_frame_.size() );
	for ( Size ii = 1; ii <= restype_atnos_to_move_early.size(); ++ii ) {
		Size iiatom = restype_atnos_to_move_early[ ii ];
		if ( restype_id_2_collision_check_id_[ iiatom ] == 0 ) continue;
		selected[ iiatom ] = true;
		++start_from;

		collision_check_id_2_restype_id_[ ii ] = iiatom;
		restype_id_2_collision_check_id_[ iiatom ] = ii;

	}
	for ( Size ii = 1; ii <= 3; ++ii ) {
		if ( selected[ atoms_123_[ ii ]] ) {
			--start_from;
		} else {
			selected[ atoms_123_[ ii ]] = true;
		}
	}
	create_collcheck_ordering( selected, start_from );
}


void
LigandConformer::get_global_coords_as_FArray2D(
	ObjexxFCL::FArray2D< numeric::Real > & coords,
	HTReal const & orientation_frame,
	utility::vector1< core::Size > const & restype_atomnos
) const
{
	for ( core::Size ii = 1; ii <= restype_atomnos.size(); ++ii ) {
		Vector coord( orientation_frame * points_in_global_orintation_frame_[ restype_atomnos[ ii ] ] );
		coords( 1, ii ) = coord[0];
		coords( 2, ii ) = coord[1];
		coords( 3, ii ) = coord[2];
	}
}


void
LigandConformer::create_collcheck_ordering(
	utility::vector1< bool > selected,
	Size count_from
)
{
	Size const natoms = ligand_restype_->natoms();
	Size const n_atoms_for_collision_check = collision_check_id_2_restype_id_.size();

	/// Order the other atoms by their distance from D3 and the atoms already selected.
	ObjexxFCL::FArray2D< Real > atdists( natoms, natoms );
	utility::vector1< Real >    atdist_sums( natoms, 0.0 );
	utility::vector1< bool >    ignore( natoms, false );
	for ( Size ii = 1; ii <= natoms; ++ii ) {
		if ( ignore_h_collisions_ && ligand_restype_->atom_type( ii ).element() == "H" ) ignore[ ii ] = true;
		if ( ligand_restype_->atom_type( ii ).element() == "X" ) ignore[ ii ] = true;
		if ( ignore[ ii ] ) continue;

		for ( Size jj = ii + 1; jj <= natoms; ++jj ) {
			atdists( ii, jj ) = atdists( jj, ii ) = points_in_D3_frame_[ ii ].distance( points_in_D3_frame_[ jj ] );
		}
	}

	/// O( N^3 )... but a one-time cost and way fast.  The question is, though,
	/// whether this shuffling of atoms yeilds a speedup and if so, how much?
	for ( Size ii = count_from; ii <= n_atoms_for_collision_check; ++ii ) {
		/// set the atdist sums to zero.
		std::fill( atdist_sums.begin(), atdist_sums.end(), 0.0);
		for ( Size jj = 1; jj <= natoms; ++jj ) {
			if ( selected[ jj ] || ignore[ jj ] ) {
				atdist_sums[ jj ] = -1.0; // make sure this atom is not the furthest atom.
			} else {
				for ( Size kk = 1; kk <= natoms; ++kk ) {
					if ( ignore[ kk ] ) continue;
					if ( selected[ kk ] ) atdist_sums[ jj ] += atdists( kk, jj );
				}
			}
			//std::cout << ii << " " << jj << " " << atdist_sums[ jj ] << std::endl;
		}

		Size furthest = utility::arg_max( atdist_sums );
		//std::cout << ii << " furthest: " << furthest << " " << ligand_restype_->atom_name( furthest ) << std::endl;
		collision_check_id_2_restype_id_[ ii ] = furthest;
		restype_id_2_collision_check_id_[ furthest ] = ii;
		selected[ furthest ] = true;
	}

}

core::chemical::ResidueType
LigandConformer::get_lig_restype() const{
	return * ligand_restype_;
}


}
}
}
