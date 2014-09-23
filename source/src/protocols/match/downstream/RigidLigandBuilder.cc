// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/downstream/RigidLigandBuilder.cc
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini
/// @author Florian Richter (floric@u.washington.edu)

/// Unit headers
#include <protocols/match/downstream/RigidLigandBuilder.hh>

/// Package headers
#include <protocols/match/BumpGrid.hh>
#include <protocols/match/OccupiedSpaceHash.hh>
#include <protocols/match/downstream/ActiveSiteGrid.hh>
#include <protocols/toolbox/match_enzdes_util/LigandConformer.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <basic/Tracer.hh>

// ObjexxFCL headers
// AUTO-REMOVED #include <ObjexxFCL/FArray2D.hh>

// Utility headers
// AUTO-REMOVED #include <utility/vector1.functions.hh>

#include <core/chemical/AtomType.hh>
#include <core/id/AtomID.hh>
#include <protocols/match/Hit.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace downstream {

static thread_local basic::Tracer TR( "protocols.match.downstream.RigidLigandBuilder" );

RigidLigandBuilder::RigidLigandBuilder() :
	parent(),
	ignore_h_collisions_( false )
{
	//TR << "Initializing RigidLigandBuilder." << std::endl;
}

RigidLigandBuilder::RigidLigandBuilder( RigidLigandBuilder const & other ) :
	parent( other ),
	downstream_restype_( other.downstream_restype_ ),
	upstream_restype_( other.upstream_restype_ ),
	ignore_h_collisions_( other.ignore_h_collisions_ ),
	orientation_atoms_( other.orientation_atoms_ ),
	atoms_123_( other.atoms_123_ ),
	radii_123_( other.radii_123_ ),
	atom_radii_( other.atom_radii_ ),
	atom_required_in_active_site_( other.atom_required_in_active_site_ ),
	non_collision_detection_atoms_reqd_in_active_site_( other.non_collision_detection_atoms_reqd_in_active_site_ ),
	lig_conformers_( other.lig_conformers_.size() ),
	min_sep_d2_from_upstream_atoms_( other.min_sep_d2_from_upstream_atoms_ )
{
  for ( Size ii = 1; ii <= lig_conformers_.size(); ++ii ) {
    lig_conformers_[ ii ] = new toolbox::match_enzdes_util::LigandConformer( * other.lig_conformers_[ ii ] );
  }
}

//RigidLigandBuilder::RigidLigandBuilder( RigidLigandBuilder const & other, core::chemical::ResidueTypeCOP upstream_restype ) :
//	parent(other),
//	downstream_restype_( other.downstream_restype_ ),
//	upstream_restype_( upstream_restype ),
//	ignore_h_collisions_( other.ignore_h_collisions_ ),
//	atom_radii_( other.atom_radii_ ),
//  atom_required_in_active_site_( other.atom_required_in_active_site_ ),
//  non_collision_detection_atoms_reqd_in_active_site_( other.non_collision_detection_atoms_reqd_in_active_site_ ),
//  lig_conformers_( other.lig_conformers_.size() )
//{
//  for ( Size ii = 1; ii <= lig_conformers_.size(); ++ii ) {
//    lig_conformers_[ ii ] = new LigandConformer( * other.lig_conformers_[ ii ] );
//		//TR << "natoms " << lig_conformers_[ii]->n_collision_check_atoms() << std::endl;
//  }
//
//	initialize_upstream_residue( upstream_restype );
//}


RigidLigandBuilder::~RigidLigandBuilder() {}

DownstreamBuilderOP
RigidLigandBuilder::clone() const
{
	return new RigidLigandBuilder( *this );
}


std::list< Hit >
RigidLigandBuilder::build(
	HTReal const & atom3_frame,
	Size const scaffold_build_point_id,
	Size const upstream_conf_id,
	Size const external_geometry_id,
	core::conformation::Residue const & upstream_residue
) const
{
	assert( downstream_restype_ );
	assert( upstream_restype_ );
	assert( bbgrid_set() );
	assert( (& upstream_residue.type()) == upstream_restype_.get() );

	std::list< Hit > hitlist;

	//std::cout << "RigidLigandBuilder::build" << std::endl;
	//for ( Size ii = 1; ii <= 3; ++ii ) {
		//Vector const iiloc = atom3_frame * ats123_in_atom3_frame_[ ii ];
		//std::cout << "Atom D" << ii << " coordinate: " << iiloc.x() << " " << iiloc.y() << " " << iiloc.z() << std::endl;
	//}

	/// collision detection and active-site containment enforcement.
	for ( Size ii = 1; ii <= lig_conformers_[1]->n_collision_check_atoms(); ++ii ) {
		Size ii_restype_id = lig_conformers_[1]->collision_check_id_2_restype_id( ii );

		Vector const iiloc = lig_conformers_[1]->coordinate_in_D3_frame( ii_restype_id, atom3_frame );
		//std::cout << "   " << downstream_restype_->atom_name( at3_frame_id_2_restype_id_[ ii ] ) << " ";
			//std::cout << iiloc.x() << "  " << iiloc.y() << " " << iiloc.z() << std::endl;
		if ( atom_radii_[ ii_restype_id ] > ZERO && bbgrid().occupied( atom_radii_[ ii_restype_id ], iiloc ) ) {
			return hitlist;
		}
		if ( atom_required_in_active_site_[ ii_restype_id ] && ! active_site_grid().occupied( iiloc ) ) {
			return hitlist;
		}

		for ( Size jj = 1; jj <= min_sep_d2_from_upstream_atoms_[ ii_restype_id ].size(); ++jj ) {
			if ( iiloc.distance_squared( upstream_residue.xyz( min_sep_d2_from_upstream_atoms_[ ii_restype_id ][ jj ].first ))
					< min_sep_d2_from_upstream_atoms_[ ii_restype_id ][ jj ].second ) {
				//std::cout << "collision between " << downstream_restype_->atom_name( at3_frame_id_2_restype_id_[ ii ] );
				//std::cout << " on " << downstream_restype_->name() << " with ";
				//std::cout << upstream_restype_->atom_name( min_sep_d2_from_upstream_atoms_[ ii ][ jj ].first ) << " on " << upstream_restype_->name() << std::endl;
				return hitlist;
			}
		}
	}

	/// Check the atoms we require to be in the active site, but which are not used in
	/// collision detection
	for ( Size ii = 1; ii <= non_collision_detection_atoms_reqd_in_active_site_.size(); ++ii ) {
		Size ii_restype_id = non_collision_detection_atoms_reqd_in_active_site_[ ii ];
		Vector const iiloc = lig_conformers_[1]->coordinate_in_D3_frame( ii_restype_id, atom3_frame );

		if ( ! active_site_grid().occupied( iiloc ) ) {
			return hitlist;
		}
	}

	Real6 global_coordinate = lig_conformers_[1]->global_orientation_from_frame3( atom3_frame );

	/// Check, if we're past the first round of hit building, that this orientation's
	/// bin in 6D is not empty.  If the bin is empty, then there is no way this
	/// ligand placement could result in a match.  Do not return this orientation as a hit.
	//std::cout << "global coordinate: ";
	//for ( Size ii = 1; ii <= 6; ++ii ) { std::cout << global_coordinate[ ii ] << " ";}
	//std::cout << std::endl;

	if ( occ_space_set() && ! occ_space().match_possible_for_hit_geometry( global_coordinate ) ) {
		return hitlist;
	}

	/// We have a hit!
	//std::cout << "HIT!" << std::endl;

	Hit hit;
	hit.first()[ 1 ] = scaffold_build_point_id;
	hit.first()[ 2 ] = upstream_conf_id;
	hit.first()[ 3 ] = external_geometry_id;
	hit.first()[ 4 ] = 1; // my state is always 1
	hit.second() = global_coordinate;


	hitlist.push_back( hit ); /// new called here -- otherwise, nothing has been allocated on the heap since the Dunbrack rotamers

	return hitlist;
}

void
RigidLigandBuilder::set_bb_grid(
	BumpGridCOP bbgrid
)
{
	parent::set_bb_grid( bbgrid );
	if ( upstream_restype_ ) {
		initialize_upstream_nonbonded_min_separation_d2();
	}
}

bool
RigidLigandBuilder::compatible(
	Hit const & my_hit,
	DownstreamBuilder const & other,
	Hit const & other_hit,
	bool //first_dispatch
) const
{
	return other.compatible( other_hit, *this, my_hit, false );
}

/// @details RigidLigandBuilder hits are always compatible, because there is only
/// one ligand conformation.
bool
RigidLigandBuilder::compatible(
	Hit const &, // my_hit,
	RigidLigandBuilder const &,// other,
	Hit const &,// other_hit,
	bool// first_dispatch
) const
{
	return true;
}


void
RigidLigandBuilder::require_atom_to_reside_in_active_site(
	core::id::AtomID const & id
)
{
	runtime_assert( id.rsd() == 1 );
	runtime_assert( id.atomno() < downstream_restype_->natoms() );
	if ( lig_conformers_[1]->restype_id_2_collision_check_id( id.atomno() ) == 0 ) {
		non_collision_detection_atoms_reqd_in_active_site_.push_back( id.atomno() );
	} else {
		atom_required_in_active_site_[ id.atomno() ] = true;
		for ( Size ii = 1; ii <= 3; ++ii ) {
			if ( atoms_123_[ ii ] == id.atomno() ) {
				ats123_reqd_in_active_site_[ ii ] = true;
				break;
			}
		}
	}
}

ProbeRadius
RigidLigandBuilder::atom1_radius() const
{
	assert( downstream_restype_ );
	//return ZERO;
	return radii_123_[ 1 ];
}


ProbeRadius
RigidLigandBuilder::atom2_radius() const
{
	assert( downstream_restype_ );
	return radii_123_[ 2 ];
}


ProbeRadius
RigidLigandBuilder::atom3_radius() const
{
	assert( downstream_restype_ );
	return radii_123_[ 3 ];
}


bool
RigidLigandBuilder::atom1_belongs_in_active_site() const
{
	return ats123_reqd_in_active_site_[ 1 ];
}

bool
RigidLigandBuilder::atom2_belongs_in_active_site() const
{
	return ats123_reqd_in_active_site_[ 2 ];
}

bool
RigidLigandBuilder::atom3_belongs_in_active_site() const
{
	return ats123_reqd_in_active_site_[ 3 ];
}

RigidLigandBuilder::Real
RigidLigandBuilder::atom1_atom2_distance() const
{
	assert( downstream_restype_ );
	return lig_conformers_[1]->atom1_atom2_distance();
}


RigidLigandBuilder::Real
RigidLigandBuilder::atom2_atom3_distance() const
{
	assert( downstream_restype_ );
	return lig_conformers_[1]->atom2_atom3_distance();
}

/// @brief Returns an angle in degrees between the three downstream atoms.

RigidLigandBuilder::Real
RigidLigandBuilder::atom1_atom2_atom3_angle() const
{
	assert( downstream_restype_ );
	return lig_conformers_[1]->atom1_atom2_atom3_angle();
}


void
RigidLigandBuilder::coordinates_from_hit(
	Hit const & hit,
	utility::vector1< AtomID > const & atom_indices,
	utility::vector1< Vector > & atom_coords
) const
{
	/*HTReal global_frame = frame_from_global_orientation( hit.second() );
	for ( Size ii = 1; ii <= atom_indices.size(); ++ii ) {
		assert( atom_indices[ ii ].rsd() == 1 );
		atom_coords[ ii ] = global_frame * points_in_global_orintation_frame_[ atom_indices[ ii ].atomno() ];
	}*/
	lig_conformers_[1]->coordinates_from_orientation( hit.second(), atom_indices, atom_coords );
}



core::pose::PoseCOP
RigidLigandBuilder::downstream_pose_from_hit(
	Hit const & hit
) const
{

	core::conformation::Residue lig_res( *downstream_restype_, false );

	HTReal global_frame = lig_conformers_[1]->frame_from_global_orientation( hit.second() );
	for ( Size ii = 1; ii <= lig_res.natoms(); ++ii ) {
		lig_res.set_xyz( ii, lig_conformers_[1]->coordinate_in_global_frame( ii, global_frame ) );
	}

	core::pose::PoseOP pose = new core::pose::Pose();
	pose->append_residue_by_jump( lig_res, 1 );

	//we should also set a different chain for the downstream pose
	core::pose::PDBInfoOP pdbinf = new core::pose::PDBInfo( *pose );
	pose->pdb_info( pdbinf );
	pose->pdb_info()->chain( 1, 'X' );

	return pose;

}

RigidLigandBuilder::Size
RigidLigandBuilder::n_possible_hits_per_at3frame() const
{
	return 1;
}



void
RigidLigandBuilder::initialize_from_residue(
	Size atom1,
	Size atom2,
	Size atom3,
	Size orientation_atom1,
	Size orientation_atom2,
	Size orientation_atom3,
	core::conformation::Residue const & residue
)
{
	atoms_123_[ 1 ] = atom1;
	atoms_123_[ 2 ] = atom2;
	atoms_123_[ 3 ] = atom3;

	orientation_atoms_[ 1 ] = orientation_atom1;
	orientation_atoms_[ 2 ] = orientation_atom2;
	orientation_atoms_[ 3 ] = orientation_atom3;

	Size const natoms = residue.natoms();
	if ( natoms < 3 ) {
		utility_exit_with_message( "ERROR in RigidLigandBuilder: cannot build a residue with fewer than three atoms" );
	}
	downstream_restype_ = residue.type().get_self_ptr();
	atom_radii_.resize( natoms );
	atom_required_in_active_site_.resize( natoms, false );

	for ( Size ii = 1; ii <= natoms; ++ii ) {
		atom_radii_[ ii ] = probe_radius_for_atom_type( residue.atom( ii ).type() );
	}

	lig_conformers_.push_back( new toolbox::match_enzdes_util::LigandConformer );
	lig_conformers_[1]->ignore_h_collisions( ignore_h_collisions_ );
	lig_conformers_[1]->initialize_from_residue( atom1, atom2, atom3,
		orientation_atom1, orientation_atom2, orientation_atom3, residue );

	for ( Size ii = 1; ii <= 3; ++ii ) {
		if ( ignore_h_collisions_ && residue.atom_type( atoms_123_[ ii ] ).element() == "H" ) {
			radii_123_[ ii ] = ZERO;
		} else {
			radii_123_[ ii ] = probe_radius_for_atom_type( residue.atom( atoms_123_[ ii ] ).type() );
		}
	}
}

void
RigidLigandBuilder::initialize_upstream_residue(
	core::chemical::ResidueTypeCOP upstream_res,
	core::scoring::etable::count_pair::CountPairFunctionCOP count_pair
)
{
	assert( downstream_restype_ );
	assert( upstream_res );

	upstream_restype_ = upstream_res;

	Size const natoms = downstream_restype_->natoms();

	min_sep_d2_from_upstream_atoms_.clear();
	min_sep_d2_from_upstream_atoms_.resize( natoms );

	for ( Size ii = 1; ii <= lig_conformers_[1]->n_collision_check_atoms(); ++ii ) {
		Size ii_restype_id = lig_conformers_[1]->collision_check_id_2_restype_id( ii );

		Size n_to_count( 0 );
		for ( Size jj = upstream_restype_->first_sidechain_atom();
				jj <= upstream_restype_->nheavyatoms(); ++jj ) {
			Real weight( 1.0 );
			Size path_dist( 0 );
			if ( ! count_pair || ( count_pair->count( jj, ii_restype_id, weight, path_dist ) && weight == 1.0 ) ) {
				++n_to_count;
				//std::cout << "Bump check " << downstream_restype_->atom_name( ii_restype_id );
				//std::cout << " on " << downstream_restype_->name() << " with ";
				//std::cout << upstream_restype_->atom_name( jj ) << " on " << upstream_restype_->name() << std::endl;
			}
		}
		/// Now make sure that if we're within count-pair striking distance of the backbone so that we don't
		/// reject a conformation due to this atom as registering a collision in bump_grid.occupied()
		for ( Size jj = 1; jj < upstream_restype_->first_sidechain_atom(); ++jj ) {
			if ( jj > upstream_restype_->natoms() ) break;
			Real weight( 1.0 );
			Size path_dist( 0 );
			if ( count_pair && ( ! count_pair->count( jj, ii_restype_id, weight, path_dist ) || weight != 1.0 )) {
				/// WITHIN STRIKING DISTANCE OF BACKBONE!  DO NOT COLLISON-CHECK THIS ATOM
				atom_radii_[ ii_restype_id ] = ZERO;
			}
		}
		min_sep_d2_from_upstream_atoms_[ ii_restype_id ].resize( n_to_count );
		n_to_count = 0;
		for ( Size jj = upstream_restype_->first_sidechain_atom();
				jj <= upstream_restype_->nheavyatoms(); ++jj ) {
			Real weight( 1.0 );
			Size path_dist( 0 );
			if ( ! count_pair || ( count_pair->count( jj, ii_restype_id, weight, path_dist ) && weight == 1.0 ) ) {
				min_sep_d2_from_upstream_atoms_[ ii_restype_id ][ ++n_to_count ].first = jj;
			}
		}
	}

	/// Now check atoms D1, D2, and D3 and makes sure they are not within striking distance of the backbone
	/// If they are, then set their radii to ZERO
	for ( Size ii = 1; ii <= 3; ++ii ) {
		if ( count_pair ) {
			Size ii_id = atoms_123_[ ii ];
			for ( Size jj = 1; jj < upstream_restype_->first_sidechain_atom(); ++jj ) {
				Real weight( 1.0 );
				Size path_dist( 0 );
				if ( ! count_pair->count( jj, ii_id, weight, path_dist ) || weight != 1.0 ) {
					radii_123_[ ii ] = ZERO;
					break;
				}
			}
		}
	}

	if ( bbgrid_set() ) {
		initialize_upstream_nonbonded_min_separation_d2();
	}

	//if ( count_pair ) {
	//	for ( Size ii = 1; ii <= upstream_restype_->n_
	//}
}

/*Real6
RigidLigandBuilder::global_orientation_from_frame3(
	HTReal const & frame3
) const
{
	HTReal global_frame = frame3 * oframe_in_at3frame_;
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



	if ( false ) {

		Vector oat1_coord = frame3 * oats_in_at3_frame_[ 1 ];
		Vector oat2_coord = frame3 * oats_in_at3_frame_[ 2 ];
		Vector oat3_coord = frame3 * oats_in_at3_frame_[ 3 ];
		HTReal global_frame2( oat1_coord, oat2_coord, oat3_coord );
		Vector euler_angles2 = global_frame2.euler_angles_deg();
		std::cout << "Euler angle comparison.";
		std::cout << " 1: " << euler_angles( 1 ) << " vs " << euler_angles2( 1 ) << " ";
		std::cout << " 2: " << euler_angles( 2 ) << " vs " << euler_angles2( 2 ) << " ";
		std::cout << " 3: " << euler_angles( 3 ) << " vs " << euler_angles2( 3 ) << std::endl;

		for ( Size ii = 1; ii <= points_in_atom3_frame_.size(); ++ii ) {
			Vector ii3loc = frame3 * points_in_atom3_frame_[ ii ];
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
	}

	return global_coords;
}

RigidLigandBuilder::HTReal
RigidLigandBuilder::frame_from_global_orientation(
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

	//Vector euler_deg2 = oframe.euler_angles_deg();
	//std::cout << "reverse euler angles: " << euler_deg2(1) << " " << euler_deg2(2) << " " << euler_deg2(3) << std::endl;

	return oframe;
}
*/

void
RigidLigandBuilder::ignore_h_collisions( bool setting )
{
	if ( downstream_restype_ != 0 ) {
		utility_exit_with_message( "ERROR: ignore_h_collisions_ must be set before the downstream restype is initialized" );
	} else {
		ignore_h_collisions_ = setting;
	}
}


void
RigidLigandBuilder::initialize_upstream_nonbonded_min_separation_d2()
{
	assert( bbgrid_set() );
	for ( Size ii = 1; ii <= lig_conformers_[1]->n_collision_check_atoms(); ++ii ) {
		Size ii_restype_id = lig_conformers_[1]->collision_check_id_2_restype_id( ii );

		ProbeRadius ii_rad = atom_radii_[ ii_restype_id ];
		for ( Size jj = 1; jj <= min_sep_d2_from_upstream_atoms_[ ii_restype_id ].size(); ++jj ) {
			Size upstream_atom_id = min_sep_d2_from_upstream_atoms_[ ii_restype_id ][ jj ].first;
			ProbeRadius jj_rad = probe_radius_for_atom_type( upstream_restype_->atom( upstream_atom_id ).atom_type_index() );

			Real dis = bbgrid().required_separation_distance( ii_rad, jj_rad );
			min_sep_d2_from_upstream_atoms_[ ii_restype_id ][ jj ].second = dis*dis;
		}
	}

}

toolbox::match_enzdes_util::LigandConformerOP
RigidLigandBuilder::get_lig_conformers(core::Size conf_id) const
{
 return lig_conformers_[ conf_id ];
}

//utility::vector1< utility::vector1< std::pair< core::Size, core::Real > > >
//RigidLigandBuilder::get_min_sep_d2_from_upstream_atoms() const
//{
//	return min_sep_d2_from_upstream_atoms_;
//}

core::chemical::ResidueTypeCOP
RigidLigandBuilder::get_upstream_restype() const
{
	return upstream_restype_;
}


}
}
}
