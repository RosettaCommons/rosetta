// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/downstream/LigandConformerBuilder.cc
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

/// Unit headers
#include <protocols/match/downstream/LigandConformerBuilder.hh>

/// Package headers
#include <protocols/match/BumpGrid.hh>
#include <protocols/match/OccupiedSpaceHash.hh>
#include <protocols/match/downstream/ActiveSiteGrid.hh>

// Project headers
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <protocols/toolbox/match_enzdes_util/LigandConformer.hh>


#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/pack/rotamers/SingleLigandRotamerLibrary.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>

#include <protocols/idealize/IdealizeMover.hh>

#include <basic/Tracer.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/string_util.hh>

//numeric headers
#include <numeric/model_quality/rms.hh>

#include <core/id/AtomID.hh>
#include <protocols/match/Hit.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace downstream {

static THREAD_LOCAL basic::Tracer TR( "protocols.match.downstream.LigandConformerBuilder" );

LigandConformerBuilder::LigandConformerBuilder() :
	parent(),
	ignore_h_collisions_( false ),
	idealize_conformers_( true ),
	rmsd_unique_cutoff_(0.05)
{
	//std::cout << "APL DEBUG LigandConformerBuilder ctor " << this << std::endl;
}

LigandConformerBuilder::LigandConformerBuilder( LigandConformerBuilder const & other ) :
	parent( other ),
	downstream_restype_( other.downstream_restype_ ),
	upstream_restype_( other.upstream_restype_ ),
	ignore_h_collisions_( other.ignore_h_collisions_ ),
	idealize_conformers_( other.idealize_conformers_ ),
	orientation_atoms_( other.orientation_atoms_ ),
	atoms_123_( other.atoms_123_ ),
	radii_123_( other.radii_123_ ),
	ats123_reqd_in_active_site_( other.ats123_reqd_in_active_site_ ),
	atom_radii_( other.atom_radii_ ),
	atom_required_in_active_site_( other.atom_required_in_active_site_ ),
	non_collision_detection_atoms_reqd_in_active_site_( other.non_collision_detection_atoms_reqd_in_active_site_ ),
	rmsd_unique_cutoff_( other.rmsd_unique_cutoff_ ),
	conformer_group_indices_( other.conformer_group_indices_ ),
	lig_conformers_( other.lig_conformers_.size() ),
	min_sep_d2_from_upstream_atoms_( other.min_sep_d2_from_upstream_atoms_ )
{
	for ( Size ii = 1; ii <= lig_conformers_.size(); ++ii ) {
		lig_conformers_[ ii ] = toolbox::match_enzdes_util::LigandConformerOP( new toolbox::match_enzdes_util::LigandConformer( * other.lig_conformers_[ ii ] ) );
	}
	//std::cout << "APL DEBUG LigandConformerBuilder copy ctor " << this << std::endl;
}

//LigandConformerBuilder::LigandConformerBuilder( LigandConformerBuilder const & other, core::chemical::ResidueTypeCOP upstream_restype ) :
//  parent(other),
//  downstream_restype_( other.downstream_restype_ ),
//  upstream_restype_( upstream_restype ),
//  ignore_h_collisions_( other.ignore_h_collisions_ ),
//  atom_radii_( other.atom_radii_ ),
//  atom_required_in_active_site_( other.atom_required_in_active_site_ ),
//  non_collision_detection_atoms_reqd_in_active_site_( other.non_collision_detection_atoms_reqd_in_active_site_ ),
//  lig_conformers_( other.lig_conformers_.size() )
//{
//  for ( Size ii = 1; ii <= lig_conformers_.size(); ++ii ) {
//    lig_conformers_[ ii ] = new LigandConformer( * other.lig_conformers_[ ii ] );
//  //TR << "Nr of n_collision_check_atoms: " << lig_conformers_[ ii ]->n_collision_check_atoms()  << std::endl;
//  runtime_assert(lig_conformers_[ ii ]);
//  }
//
//  initialize_upstream_residue( upstream_restype );
//}


LigandConformerBuilder::~LigandConformerBuilder() {}

DownstreamBuilderOP
LigandConformerBuilder::clone() const
{
	return DownstreamBuilderOP( new LigandConformerBuilder( *this ) );
}


std::list< Hit >
LigandConformerBuilder::build(
	HTReal const & atom3_frame,
	Size const scaffold_build_point_id,
	Size const upstream_conf_id,
	Size const external_geometry_id,
	core::conformation::Residue const & upstream_residue
) const
{
	std::list< Hit > hitlist;
	for ( Size ii = 1; ii <= conformer_group_indices_.size(); ++ii ) {
		std::list< Hit > ii_hits = build_conformer_group( ii, atom3_frame, scaffold_build_point_id,
			upstream_conf_id, external_geometry_id, upstream_residue );
		hitlist.splice( hitlist.end(), ii_hits );
	}
	return hitlist;
}

void
LigandConformerBuilder::determine_redundant_conformer_groups(
	utility::vector1< core::Size > const & relevant_atom_indices
)
{
	Size num_relevant_atoms( relevant_atom_indices.size() );
	conformer_group_indices_.clear();
	conformer_group_indices_.push_back( utility::vector1< Size >() );
	conformer_group_indices_[1].push_back( 1 ); //first one is always unique
	conformer_group_for_conformer_.clear();
	conformer_group_for_conformer_.push_back( 1 ); //the conformer group for the first conformer is always 1
	HTReal identity_ht;

	for ( Size ii = 2; ii <= lig_conformers_.size(); ++ii ) {

		ObjexxFCL::FArray2D< numeric::Real > queryconf_coord( 3, num_relevant_atoms );
		lig_conformers_[ii]->get_global_coords_as_FArray2D( queryconf_coord, identity_ht, relevant_atom_indices );

		bool conformer_unique(true);
		for ( Size jj = 1; jj <= conformer_group_indices_.size(); ++jj ) {

			ObjexxFCL::FArray2D< numeric::Real > uniqueconf_coord( 3, num_relevant_atoms );
			lig_conformers_[ conformer_group_indices_[jj][1] ]->get_global_coords_as_FArray2D( uniqueconf_coord, identity_ht, relevant_atom_indices );

			Real rmsd_this_pair = numeric::model_quality::rms_wrapper( num_relevant_atoms, queryconf_coord, uniqueconf_coord );
			if ( rmsd_this_pair < rmsd_unique_cutoff_ ) {
				conformer_group_indices_[jj].push_back( ii );
				conformer_group_for_conformer_.push_back( jj );
				//std::cerr << "lig conf " << ii << " identical to conf " << conformer_group_indices_[jj][1] << std::endl;
				conformer_unique = false;
				break;
			}
		} //loop over conformer groups

		if ( conformer_unique ) {
			conformer_group_indices_.push_back( utility::vector1< Size >() );
			conformer_group_indices_[ conformer_group_indices_.size() ].push_back( ii );
			conformer_group_for_conformer_.push_back( conformer_group_indices_.size() );
			//std::cerr << "lig conf " << ii << "detected to be unique" << std::endl;
		}
	} //loop over all lig_conformers 2->n
	TR << "Ligand conformers were split up into " << conformer_group_indices_.size() << " match-redundant groups." << std::endl;
}


/// @details
/// For results to make sense, the relevant atom indices being passed in
/// should be identical to the ones that the conformer groups were determined
/// with in the above function, although this isn't being enforced
core::Size
LigandConformerBuilder::assign_conformer_group_to_residue(
	core::conformation::Residue const & residue,
	utility::vector1< core::Size > const & relevant_atom_indices
) const
{
	core::Size num_relevant_atoms( relevant_atom_indices.size() );
	HTReal identity_ht;

	//first we have to convert the coordinates to stupid FArray2D format
	ObjexxFCL::FArray2D< numeric::Real > queryconf_coord( 3, num_relevant_atoms ), uniqueconf_coord( 3, num_relevant_atoms );
	for ( core::Size i(1); i <= num_relevant_atoms; ++i ) {
		queryconf_coord(1,i) = residue.atom( relevant_atom_indices[i] ).xyz().x();
		queryconf_coord(2,i) = residue.atom( relevant_atom_indices[i] ).xyz().y();
		queryconf_coord(3,i) = residue.atom( relevant_atom_indices[i] ).xyz().z();
	}
	lig_conformers_[ conformer_group_indices_[1][1] ]->get_global_coords_as_FArray2D( uniqueconf_coord, identity_ht, relevant_atom_indices );
	core::Real low_rms = numeric::model_quality::rms_wrapper( num_relevant_atoms, queryconf_coord, uniqueconf_coord );
	core::Size low_conf_group = 1;

	for ( Size ii = 2; ii <= conformer_group_indices_.size(); ++ii ) {
		lig_conformers_[ conformer_group_indices_[ii][1] ]->get_global_coords_as_FArray2D( uniqueconf_coord, identity_ht, relevant_atom_indices );
		core::Real this_rms = numeric::model_quality::rms_wrapper( num_relevant_atoms, queryconf_coord, uniqueconf_coord );
		if ( this_rms < low_rms ) {
			low_rms = this_rms;
			low_conf_group = ii;
		}
	}
	return low_conf_group;
}

void
LigandConformerBuilder::set_bb_grid(
	BumpGridCOP bbgrid
)
{
	parent::set_bb_grid( bbgrid );
	if ( upstream_restype_ ) {
		initialize_upstream_nonbonded_min_separation_d2();
	}
}

bool
LigandConformerBuilder::hits_potentially_incompatible() const
{
	if ( conformer_group_indices_.size() > 1 ) return true;
	return false;
}

bool
LigandConformerBuilder::compatible(
	Hit const & my_hit,
	DownstreamBuilder const & other,
	Hit const & other_hit,
	bool //first_dispatch
) const
{
	return other.compatible( other_hit, *this, my_hit, false );
}

/// @details LigandConformerBuilder checks whether the ligand conformers
/// in both hits are in the same conformer group
bool
LigandConformerBuilder::compatible(
	Hit const & my_hit,
	LigandConformerBuilder const & other,
	Hit const & other_hit,
	bool //first_dispatch
) const
{
	if ( conformer_group_for_conformer_[ my_hit.downstream_conf_id() ] == other.conformer_group_for_conformer_[ other_hit.downstream_conf_id() ] ) return true;
	return false;
}


void
LigandConformerBuilder::require_atom_to_reside_in_active_site(
	core::id::AtomID const & id
)
{
	runtime_assert( id.rsd() == 1 );
	runtime_assert( id.atomno() < downstream_restype_->natoms() );
	if ( lig_conformers_[ 1 ]->restype_id_2_collision_check_id( id.atomno() ) == 0 ) {
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
LigandConformerBuilder::atom1_radius() const
{
	debug_assert( downstream_restype_ );
	//return ZERO;
	return radii_123_[ 1 ];
}


ProbeRadius
LigandConformerBuilder::atom2_radius() const
{
	debug_assert( downstream_restype_ );
	return radii_123_[ 2 ];
}


ProbeRadius
LigandConformerBuilder::atom3_radius() const
{
	debug_assert( downstream_restype_ );
	return radii_123_[ 3 ];
}


bool
LigandConformerBuilder::atom1_belongs_in_active_site() const
{
	return ats123_reqd_in_active_site_[ 1 ];
}

bool
LigandConformerBuilder::atom2_belongs_in_active_site() const
{
	return ats123_reqd_in_active_site_[ 2 ];
}

bool
LigandConformerBuilder::atom3_belongs_in_active_site() const
{
	return ats123_reqd_in_active_site_[ 3 ];
}

LigandConformerBuilder::Real
LigandConformerBuilder::atom1_atom2_distance() const
{
	debug_assert( downstream_restype_ );
	return lig_conformers_[ 1 ]->atom1_atom2_distance();
}


LigandConformerBuilder::Real
LigandConformerBuilder::atom2_atom3_distance() const
{
	debug_assert( downstream_restype_ );
	return lig_conformers_[ 1 ]->atom2_atom3_distance();
}

/// @brief Returns an angle in degrees between the three downstream atoms.

LigandConformerBuilder::Real
LigandConformerBuilder::atom1_atom2_atom3_angle() const
{
	debug_assert( downstream_restype_ );
	return lig_conformers_[ 1 ]->atom1_atom2_atom3_angle();
}


void
LigandConformerBuilder::coordinates_from_hit(
	Hit const & hit,
	utility::vector1< AtomID > const & atom_indices,
	utility::vector1< Vector > & atom_coords
) const
{
	//std::cout << "APL DEBUG coordinates_from_hit LigandConformerBuilder" << this << std::endl;

	lig_conformers_[ hit.downstream_conf_id() ]->coordinates_from_orientation( hit.second(), atom_indices, atom_coords );
}


core::pose::PoseCOP
LigandConformerBuilder::downstream_pose_from_hit(
	Hit const & hit
) const
{

	core::conformation::Residue lig_res( *downstream_restype_, false );

	HTReal global_frame = lig_conformers_[ hit.downstream_conf_id() ]->frame_from_global_orientation( hit.second() );
	for ( Size ii = 1; ii <= lig_res.natoms(); ++ii ) {
		lig_res.set_xyz( ii, lig_conformers_[ hit.downstream_conf_id() ]->
			coordinate_in_global_frame( ii, global_frame ) );
	}

	core::pose::PoseOP pose( new core::pose::Pose() );
	pose->append_residue_by_jump( lig_res, 1 );

	//we should also set a different chain for the downstream pose
	core::pose::PDBInfoOP pdbinf( new core::pose::PDBInfo( *pose ) );
	pose->pdb_info( pdbinf );
	pose->pdb_info()->chain( 1, 'X' );

	return pose;

}

LigandConformerBuilder::Size
LigandConformerBuilder::n_possible_hits_per_at3frame() const
{
	return lig_conformers_.size();
}


void
LigandConformerBuilder::initialize_from_residue(
	Size atom1,
	Size atom2,
	Size atom3,
	Size orientation_atom1,
	Size orientation_atom2,
	Size orientation_atom3,
	core::conformation::Residue const & residue
)
{
	TR << "Initializing from residue " << residue.name() << std::endl;

	atoms_123_[ 1 ] = atom1;
	atoms_123_[ 2 ] = atom2;
	atoms_123_[ 3 ] = atom3;

	orientation_atoms_[ 1 ] = orientation_atom1;
	orientation_atoms_[ 2 ] = orientation_atom2;
	orientation_atoms_[ 3 ] = orientation_atom3;

	Size const natoms = residue.natoms();
	if ( natoms < 3 ) {
		utility_exit_with_message( "ERROR in LigandConformerBuilder: cannot build a residue with fewer than three atoms" );
	}
	downstream_restype_ = residue.type_ptr();
	atom_radii_.resize( natoms );
	atom_required_in_active_site_.resize( natoms, false );

	for ( Size ii = 1; ii <= natoms; ++ii ) {
		atom_radii_[ ii ] = probe_radius_for_atom_type( residue.atom( ii ).type() );
	}

	initialize_conformers( residue );

	for ( Size ii = 1; ii <= 3; ++ii ) {
		if ( ignore_h_collisions_ && residue.atom_type( atoms_123_[ ii ] ).element() == "H" ) {
			radii_123_[ ii ] = ZERO;
		} else {
			radii_123_[ ii ] = probe_radius_for_atom_type( residue.atom( atoms_123_[ ii ] ).type() );
		}
	}
}


void
LigandConformerBuilder::initialize_upstream_residue(
	core::chemical::ResidueTypeCOP upstream_res,
	core::scoring::etable::count_pair::CountPairFunctionCOP count_pair
)
{
	debug_assert( downstream_restype_ );
	debug_assert( upstream_res );

	upstream_restype_ = upstream_res;

	Size const natoms = downstream_restype_->natoms();

	min_sep_d2_from_upstream_atoms_.clear();
	min_sep_d2_from_upstream_atoms_.resize( natoms );

	for ( Size ii = 1; ii <= lig_conformers_[ 1 ]->n_collision_check_atoms(); ++ii ) {
		Size ii_restype_id = lig_conformers_[ 1 ]->collision_check_id_2_restype_id( ii );

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
			if ( count_pair && ( ! count_pair->count( jj, ii_restype_id, weight, path_dist ) || weight != 1.0 ) ) {
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

}

void
LigandConformerBuilder::ignore_h_collisions( bool setting )
{
	if ( downstream_restype_ != 0 ) {
		utility_exit_with_message( "ERROR: ignore_h_collisions_ must be set before the downstream restype is initialized" );
	} else {
		ignore_h_collisions_ = setting;
	}
}

void
LigandConformerBuilder::set_idealize_conformers( bool setting )
{
	idealize_conformers_ = setting;
}

void
LigandConformerBuilder::set_rmsd_unique_cutoff( core::Real setting )
{
	rmsd_unique_cutoff_ = setting;
}

std::list< Hit >
LigandConformerBuilder::build_conformer_group(
	Size const confgrp_id,
	HTReal const & atom3_frame,
	Size const scaffold_build_point_id,
	Size const upstream_conf_id,
	Size const external_geometry_id,
	core::conformation::Residue const & upstream_residue
) const
{
	debug_assert( downstream_restype_ );
	debug_assert( upstream_restype_ );
	debug_assert( bbgrid_set() );
	debug_assert( (& upstream_residue.type()) == upstream_restype_.get() );

	std::list< Hit > hitlist;

	for ( Size ii = 1; ii <= conformer_group_indices_[ confgrp_id ].size(); ++ii ) {
		Size ii_conf_id = conformer_group_indices_[ confgrp_id ][ ii ];
		//std::cout << "LigandConformerBuilder::build" << std::endl;
		//for ( Size jj = 1; jj <= 3; ++jj ) {
		//Vector const jjloc = atom3_frame * ats123_in_atom3_frame_[ jj ];
		//std::cout << "Atom D" << jj << " coordinate: " << jjloc.x() << " " << jjloc.y() << " " << jjloc.z() << std::endl;
		//}

		/// collision detection and active-site containment enforcement.
		bool ii_good( true );
		for ( Size jj = 1; jj <= lig_conformers_[ ii_conf_id ]->n_collision_check_atoms(); ++jj ) {
			Size jj_restype_id = lig_conformers_[ ii_conf_id ]->collision_check_id_2_restype_id( jj );

			Vector const jjloc = lig_conformers_[ ii_conf_id ]->coordinate_in_D3_frame( jj_restype_id, atom3_frame );
			//std::cout << "   " << downstream_restype_->atom_name( at3_frame_id_2_restype_id_[ jj ] ) << " ";
			//std::cout << jjloc.x() << "  " << jjloc.y() << " " << jjloc.z() << std::endl;
			if ( atom_radii_[ jj_restype_id ] > ZERO && bbgrid().occupied( atom_radii_[ jj_restype_id ], jjloc ) ) {
				ii_good = false;
				break;
			}
			if ( atom_required_in_active_site_[ jj_restype_id ] && ! active_site_grid().occupied( jjloc ) ) {
				ii_good = false;
				break;
			}

			for ( Size kk = 1; kk <= min_sep_d2_from_upstream_atoms_[ jj_restype_id ].size(); ++kk ) {
				if ( jjloc.distance_squared( upstream_residue.xyz( min_sep_d2_from_upstream_atoms_[ jj_restype_id ][ kk ].first ))
						< min_sep_d2_from_upstream_atoms_[ jj_restype_id ][ kk ].second ) {
					//std::cout << "collision between " << downstream_restype_->atom_name( at3_frame_id_2_restype_id_[ jj ] );
					//std::cout << " on " << downstream_restype_->name() << " with ";
					//std::cout << upstream_restype_->atom_name( min_sep_d2_from_upstream_atoms_[ jj ][ kk ].first ) << " on " << upstream_restype_->name() << std::endl;
					ii_good = false;
					break;
				}
			}
		}

		if ( ! ii_good ) continue;

		/// Check the atoms we require to be in the active site, but which are not used in
		/// collision detection
		for ( Size jj = 1; jj <= non_collision_detection_atoms_reqd_in_active_site_.size(); ++jj ) {
			Size jj_restype_id = non_collision_detection_atoms_reqd_in_active_site_[ jj ];
			Vector const jjloc = lig_conformers_[ ii_conf_id ]->coordinate_in_D3_frame( jj_restype_id, atom3_frame );

			if ( ! active_site_grid().occupied( jjloc ) ) {
				ii_good = false;
				break;
			}
		}

		if ( ! ii_good ) continue;

		Real6 global_coordinate = lig_conformers_[ ii_conf_id ]->global_orientation_from_frame3( atom3_frame );

		/// Check, if we're past the first round of hit building, that this orientation's
		/// bin in 6D is not empty.  If the bin is empty, then there is no way this
		/// ligand placement could result in a match.  Do not return this orientation as a hit.
		//std::cout << "global coordinate: ";
		//for ( Size jj = 1; jj <= 6; ++jj ) { std::cout << global_coordinate[ jj ] << " ";}
		//std::cout << std::endl;

		if ( occ_space_set() && ! occ_space().match_possible_for_hit_geometry( global_coordinate ) ) {
			continue;
		}

		/// We have a hit!
		//std::cout << "HIT!" << std::endl;

		Hit hit;
		hit.first()[ 1 ] = scaffold_build_point_id;
		hit.first()[ 2 ] = upstream_conf_id;
		hit.first()[ 3 ] = external_geometry_id;
		hit.first()[ 4 ] = ii_conf_id;
		hit.second() = global_coordinate;


		hitlist.push_back( hit ); /// new called here -- otherwise, nothing has been allocated on the heap since the Dunbrack rotamers
		//break at the first hit
		break;
	}

	return hitlist;
}


void
LigandConformerBuilder::initialize_upstream_nonbonded_min_separation_d2()
{
	runtime_assert( bbgrid_set() );
	runtime_assert( lig_conformers_[ 1 ] != 0 );

	for ( Size ii = 1; ii <= lig_conformers_[ 1 ]->n_collision_check_atoms(); ++ii ) {
		Size ii_restype_id = lig_conformers_[ 1 ]->collision_check_id_2_restype_id( ii );

		ProbeRadius ii_rad = atom_radii_[ ii_restype_id ];
		for ( Size jj = 1; jj <= min_sep_d2_from_upstream_atoms_[ ii_restype_id ].size(); ++jj ) {
			Size upstream_atom_id = min_sep_d2_from_upstream_atoms_[ ii_restype_id ][ jj ].first;
			ProbeRadius jj_rad = probe_radius_for_atom_type( upstream_restype_->atom( upstream_atom_id ).atom_type_index() );

			Real dis = bbgrid().required_separation_distance( ii_rad, jj_rad );
			min_sep_d2_from_upstream_atoms_[ ii_restype_id ][ jj ].second = dis*dis;
		}
	}

}

void
LigandConformerBuilder::initialize_conformers( core::conformation::Residue const & residue )
{
	// Retrieve the rotamer library for this ligand, and create a LigandConformer object for each
	// rotamer in the library.  Make sure that the orientation-atom geometry and the
	// D1, D2, and D3 geometry is the same for all input conformations.

	using namespace core;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::pack::dunbrack;
	using namespace core::pack::rotamers;

	SingleResidueRotamerLibraryFactory const & rotlib( *SingleResidueRotamerLibraryFactory::get_instance() );
	SingleResidueRotamerLibraryCOP res_rotlib( rotlib.get( residue.type() ) );

	if ( res_rotlib != 0 ) {
		SingleLigandRotamerLibraryCOP lig_rotlib( utility::pointer::dynamic_pointer_cast< SingleLigandRotamerLibrary const > ( res_rotlib ));

		// RM: Dependance on SingleLigandRotamerLibrary here is poor - ideally this should be generalized to work with
		// any type of SingleResidueRotamerLibrary
		if ( lig_rotlib == 0 ) {
			utility_exit_with_message( "Failed to retrieve a ligand rotamer library for "
				+ residue.name() + ". Did you mean to remove the -match::enumerate_ligand_rotamers flag from your command line?");
		}

		RotamerVector rot_vector;
		lig_rotlib->build_base_rotamers( residue.type(), rot_vector );
		Size const nligrots = rot_vector.size();
		if ( nligrots == 0 ) {
			utility_exit_with_message( "Ligand rotamer library for " + residue.name() + " has 0 rotamers." );
		}

		TR << "Found " << nligrots << " rotamers for " << residue.name() << std::endl;
		//std::cout << "APL DEBUG LigandConformerBuilder" << this << std::endl;

		conformer_group_indices_.resize( nligrots );
		conformer_group_for_conformer_.resize( nligrots );
		lig_conformers_.resize( nligrots );

		if ( idealize_conformers_ ) {

			TR << "Idealizing ligand rotamers" << std::endl;

			pose::Pose ligpose;
			ligpose.append_residue_by_jump( residue, 1 );

			for ( Size ii = 1; ii <= nligrots; ++ii ) {
				conformer_group_indices_[ ii ].resize( 1 );
				conformer_group_indices_[ ii ][ 1 ] = ii;
				conformer_group_for_conformer_[ii] = ii;
				lig_conformers_[ ii ] = protocols::toolbox::match_enzdes_util::LigandConformerOP( new toolbox::match_enzdes_util::LigandConformer );
				lig_conformers_[ ii ]->ignore_h_collisions( ignore_h_collisions_ );

				ligpose.replace_residue( 1, *rot_vector[ ii ], false );
				idealize::IdealizeMover idealizer;
				idealizer.report_CA_rmsd( false );
				idealizer.apply( ligpose );

				Real rms( 0.0 );
				for ( Size jj = 1; jj <= residue.nheavyatoms(); ++jj ) {
					rms += rot_vector[ ii ]->xyz(jj).distance_squared( ligpose.residue(1).xyz(jj) );
				}
				Real rms_this_rot( std::sqrt( rms ) / residue.nheavyatoms() );
				if ( rms_this_rot >= 0.1 )  TR.Warning << "Ligand rotamer " << ii << " has idealized RMS of " << rms_this_rot <<". Usually this number is < 0.1. Check whether ligand rotamers have the same bond lengths/bond angles as specified in the ligand .params file."  << std::endl;

				lig_conformers_[ ii ]->initialize_from_residue(
					atoms_123_[ 1 ], atoms_123_[ 2 ], atoms_123_[ 3 ],
					orientation_atoms_[ 1 ], orientation_atoms_[ 2 ], orientation_atoms_[ 3 ],
					ligpose.residue( 1 ) );
			}
			/// No error checking since we have idealized everything.  The bond lengths and angles
			/// that we need to be identical in all conformers (specifically, rotamers) so we're
			/// safe.

		} else {

			for ( Size ii = 1; ii <= nligrots; ++ii ) {
				conformer_group_indices_[ ii ].resize( 1 );
				conformer_group_indices_[ ii ][ 1 ] = ii;
				conformer_group_for_conformer_[ii] = ii;
				lig_conformers_[ ii ] = protocols::toolbox::match_enzdes_util::LigandConformerOP( new toolbox::match_enzdes_util::LigandConformer );
				lig_conformers_[ ii ]->ignore_h_collisions( ignore_h_collisions_ );
				lig_conformers_[ ii ]->initialize_from_residue(
					atoms_123_[ 1 ], atoms_123_[ 2 ], atoms_123_[ 3 ],
					orientation_atoms_[ 1 ], orientation_atoms_[ 2 ], orientation_atoms_[ 3 ],
					* rot_vector[ ii ] );
			}

			/// Error checking.  Make sure that the geometries specified in the input pdb files
			/// are in agreement between the various rotamers for atoms D1 D2 and D3 as well as
			/// for orientation atoms 1, 2 and 3.  PDB files are terribly low resolution, so we
			/// there has to be sufficient tolerance in the error thresholds.  What's a good amount?
			/// 2 thousandths of an angstrom for distances should be enough; it's harder to say
			/// with degrees.  Currently, it's coded to tolerate only a tenth of a degree disagreement.

			utility::vector1< Real > const conf1_d12_vec(    lig_conformers_[ 1 ]->atom1_atom2_distance());
			utility::vector1< Real > const conf1_d23_vec(    lig_conformers_[ 1 ]->atom2_atom3_distance());
			utility::vector1< Real > const conf1_ang123_vec( lig_conformers_[ 1 ]->atom1_atom2_atom3_angle());

			if ( conf1_d12_vec.size() > 1 || conf1_d23_vec.size() > 1 || conf1_ang123_vec.size() > 1 ) {
				utility_exit_with_message( "Somehow there is more than allowable position for the ligand rotamers; this behavior is not fully supported at the moment." );
			}

			Real const conf1_d12( conf1_d12_vec[ 1 ] );
			Real const conf1_d23( conf1_d23_vec[ 1 ] );
			Real const conf1_ang123( conf1_ang123_vec[ 1 ] );

			Real const conf1_oat_d12(    lig_conformers_[ 1 ]->oatom1_oatom2_distance());
			Real const conf1_oat_d23(    lig_conformers_[ 1 ]->oatom2_oatom3_distance());
			Real const conf1_oat_ang123( lig_conformers_[ 1 ]->oatom1_oatom2_oatom3_angle());

			Real const distance_tolerance = 2e-3;
			Real const angle_tolerance = 1e-1;

			for ( Size ii = 2; ii <= nligrots; ++ii ) {
				if ( std::abs( lig_conformers_[ ii ]->atom1_atom2_distance() - conf1_d12 ) > distance_tolerance ) {
					utility_exit_with_message( "Ligand rotamers disagree on distances between atoms "
						+ utility::trim( residue.atom_name( atoms_123_[ 1 ] ) ) + " and "
						+ utility::trim( residue.atom_name( atoms_123_[ 2 ] ) ) + ": "
						+ utility::to_string( conf1_d12 ) + " vs "
						+ utility::to_string( lig_conformers_[ ii ]->atom1_atom2_distance() )
						+ " for ligand rotamers #1 vs #" + utility::to_string( ii ) );
				}
				if ( std::abs( lig_conformers_[ ii ]->atom2_atom3_distance() - conf1_d23 ) > distance_tolerance ) {
					utility_exit_with_message( "Ligand rotamers disagree on distances between atoms "
						+ utility::trim( residue.atom_name( atoms_123_[ 2 ] ) ) + " and "
						+ utility::trim( residue.atom_name( atoms_123_[ 3 ] ) ) + ": "
						+ utility::to_string( conf1_d23 ) + " vs "
						+ utility::to_string( lig_conformers_[ ii ]->atom2_atom3_distance() )
						+ " for ligand rotamers #1 vs #" + utility::to_string( ii ) );
				}
				if ( std::abs( lig_conformers_[ ii ]->atom1_atom2_atom3_angle() - conf1_ang123 ) > angle_tolerance ) {
					utility_exit_with_message( "Ligand rotamers disagree on the angle between atoms "
						+ utility::trim( residue.atom_name( atoms_123_[ 1 ] ) ) + ", "
						+ utility::trim( residue.atom_name( atoms_123_[ 2 ] ) ) + " and "
						+ utility::trim( residue.atom_name( atoms_123_[ 3 ] ) ) + ": "
						+ utility::to_string( conf1_ang123 ) + " vs "
						+ utility::to_string( lig_conformers_[ ii ]->atom1_atom2_atom3_angle() )
						+ " for ligand rotamers #1 vs #" + utility::to_string( ii ) );
				}

				if ( std::abs( lig_conformers_[ ii ]->oatom1_oatom2_distance() - conf1_oat_d12 ) > distance_tolerance ) {
					utility_exit_with_message( "Ligand rotamers disagree on distances between atoms "
						+ utility::trim( residue.atom_name( orientation_atoms_[ 1 ] ) ) + " and "
						+ utility::trim( residue.atom_name( orientation_atoms_[ 2 ] ) ) + ": "
						+ utility::to_string( conf1_oat_d12 ) + " vs "
						+ utility::to_string( lig_conformers_[ ii ]->oatom1_oatom2_distance() )
						+ " for ligand rotamers #1 vs #" + utility::to_string( ii ) );
				}
				if ( std::abs( lig_conformers_[ ii ]->oatom2_oatom3_distance() - conf1_oat_d23 ) > distance_tolerance ) {
					utility_exit_with_message( "Ligand rotamers disagree on distances between atoms "
						+ utility::trim( residue.atom_name( orientation_atoms_[ 2 ] ) ) + " and "
						+ utility::trim( residue.atom_name( orientation_atoms_[ 3 ] ) ) + ": "
						+ utility::to_string( conf1_oat_d23 ) + " vs "
						+ utility::to_string( lig_conformers_[ ii ]->oatom2_oatom3_distance() )
						+ " for ligand rotamers #1 vs #" + utility::to_string( ii ) );
				}
				if ( std::abs( lig_conformers_[ ii ]->oatom1_oatom2_oatom3_angle() - conf1_oat_ang123 ) > angle_tolerance ) {
					utility_exit_with_message( "Ligand rotamers disagree on the angle between atoms "
						+ utility::trim( residue.atom_name( orientation_atoms_[ 1 ] ) ) + ", "
						+ utility::trim( residue.atom_name( orientation_atoms_[ 2 ] ) ) + " and "
						+ utility::trim( residue.atom_name( orientation_atoms_[ 3 ] ) ) + ": "
						+ utility::to_string( conf1_oat_ang123 ) + " vs "
						+ utility::to_string( lig_conformers_[ ii ]->oatom1_oatom2_oatom3_angle() )
						+ " for ligand rotamers #1 vs #" + utility::to_string( ii ) );
				}
			} // else ! idealize_conformers_
		}
	} else {
		TR << "No ligand rotamer library found, matching with geometry specified form .params file." << std::endl;

		conformer_group_indices_.resize( 1 );
		conformer_group_for_conformer_.resize(1);
		lig_conformers_.resize( 1 );
		conformer_group_indices_[ 1 ].resize( 1 );
		conformer_group_indices_[ 1 ][ 1 ] = 1;
		conformer_group_for_conformer_[1] = 1;
		lig_conformers_[ 1 ] = protocols::toolbox::match_enzdes_util::LigandConformerOP( new toolbox::match_enzdes_util::LigandConformer );
		lig_conformers_[ 1 ]->ignore_h_collisions( ignore_h_collisions_ );
		lig_conformers_[ 1 ]->initialize_from_residue(
			atoms_123_[ 1 ], atoms_123_[ 2 ], atoms_123_[ 3 ],
			orientation_atoms_[ 1 ], orientation_atoms_[ 2 ], orientation_atoms_[ 3 ],
			residue );

		//utility_exit_with_message( "Failed to find ligand rotamer library for " +
		//residue.name() + ". Did you mean to remove the -match::enumerate_ligand_rotamers flag from your command line?" );
	}

}

toolbox::match_enzdes_util::LigandConformerOP
LigandConformerBuilder::get_lig_conformers(core::Size conf_id) const
{
	return  lig_conformers_[ conf_id ] ;
}

//utility::vector1< utility::vector1< std::pair< core::Size, core::Real > > >
//LigandConformerBuilder::get_min_sep_d2_from_upstream_atoms() const
//{
// return min_sep_d2_from_upstream_atoms_;
//}

core::chemical::ResidueTypeCOP
LigandConformerBuilder::get_upstream_restype() const
{
	return upstream_restype_;
}

}
}
}
