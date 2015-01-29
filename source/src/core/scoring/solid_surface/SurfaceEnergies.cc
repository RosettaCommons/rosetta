// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/solid_surface/SurfaceEnergies.cc
/// @brief  SurfaceEnergies class avoids calculating energies between surface residues, and detecting their
/// neighbor relationships
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Mike Pacella (mpacella88@gmail.com)

// Unit Headers
#include <core/scoring/solid_surface/SurfaceEnergies.hh>

#include <basic/Tracer.hh>

// Project Headers
#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/scoring/ScoreFunctionInfo.hh>

// Numeric headers
#include <numeric/numeric.functions.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/string_util.hh>

static thread_local basic::Tracer tr( "core.scoring.solid_surface.SurfaceEnergies" );

namespace core {
namespace scoring {
namespace solid_surface {

SurfaceEnergies::SurfaceEnergies() :
	parent(),
	total_residue_( 0 ),
	neighbor_cutoff_( -1 )
{}


/// copy ctor -- deep copy
SurfaceEnergies::SurfaceEnergies( SurfaceEnergies const & other ) :
	parent( other ),
	total_residue_( other.total_residue_ ),
	non_surface_ranges_( other.non_surface_ranges_ ),
	is_surface_( other.is_surface_ ),
	neighbor_cutoff_( other.neighbor_cutoff_ )
	//surface_grid_( other.surface_grid_ ),
	//surf_bb_( other.surf_bb_ ),
	//surf_dim_( other.surf_dim_ )
{}

/// assignment operator -- deep copy
Energies const &
SurfaceEnergies::operator = ( Energies const & rhs )
{
debug_assert( dynamic_cast< SurfaceEnergies const * > ( & rhs ) );
	if ( this == &rhs ) return *this;

	SurfaceEnergies const & surf_rhs( static_cast< SurfaceEnergies const & > ( rhs ) );
	total_residue_ = surf_rhs.total_residue_;
	non_surface_ranges_ = surf_rhs.non_surface_ranges_;
	is_surface_ = surf_rhs.is_surface_;

	neighbor_cutoff_ = surf_rhs.neighbor_cutoff_;
	//surface_grid_ = surf_rhs.surface_grid_;
	//surf_bb_ = surf_rhs.surf_bb_;
	//surf_dim_ = surf_rhs.surf_dim_;

	parent::operator = ( rhs );

	return *this;
}

///@details If recurse is true, then this is the first call to same_type_as_me;
// determine if the other object is also a SurfaceEnergies object.  If recurse is false
// then the other object is also of type SurfaceEnergies, so return true, otherwise,
// ask the other objec to make sure it's the same type as me
bool
SurfaceEnergies::same_type_as_me( Energies const & other, bool recurse /* = true */ ) const
{
	if ( dynamic_cast< SurfaceEnergies const * > (&other) ) {
		if ( ! recurse ) {
			return true;
		} else {
			return other.same_type_as_me( *this, false );
		}
	} else {
		return false;
	}
}

SurfaceEnergies::~SurfaceEnergies() {}

///@details make a copy of this Energies( allocate actual memory for it )
EnergiesOP
SurfaceEnergies::clone() const
{
  return EnergiesOP( new SurfaceEnergies( *this ) );
}

	
/// @brief The SurfaceEnergies object needs to know how many residues are in its pose;
/// it also has to be told which residues are considered part of the surface and which
/// residues are not part of the surface.
void
SurfaceEnergies::set_total_residue( Size total_residue ) {
	total_residue_ = total_residue;
	reset_surface_residue_information();
}
	
/// @brief Wipe away all the information in the SurfacEnergies object describing which
/// residues are considered part of the surface and which residues are not part of the
/// surface.  (Afterwards, all residues are going to be considered part of the surface).
void
SurfaceEnergies::reset_surface_residue_information() {
	non_surface_ranges_.clear();
	is_surface_.resize( total_residue_ );
	std::fill( is_surface_.begin(), is_surface_.end(), true );	
}

/// @brief Tell the SurfacEnergies that the following residues are considered not
/// part of the surface.
void
SurfaceEnergies::set_residue_range_not_surface( Size seqpos_begin, Size seqpos_end ) {
	// 1. make sure the input is valid:
debug_assert( seqpos_begin <= seqpos_end );                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
debug_assert( seqpos_begin <= total_residue_ );
debug_assert( seqpos_end   <= total_residue_ );
	for ( Size ii = 1; ii <= non_surface_ranges_.size(); ++ii ) {
	debug_assert( non_surface_ranges_[ ii ].first < seqpos_end || non_surface_ranges_[ ii ].second > seqpos_begin );
	}
	non_surface_ranges_.push_back( std::make_pair( seqpos_begin, seqpos_end ));
	for ( Size ii = seqpos_begin; ii <= seqpos_end; ++ii ) {
	debug_assert( is_surface_[ ii ] );
		is_surface_[ ii ] = false;
	}
}

/// @brief Does the SurfaceEnergies object consider a particular residue to be part of the surface?
bool
SurfaceEnergies::residue_is_surface( Size seqpos ) const {
	return is_surface_[ seqpos ];
}

/// @brief determine distance cutoff threshold based on scorefxn_info_ and
/// then add edges to the PointGraph class
void
SurfaceEnergies::fill_point_graph( pose::Pose const & pose, conformation::PointGraphOP pg ) const {

	core::conformation::residue_point_graph_from_conformation( pose.conformation(), *pg );

	Distance const max_pair_radius = pose::pose_max_nbr_radius( pose );
	Distance const energy_neighbor_cutoff = 2 * max_pair_radius + get_scorefxn_info().max_atomic_interaction_distance();

	Distance const context_cutoff = max_context_neighbor_cutoff();

	Distance const neighbor_cutoff = numeric::max( energy_neighbor_cutoff, context_cutoff );

	//if ( neighbor_cutoff != neighbor_cutoff_ ) { // the surface grid is out of date!
	//	prepare_surface_grid( pg, neighbor_cutoff );
	//}

	core::conformation::find_neighbors_naive_surface<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>( pg, neighbor_cutoff, non_surface_ranges_, is_surface_ );


	//core::conformation::find_neighbors_octree_surface<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>(
		//pg, neighbor_cutoff, non_surface_ranges_, is_surface_, surface_grid_, surf_bb_, surf_dim_ );
}

/// @brief Create a new surface grid, bounding box and dimension for the surface residues.
/// This might get called if the range on the score function changes or if the surface 
/// moves.  (NOTE: currently not prepared to handle the case where the surface moves!)
/*void
SurfaceEnergies::prepare_surface_grid( conformation::PointGraphOP pg, Real neighbor_cutoff ) const
{
	neighbor_cutoff_ = neighbor_cutoff;

	core::Size const n_points( pg->num_vertices() );

	//local copy
	utility::vector1< PointPosition > points( n_points );
	for ( core::Size ii = 1; ii <= n_points; ++ii ) { points[ ii ] = pg->get_vertex( ii ).data().xyz(); }

	bool first_surface_residue_found = false;
	for ( Size ii = 1; ii <= total_residue_; ++ii ) {
		if ( ! is_surface_[ ii ] ) continue;
		if ( first_surface_residue_found ) {
			surf_bb_.add( points[ ii ] );
		} else {
			surf_bb_.set_lower( points[ ii ] );
			surf_bb_.set_upper( points[ ii ] );
			first_surface_residue_found = true;
		}
	}

	core::Size const epsilon_multiplier( 10 ); // Increase this if assert failures hit in finding a point's cube
	core::Real const epsilon( epsilon_multiplier * std::numeric_limits< core::Real >::epsilon() );
	surf_bb_.set_lower( surf_bb_.lower() - epsilon ); // Expand bounding box to assure all points get assigned cubes in it
	surf_bb_.set_upper( surf_bb_.upper() + epsilon );

	// Set cube size and dimensions within bounding box
	core::Size const side_factor( 1 ); // 1 factor => Check <= 27 adjacent cubes // 2 factor => Check <= 8 adjacent cubes
	// Might gain some speed by replacing max_residue_pair_cutoff below with the max cutoff for pairs present
	core::Real const side( side_factor * neighbor_cutoff );
debug_assert( side > core::Real( 0 ) );
	core::Real const side_inv( core::Real( 1 ) / side );
	surf_dim_ = core::conformation::CubeKey(
		core::Size( std::ceil( ( surf_bb_.upper().x() - surf_bb_.lower().x() ) * side_inv ) ),
		core::Size( std::ceil( ( surf_bb_.upper().y() - surf_bb_.lower().y() ) * side_inv ) ),
		core::Size( std::ceil( ( surf_bb_.upper().z() - surf_bb_.lower().z() ) * side_inv ) ));
	surface_grid_.clear();
	for ( Size ii = 1; ii <= total_residue_; ++ii ) {
		if ( ! is_surface_[ ii ] ) continue;
		PointPosition const pp( points[ ii ]);

		// Find the residue's cube: Cube coords are indexed from 0 to cube_dim -1
		core::conformation::CubeKey const cube_key(
			core::Size( ( pp.x() - surf_bb_.lower().x() ) * side_inv ),
			core::Size( ( pp.y() - surf_bb_.lower().y() ) * side_inv ),
			core::Size( ( pp.z() - surf_bb_.lower().z() ) * side_inv )
		);

		// Check that it is within the expanded bounding box
	debug_assert( cube_key.x() < surf_dim_.x() );
	debug_assert( cube_key.y() < surf_dim_.y() );
	debug_assert( cube_key.z() < surf_dim_.z() );

		// Add the point's position to the cube's collection
		surface_grid_[ cube_key ].push_back( ii ); // Creates the cube if it doesn't exist yet
	}

}

/// @brief Wipe everything held in the surface grid, ensuring that in the next score function
/// evaluation, a new grid will be computed.
void
SurfaceEnergies::reset_surface_grid() const
{

} */



} // solid_surface
} // scoring
} // core

