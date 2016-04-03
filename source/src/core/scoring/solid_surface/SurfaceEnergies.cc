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

static THREAD_LOCAL basic::Tracer tr( "core.scoring.solid_surface.SurfaceEnergies" );

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
Energies &
SurfaceEnergies::operator = ( Energies const & rhs )
{
	debug_assert( dynamic_cast< SurfaceEnergies const * > ( & rhs ) );
	if ( this == &rhs ) return *this;

	SurfaceEnergies const & surf_rhs( static_cast< SurfaceEnergies const & > ( rhs ) );
	total_residue_ = surf_rhs.total_residue_;
	non_surface_ranges_ = surf_rhs.non_surface_ranges_;
	is_surface_ = surf_rhs.is_surface_;
	neighbor_cutoff_ = surf_rhs.neighbor_cutoff_;

	parent::operator = ( rhs );
	return *this;
}

/// @details If recurse is true, then this is the first call to same_type_as_me;
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

/// @details make a copy of this Energies( allocate actual memory for it )
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

	core::conformation::find_neighbors_naive_surface<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>( pg, neighbor_cutoff, non_surface_ranges_, is_surface_ );
}


} // solid_surface
} // scoring
} // core

