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


// Unit headers
#include <core/scoring/NeighborList.hh>

// Package headers
//#include <core/scoring/etable/EtableEnergy.hh>
#include <basic/options/option.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/run.OptionKeys.gen.hh>

#include <utility/vector1.hh>


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>
#include <utility/serialization/ObjexxFCL/FArray1D.srlz.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {

static THREAD_LOCAL basic::Tracer tr( "core.scoring.NeighborList" );

NeighborList::NeighborList(
	kinematics::DomainMap const & domain_map_in,
	Real const XX_cutoff,
	Real const XH_cutoff,
	Real const HH_cutoff
) :
	auto_update_( false ),
	move_tolerance_sqr_( 0.0 ),
	// add an extra (2.0 - move_tolerance_) angstroms of fudge
	wide_nblist_extension_(
	basic::options::option[basic::options::OptionKeys::run::wide_nblist_extension].user() ?
	basic::options::option[basic::options::OptionKeys::run::wide_nblist_extension]() :
	basic::options::option[basic::options::OptionKeys::run::nblist_autoupdate_wide]()
	),
	wide_move_tolerance_sqr_( 0.0 ), // this is updated when move_tolerance is set
	domain_map_( domain_map_in ), // make a copy of the domain_map
	XX_cutoff_( XX_cutoff ),
	XH_cutoff_( XH_cutoff ),
	HH_cutoff_( HH_cutoff ),
	sqrt_XX_cutoff_( std::sqrt(XX_cutoff) ),
	sqrt_XH_cutoff_( std::sqrt(XH_cutoff) ),
	sqrt_HH_cutoff_( std::sqrt(HH_cutoff) ),
	XX_cutoff_wide_( std::pow( sqrt_XX_cutoff_ + 2 * wide_nblist_extension_, 2 ) ),
	XH_cutoff_wide_( std::pow( sqrt_XH_cutoff_ + 2 * wide_nblist_extension_, 2 ) ),
	HH_cutoff_wide_( std::pow( sqrt_HH_cutoff_ + 2 * wide_nblist_extension_, 2 ) ),
	n_prepare_for_scorings_( 0 ),
	n_update_from_wide_( 0 ),
	n_full_updates_( 0 )
{
	//debug_assert( XX_cutoff_ >= XH_cutoff_ + 1e-3 &&
	// XH_cutoff_ >= HH_cutoff_ + 1e-3 ); //debug order of arguments
}

NeighborList::~NeighborList()
{
	if ( auto_update_ ) {
		tr.Debug << "Minimization stats: " << n_prepare_for_scorings_ << " score/deriv cals, ";
		tr.Debug << n_update_from_wide_ << " narrow-from-wide updates, ";
		tr.Debug << n_full_updates_ << " full updates." << std::endl;
	}
}


///////////////////////////////////////////////////////////////////////////////
void
NeighborList::check_domain_map(
	kinematics::DomainMap const & domain_map_in
) const
{
	for ( Size i=1; i<= nblist_.size(); ++i ) {
		if ( domain_map_in(i) != domain_map_(i) ) {
			utility_exit_with_message("domain_map mismatch in nblist");
		}
	}
}

void
NeighborList::set_auto_update( Distance move_tolerance )
{
	auto_update_ = true;
	move_tolerance_sqr_ = move_tolerance * move_tolerance;
	wide_move_tolerance_sqr_ = std::pow( wide_nblist_extension_ - move_tolerance, 2 );
}

void
NeighborList::disable_auto_update()
{
	auto_update_ = false;
}

/// @details instead of doing a full neighbor calculation, we're sure we can
/// update the scoring neighborlist from the wide neighborlist
///
void
NeighborList::update_from_wide_nblist( pose::Pose const & pose ) const
{
	Size const natoms_moved_more_than_move_tolerance( atoms_to_update_.size() );

	/// Move the reference coordinate for an atom that has moved more than move_tolerance_ from
	/// its reference coordinate, but do not move the reference coordinate for the other atoms
	/// that have to have their nblist_s updated.  Do this step first so that all of the reference
	/// coordinates are up-to-date by the time the next loop begins.
	for ( Size ii = 1; ii <= natoms_moved_more_than_move_tolerance; ++ii ) {
		id::AtomID const ii_atom( atoms_to_update_[ ii ] ); // weird bug if this is a reference...

		debug_assert( pose.residue( ii_atom.rsd() ).xyz( ii_atom.atomno()).distance_squared(
			wide_reference_coords_[ ii_atom.rsd() ][ ii_atom.atomno() ] )
			<= wide_move_tolerance_sqr_ );

		reference_coords_[ ii_atom.rsd() ][ ii_atom.atomno() ] = pose.residue( ii_atom.rsd() ).xyz( ii_atom.atomno() );
	}

	/// Iterate across the subset of atoms that need to be updated because they, or their
	/// neighbors have moved far enough to invalidate their nblist_ data.
	/// Note: atoms_need_updating_ grows over the course of this loop
	for ( Size ii = 1; ii <= atoms_to_update_.size(); ++ii ) {
		id::AtomID const ii_atom( atoms_to_update_[ ii ] ); // weird bug on Will's laptop if this is a reference...
		bool const ii_is_hydrogen = pose.residue( ii_atom.rsd() ).atom_is_hydrogen( ii_atom.atomno() );

		bool const ii_original_atom_moved( ii <= natoms_moved_more_than_move_tolerance );
		if ( ii_original_atom_moved ) {
			AtomNeighbors const & ii_nblist = nblist_[ ii_atom.rsd() ][ ii_atom.atomno() ];
			for ( Size jj = 1, jje = ii_nblist.size(); jj <= jje; ++jj ) {
				int const nbr_rsd = ii_nblist[ jj ].rsd();
				int const nbr_atom = ii_nblist[ jj ].atomno();
				if ( atom_needs_update_from_wide_[ nbr_rsd ][ nbr_atom ] == 0 ) {
					// add the neighbor atom to the to-do list
					atom_needs_update_from_wide_[ nbr_rsd ][ nbr_atom ] = 1;
					atoms_to_update_.push_back( id::AtomID( nbr_atom, nbr_rsd ) );
				}
			}
		}

		/// Clear out the stale nblist data
		//AtomNeighbors & ii_nblist = nblist_[ ii_atom.rsd() ][ ii_atom.atomno() ];
		//ii_nblist.clear();
		nblist_[ ii_atom.rsd() ][ ii_atom.atomno() ].clear();
		upper_nblist_[ ii_atom.rsd() ][ ii_atom.atomno() ].clear();
		intrares_upper_nblist_[ ii_atom.rsd() ][ ii_atom.atomno() ].clear();


		Vector const & ii_coord = reference_coords_[ ii_atom.rsd() ][ ii_atom.atomno() ];

		AtomNeighbors const & ii_wide_nblist = wide_nblist_[ ii_atom.rsd() ][ ii_atom.atomno() ];

		for ( Size jj = 1, jje = ii_wide_nblist.size(); jj <= jje; ++jj ) {
			int const nbr_rsd = ii_wide_nblist[ jj ].rsd();
			int const nbr_atom = ii_wide_nblist[ jj ].atomno();
			bool const nbr_is_hydrogen = pose.residue( nbr_rsd ).atom_is_hydrogen( nbr_atom );

			if ( ii_coord.distance_squared( reference_coords_[ nbr_rsd ][ nbr_atom ] ) <=
					atom_pair_cutoff( ii_is_hydrogen, nbr_is_hydrogen ) ) {
				//ii_nblist.push_back( ii_wide_nblist[ jj ] ); // copy from wide into narrow
				declare_atom_neighbor_1sided( ii_atom, id::AtomID( nbr_atom, nbr_rsd ), ii_wide_nblist[ jj ].path_dist(), ii_wide_nblist[ jj ].weight() );

				/// Neighbors of one of the original atoms to exceed the move_tolerance_sqr_ movement
				/// limitation need to have their nblist_'s updated.
				if ( ii_original_atom_moved && atom_needs_update_from_wide_[ nbr_rsd ][ nbr_atom ] == 0 ) {
					// add the neighbor atom to the to-do list
					atom_needs_update_from_wide_[ nbr_rsd ][ nbr_atom ] = 1;
					atoms_to_update_.push_back( id::AtomID( nbr_atom, nbr_rsd ) );
				}

			}
		}
	}
	/// Now that we've updated all the atoms that need updating,
	/// clear out the atom_needs_update_from_wide_ array so that the next
	/// invocation of this function starts with clean data
	for ( Size ii = 1; ii <= atoms_to_update_.size(); ++ii ) {
		id::AtomID const & ii_atom( atoms_to_update_[ ii ] );
		atom_needs_update_from_wide_[ ii_atom.rsd() ][ ii_atom.atomno() ] = 0;
	}
}

void
NeighborList::declare_atoms_neighbors( id::AtomID at1, id::AtomID at2, Size path_dist, Real weight, Real weight_func ) const
{
	nblist_[ at1.rsd() ][ at1.atomno() ].push_back( AtomNeighbor( at2.rsd(), at2.atomno(), path_dist, weight, weight_func ));
	nblist_[ at2.rsd() ][ at2.atomno() ].push_back( AtomNeighbor( at1.rsd(), at1.atomno(), path_dist, weight, weight_func ));
	if ( at1.rsd() != at2.rsd() ) {
		if ( at1.rsd() < at2.rsd() ) {
			upper_nblist_[ at1.rsd() ][ at1.atomno() ].push_back( AtomNeighbor( at2.rsd(), at2.atomno(), path_dist, weight, weight_func ));
		} else {
			upper_nblist_[ at2.rsd() ][ at2.atomno() ].push_back( AtomNeighbor( at1.rsd(), at1.atomno(), path_dist, weight, weight_func ));
		}
	} else {
		debug_assert( at1.atomno() != at2.atomno() ); // no atom-self interactions!
		if ( at1.atomno() < at2.atomno() ) {
			intrares_upper_nblist_[ at1.rsd() ][ at1.atomno() ].push_back( AtomNeighbor( at2.rsd(), at2.atomno(), path_dist, weight, weight_func ));
		} else {
			intrares_upper_nblist_[ at2.rsd() ][ at2.atomno() ].push_back( AtomNeighbor( at1.rsd(), at1.atomno(), path_dist, weight, weight_func ));
		}
	}
}

/// @details Add at2 to atom 1's neighbor lists, but don't touch at2's neighbor lists.
void
NeighborList::declare_atom_neighbor_1sided( id::AtomID at1, id::AtomID at2, Size path_dist, Real weight, Real weight_func ) const
{
	nblist_[ at1.rsd() ][ at1.atomno() ].push_back( AtomNeighbor( at2.rsd(), at2.atomno(), path_dist, weight, weight_func ));
	if ( at1.rsd() != at2.rsd() ) {
		if ( at1.rsd() < at2.rsd() ) {
			upper_nblist_[ at1.rsd() ][ at1.atomno() ].push_back( AtomNeighbor( at2.rsd(), at2.atomno(), path_dist, weight, weight_func ));
		}
	} else {
		debug_assert( at1.atomno() != at2.atomno() ); // no atom-self interactions!
		if ( at1.atomno() < at2.atomno() ) {
			intrares_upper_nblist_[ at1.rsd() ][ at1.atomno() ].push_back( AtomNeighbor( at2.rsd(), at2.atomno(), path_dist, weight, weight_func ));
		}
	}
}


} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::NeighborList::save( Archive & arc ) const {
	arc( CEREAL_NVP( domain_map_ ) ); // const kinematics::DomainMap
	arc( CEREAL_NVP( XX_cutoff_ ) ); // const DistanceSquared
	arc( CEREAL_NVP( XH_cutoff_ ) ); // const DistanceSquared
	arc( CEREAL_NVP( HH_cutoff_ ) ); // const DistanceSquared


	arc( CEREAL_NVP( auto_update_ ) ); // _Bool
	arc( CEREAL_NVP( move_tolerance_sqr_ ) ); // DistanceSquared
	arc( CEREAL_NVP( wide_nblist_extension_ ) ); // const Distance
	arc( CEREAL_NVP( wide_move_tolerance_sqr_ ) ); // DistanceSquared
	arc( CEREAL_NVP( nblist_ ) ); // utility::vector1<utility::vector1<AtomNeighbors> >
	arc( CEREAL_NVP( upper_nblist_ ) ); // utility::vector1<utility::vector1<AtomNeighbors> >
	arc( CEREAL_NVP( intrares_upper_nblist_ ) ); // utility::vector1<utility::vector1<AtomNeighbors> >
	arc( CEREAL_NVP( wide_nblist_ ) ); // utility::vector1<utility::vector1<AtomNeighbors> >
	arc( CEREAL_NVP( reference_coords_ ) ); // utility::vector1<utility::vector1<Vector> >
	arc( CEREAL_NVP( wide_reference_coords_ ) ); // utility::vector1<utility::vector1<Vector> >
	arc( CEREAL_NVP( sqrt_XX_cutoff_ ) ); // const Distance
	arc( CEREAL_NVP( sqrt_XH_cutoff_ ) ); // const Distance
	arc( CEREAL_NVP( sqrt_HH_cutoff_ ) ); // const Distance
	arc( CEREAL_NVP( XX_cutoff_wide_ ) ); // const DistanceSquared
	arc( CEREAL_NVP( XH_cutoff_wide_ ) ); // const DistanceSquared
	arc( CEREAL_NVP( HH_cutoff_wide_ ) ); // const DistanceSquared
	arc( CEREAL_NVP( atom_needs_update_from_wide_ ) ); // utility::vector1<utility::vector1<Size> >
	arc( CEREAL_NVP( atoms_to_update_ ) ); // utility::vector1<id::AtomID>
	arc( CEREAL_NVP( n_prepare_for_scorings_ ) ); // Size
	arc( CEREAL_NVP( n_update_from_wide_ ) ); // Size
	arc( CEREAL_NVP( n_full_updates_ ) ); // Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::NeighborList::load_and_construct( Archive & arc, cereal::construct< core::scoring::NeighborList > & construct ) {
	kinematics::DomainMap domain_map; arc( domain_map );
	DistanceSquared XX_cutoff; arc( XX_cutoff );
	DistanceSquared XH_cutoff; arc( XH_cutoff );
	DistanceSquared HH_cutoff; arc( HH_cutoff );

	construct( domain_map, XX_cutoff, XH_cutoff, HH_cutoff );
	// The line below says to the serialization validator that four variables
	// do not need to be explicitly initialized.
	// EXEMPT domain_map_ XX_cutoff_ XH_cutoff_ HH_cutoff_

	arc( construct->auto_update_ ); // _Bool
	arc( construct->move_tolerance_sqr_ ); // DistanceSquared
	arc( const_cast< Real & > ( construct->wide_nblist_extension_ ) ); // const Distance
	arc( construct->wide_move_tolerance_sqr_ ); // DistanceSquared
	arc( construct->nblist_ ); // utility::vector1<utility::vector1<AtomNeighbors> >
	arc( construct->upper_nblist_ ); // utility::vector1<utility::vector1<AtomNeighbors> >
	arc( construct->intrares_upper_nblist_ ); // utility::vector1<utility::vector1<AtomNeighbors> >
	arc( construct->wide_nblist_ ); // utility::vector1<utility::vector1<AtomNeighbors> >
	arc( construct->reference_coords_ ); // utility::vector1<utility::vector1<Vector> >
	arc( construct->wide_reference_coords_ ); // utility::vector1<utility::vector1<Vector> >

	arc( const_cast< Real & > ( construct->sqrt_XX_cutoff_ ) ); // const Distance
	arc( const_cast< Real & > ( construct->sqrt_XH_cutoff_ ) ); // const Distance
	arc( const_cast< Real & > ( construct->sqrt_HH_cutoff_ ) ); // const Distance
	arc( const_cast< Real & > ( construct->XX_cutoff_wide_ ) ); // const DistanceSquared
	arc( const_cast< Real & > ( construct->XH_cutoff_wide_ ) ); // const DistanceSquared
	arc( const_cast< Real & > ( construct->HH_cutoff_wide_ ) ); // const DistanceSquared
	arc( construct->atom_needs_update_from_wide_ ); // utility::vector1<utility::vector1<Size> >
	arc( construct->atoms_to_update_ ); // utility::vector1<id::AtomID>
	arc( construct->n_prepare_for_scorings_ ); // Size
	arc( construct->n_update_from_wide_ ); // Size
	arc( construct->n_full_updates_ ); // Size
}
SAVE_AND_LOAD_AND_CONSTRUCT_SERIALIZABLE( core::scoring::NeighborList );
CEREAL_REGISTER_TYPE( core::scoring::NeighborList )


/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::AtomNeighbor::save( Archive & arc ) const {
	arc( CEREAL_NVP( rsd_ ) ); // int
	arc( CEREAL_NVP( atomno_ ) ); // int
	arc( CEREAL_NVP( path_dist_ ) ); // Size
	arc( CEREAL_NVP( weight_ ) ); // Real
	arc( CEREAL_NVP( weight_func_ ) ); // Real
	arc( CEREAL_NVP( temp1_ ) ); // Real
	arc( CEREAL_NVP( temp2_ ) ); // Real
	arc( CEREAL_NVP( temp3_ ) ); // Real
	arc( CEREAL_NVP( temp4_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::AtomNeighbor::load( Archive & arc ) {
	arc( rsd_ ); // int
	arc( atomno_ ); // int
	arc( path_dist_ ); // Size
	arc( weight_ ); // Real
	arc( weight_func_ ); // Real
	arc( temp1_ ); // Real
	arc( temp2_ ); // Real
	arc( temp3_ ); // Real
	arc( temp4_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::AtomNeighbor );
CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_NeighborList )
#endif // SERIALIZATION
