// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/packstat/CavityBall.hh
///
/// @brief
/// @author will sheffler


#ifndef INCLUDED_core_scoring_packstat_CavityBall_hh
#define INCLUDED_core_scoring_packstat_CavityBall_hh


// Project forward headers
#include <core/scoring/packstat/types.hh>
#include <core/scoring/packstat/CavityBall.fwd.hh>

#include <numeric/xyzVector.hh>

#include <iosfwd>
#include <string>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace packstat {

class CavityBall {
public:

	CavityBall( int const id, int const sphere, XYZ const xyz, PackstatReal const r );

	std::string const hetero_atom_line( int hetresnum = 1, int chain = 0, core::Real radsub = 0.0 ) const ;
	// bool cmp( CavityBall * a, CavityBall * b );
	// inline PackstatReal  distto( CavityBall & b ) const { return distto( &b ); }
	// inline bool  touches( CavityBall & b ) const { return touches( &b ); }
	// inline PackstatReal overlap( CavityBall & b ) const { return overlap( &b ); }
	//
	// bool overlaps( CavityBall const *b) const ;
	//
	// inline PackstatReal distto( CavityBall const * b) const {
	//  return distance(this->xyz_,b->xyz_) - radius_ - b->radius();
	// }
	// inline bool touches( CavityBall * b) const {
	//  return distance_squared(this->xyz_,b->xyz_) <
	//  (radius_+b->radius())*(radius_+b->radius());
	// }
	//
	// inline PackstatReal overlap( CavityBall * b) const {
	//  PackstatReal o = std::abs(distance(this->xyz_,b->xyz_) - (radius_+b->radius()));
	//  return std::max( o/radius_, o/b->radius_ );
	// }
	// // int recursive_mark_hole_neighbors( utility::vector1<CavityBall> & holes, int const cluster );

	std::string const str() const;

	inline int    sphere()        const { return sphere_; };
	inline PackstatReal  radius()        const { return radius_; };
	// inline PackstatReal  sasa()          const { return sasa_;   };
	inline int    id()            const { return id_; };
	inline numeric::xyzVector<PackstatReal> xyz() const { return xyz_; };
	// inline void   set_cluster_id(int cid) { cluster_id_ = cid; }

	int operator<(const CavityBall &rhs) const {
		return radius() > rhs.radius();
	}
	// inline CavityBallCluster * cluster() const { return cluster_; }


	// int get_cluster_id() const { return cluster_id_; }
	//private: // make this stuf private at some point

	int id_,sphere_,cluster_;
	numeric::xyzVector<PackstatReal> xyz_;
	PackstatReal radius_;
	// PackstatReal sasa_; // if this > 0, the ball is considered pruned
	PackstatReal area,vol;
	PackstatReal evdw;
	PackstatReal exposed_radius; // largest radius that can reach ball from "outside"
	int anb;

	// utility::vector1<PackstatReal> hole_sasa_;    // these will only be computed if
	// FArray2D_float surrounding_sasa_; // compute_packing_statistics is called
	// //FArray1D_float surrounding_sasa_5A_; // withing 5A of center, like before for comparison
	// utility::vector1<int>  neighbor_count_; // # atoms within 1A, 2A, 3A...
	// utility::vector1<PackstatReal> avg_occupancy_;
	// utility::vector1<PackstatReal> avg_bfactor_;
	// utility::vector1<PackstatReal> absolute_shell_rms_;
	// utility::vector1<PackstatReal> relative_shell_rms_;
	// //int  neighbor_count_5A_;  // within 5A of center, like before
	//
	// int num_other_balls_overlap_;        // number of other "holes" we overlap with (any and all)
	// int num_buried_other_balls_overlap_; // number of non-exposed balls we overlap with
	// int num_big_other_balls_overlap_;    // number of properly sized (not too small) balls overlap
	// int num_big_buried_other_balls_overlap_; // num. properly sized and not exposed overlapers

	// utility::vector1<CavityBall*> neighboring_cavity_balls_;            // number of other "holes" we overlap with (any and all)
	// utility::vector1<CavityBall*> buried_neighboring_cavity_balls_;     // number of non-exposed balls we overlap with
	// utility::vector1<CavityBall*> big_neighboring_cavity_balls_;        // number of properly sized (not too small) balls overlap
	// utility::vector1<CavityBall*> big_buried_neighboring_cavity_balls_; // num. properly sized and not exposed overlapers
	//
	// int cluster_id_; // used for grouping holesa
	// CavityBallCluster *cluster_;

#ifdef    SERIALIZATION
public:
	/// @brief Default constructor that should only be used when deserializing
	/// %CavityBalls or containers of %CavityBalls.
	CavityBall();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

typedef utility::vector1< CavityBall > CavBalls;
typedef utility::vector1< CavityBall >::iterator CavBallIter;
typedef utility::vector1< CavityBall >::const_iterator CavBallCIter;


} // namespace packstat
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_packstat_CavityBall_HH
