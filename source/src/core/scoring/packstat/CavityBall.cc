// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/packstat/CavityBall.cc
///
/// @brief
/// @author will sheffler


// Unit header or inline function header
#include <core/scoring/packstat/CavityBall.hh>


#include <iostream>
#include <sstream>

#include <ObjexxFCL/format.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace packstat {


using namespace std;
using namespace utility;
using namespace numeric;
using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

CavityBall::CavityBall( int const id, int const sphere,
	numeric::xyzVector<PackstatReal> const xyz, PackstatReal const r ) :
	id_(id),
	sphere_(sphere),
	cluster_(id),
	xyz_(xyz),
	radius_(r),
	area(0.0),
	vol(0.0),
	exposed_radius(0.0f),
	anb(-12345)
	// sasa_(-1234.0f)//,
	// surrounding_sasa_(N_CAVBALL_DISBINS,30,0.0f),
	//surrounding_sasa_5A_(30,0.0f),
	// heavyatom_(false),
	// num_other_balls_overlap_(-1234),
	// num_buried_other_balls_overlap_(-1234),
	// num_big_other_balls_overlap_(-1234),
	// num_big_buried_other_balls_overlap_(-1234),
	// cluster_id_(-1234),
	// cluster_(NULL)
{

	// hole_sasa_     .clear();
	// hole_sasa_     .resize(N_PR_BIN,0);
	// neighbor_count_.clear();
	// avg_bfactor_   .clear();
	// avg_occupancy_ .clear();
	// absolute_shell_rms_.clear();
	// relative_shell_rms_.clear();
	// neighbor_count_.resize(N_CAVBALL_DISBINS,-1234);
	// avg_bfactor_   .resize(N_CAVBALL_DISBINS,-1234);
	// avg_occupancy_ .resize(N_CAVBALL_DISBINS,-1234);
	// absolute_shell_rms_.resize(N_CAVBALL_DISBINS,-1234);
	// relative_shell_rms_.resize(N_CAVBALL_DISBINS,-1234);

}

// bool CavityBall::cmp( CavityBall * a, CavityBall * b ) {
//  return a->radius() > b->radius();
// }
//
// bool CavityBall::overlaps( CavityBall const *b ) const {
//  return distto(b) < -0.5;
// }

string const CavityBall::str() const {
	ostringstream oss;
	oss << "CavityBall "
		<< id_ << " "
		<< sphere_ << " "
		<< radius_ << " "
		// << sasa_ << " "
		// << I(5, num_other_balls_overlap_ )
		// << I(5, num_buried_other_balls_overlap_ )
		// << I(5, num_big_other_balls_overlap_ )
		// << I(5, num_big_buried_other_balls_overlap_ )
		<< ' ' << xyz_.x() << ',' << xyz_.y() << ',' << xyz_.z()
		<< ' ';
	// for (int ii=1; ii <= (int)big_buried_neighboring_cavity_balls_.size(); ii++) {
	//  oss<< big_buried_neighboring_cavity_balls_[ii]->id_ << ' ';
	// }
	return oss.str();
}

// new convention is to use rosetta atom# for atom# line
// and 500 + rosetta res# for res number
// before, atom# was 10 x radius and res# was
// hole index in proteinSasa's cavity_balls_ array
string const CavityBall::hetero_atom_line( int const hetresnum, int const /*chain*/, core::Real radsub ) const
{
	// string CAV = ObjexxFCL::string_of(cluster_);
	// if( CAV.size() == 1 ) {
	//  CAV = "C0"+CAV;
	// } else if( CAV.size() == 2 ) {
	//  CAV = "C"+CAV;
	// } else {
	//  CAV = CAV.substr(0,3);
	// }
	std::string alpha = "ZYXWVUTSRQPONMLKJIHGFEDCBA";
	for ( int i = 1; i <= 6; ++i ) alpha += alpha; // 64 x
	return "HETATM" + I( 5, ( min( 9999, id_) ) ) + "  V   CAV "+alpha[cluster_-1]
		+ I( 4, hetresnum ) + "    "
		+ F( 8, 3, xyz_.x() ) + F( 8, 3, xyz_.y() ) + F( 8, 3, xyz_.z() )
		+ F( 6, 2, exposed_radius ) + ' ' + F( 5, 2, max(0.1,radius_-radsub) );
	//+ F( 6, 2, evdw ) + ' ' + F( 5, 2, max(0.1,radius_-radsub) );
}


// int CavityBall::recursive_mark_hole_neighbors( vector1<CavityBall> & holes, int const cluster ) {
//  if ( cluster_id_ != -1234 ) {
//   return 0;
//  }
//  cluster_id_ = cluster;
//  int count = 1;
//  for ( int ii=1; ii<=(int)big_buried_neighboring_cavity_balls_.size(); ii++) {
//   /*cerr << "PACKING: add hole "
//        << big_buried_neighboring_cavity_balls_[ii]->id_
//        << " to cluster "
//        << cluster
//        << " base atom "
//        << id_
//        << endl; */
//   count += big_buried_neighboring_cavity_balls_[ii]->recursive_mark_hole_neighbors( holes, cluster );
//  }
//  return count;
// }


} // namespace packstat
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::packstat::CavityBall::CavityBall() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::packstat::CavityBall::save( Archive & arc ) const {
	arc( CEREAL_NVP( id_ ) ); // int
	arc( CEREAL_NVP( sphere_ ) ); // int
	arc( CEREAL_NVP( cluster_ ) ); // int
	arc( CEREAL_NVP( xyz_ ) ); // numeric::xyzVector<PackstatReal>
	arc( CEREAL_NVP( radius_ ) ); // PackstatReal
	arc( CEREAL_NVP( area ) ); // PackstatReal
	arc( CEREAL_NVP( vol ) ); // PackstatReal
	arc( CEREAL_NVP( evdw ) ); // PackstatReal
	arc( CEREAL_NVP( exposed_radius ) ); // PackstatReal
	arc( CEREAL_NVP( anb ) ); // int
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::packstat::CavityBall::load( Archive & arc ) {
	arc( id_ ); // int
	arc( sphere_ ); // int
	arc( cluster_ ); // int
	arc( xyz_ ); // numeric::xyzVector<PackstatReal>
	arc( radius_ ); // PackstatReal
	arc( area ); // PackstatReal
	arc( vol ); // PackstatReal
	arc( evdw ); // PackstatReal
	arc( exposed_radius ); // PackstatReal
	arc( anb ); // int
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::packstat::CavityBall );
#endif // SERIALIZATION
