// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
	numeric::xyzVector<PackstatReal> const & xyz, PackstatReal const r ) :
	id_(id),
	sphere_(sphere),
	cluster_(id),
	xyz_(xyz),
	radius_(r),
	area(0.0),
	vol(0.0),
	evdw(0.0),
	exposed_radius(0.0f),
	anb(-12345)
{}

string const CavityBall::str() const {
	ostringstream oss;
	oss << "CavityBall "
		<< id_ << " "
		<< sphere_ << " "
		<< radius_ << " "
		<< ' ' << xyz_.x() << ',' << xyz_.y() << ',' << xyz_.z()
		<< ' ';
	return oss.str();
}

// new convention is to use rosetta atom# for atom# line
// and 500 + rosetta res# for res number
// before, atom# was 10 x radius and res# was
// hole index in proteinSasa's cavity_balls_ array
string const CavityBall::hetero_atom_line( int const hetresnum, int const /*chain*/, core::Real radsub ) const
{
	std::string alpha = "ZYXWVUTSRQPONMLKJIHGFEDCBA";
	for ( int i = 1; i <= 6; ++i ) alpha += alpha; // 64 x
	return "HETATM" + I( 5, ( min( 9999, id_) ) ) + "  V   CAV "+alpha[cluster_-1]
		+ I( 4, hetresnum ) + "    "
		+ F( 8, 3, xyz_.x() ) + F( 8, 3, xyz_.y() ) + F( 8, 3, xyz_.z() )
		+ F( 6, 2, exposed_radius ) + ' ' + F( 5, 2, max(0.1,radius_-radsub) );
}

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
