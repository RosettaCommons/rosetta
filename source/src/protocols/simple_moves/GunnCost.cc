// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief Cost computation for Gunn Moves
/// @author Oliver Lange

// Unit Headers
#include <protocols/simple_moves/GunnCost.hh>

// Package Headers

// Project Headers
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <basic/prof.hh>
#include <basic/Tracer.hh>

#include <protocols/moves/Mover.hh>


// Numeric Headers
#include <numeric/constants.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <numeric/numeric.functions.hh>

//// C++ headers
#include <cstdlib>
#include <string>

#include <cmath>

#include <core/fragment/FragData.hh>
#include <core/id/AtomID.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

using namespace core;

static thread_local basic::Tracer trGunnCost( "protocols.FragmentMover.GunnCost" );

GunnCost::GunnCost() : FragmentCost( "GunnCost", 7.0 /*cutoff*/ ), frag_cache_("GunnCost")  {}
GunnCost::GunnCost( core::Real cutoff ) :  FragmentCost( "GunnCost", cutoff ), frag_cache_("GunnCost")  {}
GunnCost::~GunnCost() {}

void GunnCost::score( core::fragment::Frame const& frame, core::pose::Pose const& pose, ScoreList &scores ) {
	runtime_assert ( frame.is_continuous() );

	scores.resize( frame.nr_frags() );

	if ( frame.length() == 1 ) {
		for ( core::Size i=1; i<=frame.nr_frags(); i++ ) {
			scores[ i ] = 5.0; //there are no scores below 2.95 .. 5.0 should give larger fragments a chance. all frags below 7 are considered
		}
		return; // jump out
	}

	GunnTuple curr_pose_gunn;
	compute_gunn( pose, frame.start(), frame.end(), curr_pose_gunn );
	for ( core::Size i=1; i<=frame.nr_frags(); i++ ) {
		GunnTuple frag_gunn;
		if ( !frag_cache_.retrieve( frame, i, frag_gunn ) ) {
			compute_gunn( frame, i, frag_gunn );
			frag_cache_.store( frame, i, frag_gunn );
		}
		scores[ i ] = score_tuple( curr_pose_gunn, frag_gunn );
	}
}

void GunnCost::compute_gunn( core::fragment::Frame const& frame, core::Size frag_num, GunnTuple &data ) {
	using namespace core::pose;

	PROF_START( basic::TEST3 );

	if ( frame.length() > various_length_poses_.size() ) {
		various_length_poses_.resize( frame.length() );
	}
	PoseOP & poseptr = various_length_poses_[ frame.length()  ];
	if ( poseptr == 0 ) {
		poseptr = PoseOP( new Pose );
		frame.fragment_as_pose( frag_num, *poseptr, chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID ) );
	} else {
		frame.fragment( frag_num ).apply( *poseptr, 1, frame.length() );
	}

	compute_gunn( *poseptr, 1, frame.length(), data );

	PROF_STOP( basic::TEST3 );
}

void GunnCost::compute_gunn( core::pose::Pose const& pose, core::Size begin, core::Size end, GunnTuple &data) {

	if ( begin == end ) return; // can't compute gunn for 1-residue fragments...

	PointPosition p1,p2,p3,temp;
	Vector x1,y1,z1; //res j
	Vector x2,y2,z2; //res k
	Vector R; // connection resj -> resk

	// setup coordinate system around res j
	chemical::ResidueType const& rtj ( pose.residue_type ( begin ) );
	id::AtomID Nj( rtj.atom_index ("N") , begin );
	id::AtomID CAj( rtj.atom_index ("CA") , begin );
	id::AtomID Cj( rtj.atom_index ("C") , begin );
	p1 = pose.xyz ( Nj );
	p2 = pose.xyz ( CAj );
	p3 = pose.xyz ( Cj );
	/*{
		Vector d1,d2,d3;
		d1 = p1-p2;
		d2 = p1-p3;
		d3 = p2-p3;
		trGunnCost.Trace << "N-CA " << d1.length() << " N-C" << d2.length() << " CA-C " << d3.length()<< std::endl;
		trGunnCost.Trace << p1[0] << " " << p1[1] << " " << p1[2] << std::endl;
		trGunnCost.Trace << p2[0] << " " << p2[1] << " " << p2[2] << std::endl;
		trGunnCost.Trace << p3[0] << " " << p3[1] << " " << p3[2] << std::endl;
	}*/
	x1 = p1 - p2;
	x1.normalize();

	cross( x1, p3 - p2, z1 );
	z1.normalize();

	cross( z1, x1, y1);

	// setup coordinate system around res k
	chemical::ResidueType const& rtk ( pose.residue_type (  end ) );
	id::AtomID Nk( rtk.atom_index ("N") , end );
	id::AtomID Ck( rtk.atom_index ("C") , end );
	id::AtomID CAk( rtk.atom_index ("CA") , end );

	p1 = pose.xyz ( Nk );
	p2 = pose.xyz ( CAk );
	p3 = pose.xyz ( Ck );
	/*
	trGunnCost.Trace  << p1[0] << " " << p1[1] << " " << p1[2] << std::endl;
	trGunnCost.Trace  << p2[0] << " " << p2[1] << " " << p2[2] << std::endl;
	trGunnCost.Trace  << p3[0] << " " << p3[1] << " " << p3[2] << std::endl;
	*/
	x2 = p1 - p2;
	x2.normalize();

	cross( x2, p3 - p2, z2 );
	z2.normalize();

	cross( z2, x2, y2 );

	//  compute inter residue vector
	R = pose.xyz( CAk ) - pose.xyz( CAj );
    	/* trGunnCost.Trace  << R[0] << " " << R[1] << " " << R[2]  << std::endl; */
	//  compute gunn quantities

	data.q6 = R.length();
	R.normalize();

	data.q1 = dot_product( x1, R );
	data.q2 = dot_product( x2, R );
	data.q3 = dot_product( x1, x2) - data.q1 * data.q2;
	data.q3 /= ( std::sqrt( (1-data.q1*data.q1) * (1-data.q2*data.q2) ) + .0001 );
	data.q3 = std::acos( numeric::sin_cos_range( data.q3 ) );
	data.q4 = std::acos( numeric::sin_cos_range( dot_product( y1, R ) / ( std::sqrt( 1-data.q1 * data.q1 ) + .0001 ) ) );
	data.q5 = std::acos( numeric::sin_cos_range( dot_product( y2, R ) / ( std::sqrt( 1-data.q2 * data.q2 ) + .0001 ) ) );

}

Real GunnCost::score_tuple( GunnTuple const& g1, GunnTuple const& g2 ) {
	using numeric::constants::d::pi_over_2;
	using numeric::constants::d::pi;
	Real c1,c2,c3,c4;
	c1 = 2.035; //magic weights
	c2 = 0.346;
	c3 = 5.72;
	c4 = 3.84;

	Real d3 ( std::abs( g1.q3 - g2.q3 ) );
	if ( d3 > pi_over_2 ) d3 = pi - d3;
	Real d4 ( std::abs( g1.q4 - g2.q4 ) );
	if ( d4 > pi_over_2 ) d4 = pi - d4;
	Real d5 ( std::abs( g1.q5 - g2.q5 ) );
	if ( d5 > pi_over_2 ) d5 = pi - d5;

	Real cost = 2.92 +
		c3 * log( 1.0+( std::abs( g1.q1 - g2.q1 ) + std::abs( g1.q2 - g2.q2 ) ) ) +
		c2 * log( 1.0+std::abs( g1.q6 - g2.q6 ) ) +
		c1 * log( 1.0+ d3 ) +
		c4 * log( 1.0+d4 + d5 );
	if ( cost < 2.95 ) cost = 100; // too similar
	return cost;
}

} //FragmentMover
} //protocols
