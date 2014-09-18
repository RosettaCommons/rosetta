// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   ./src/protocols/fldsgn/potentials/sspot/NatbiasHelixPairPotential.cc
/// @brief  calss for helix-pairing potential
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit header
#include <protocols/fldsgn/potentials/sspot/NatbiasHelixPairPotential.hh>

// Package headers

// Project headers
// AUTO-REMOVED #include <core/pose/Pose.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/fldsgn/topology/HelixPairing.hh>
#include <numeric/numeric.functions.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


static thread_local basic::Tracer TR( "protocols.fldsgn.potentials.sspot.NatbiasHelixPairPotential", basic::t_info );

namespace protocols {
namespace fldsgn {
namespace potentials {
namespace sspot {

/// @Brief default constructor
NatbiasHelixPairPotential::NatbiasHelixPairPotential():
	hpairset_( NULL )
{
	set_params();
}


/// @brief value constructor
NatbiasHelixPairPotential::NatbiasHelixPairPotential( HelixPairingSetOP const hpairset ):
	hpairset_( hpairset )
{
	set_params();
}

/// @brief set parmeters
void
NatbiasHelixPairPotential::set_params()
{
	bend_angle_ = 20.0;

	dist_wts_ = 0.5;
	mid_dist_ = 15.0;
	dist_sigma2_ = 7.5;

	angle_wts_ = 0.5;
	cross_angle_ = 45.0;
	angle_sigma2_ = 20.0;
}

/// @brief copy constructor
NatbiasHelixPairPotential::NatbiasHelixPairPotential( NatbiasHelixPairPotential const & src ):
	ReferenceCount(),
	bend_angle_( src.bend_angle_ ),
	dist_wts_( src.dist_wts_ ),
	mid_dist_( src.mid_dist_ ),
	dist_sigma2_( src.dist_sigma2_ ),
	angle_wts_( src.angle_wts_ ),
	cross_angle_( src.cross_angle_ ),
	angle_sigma2_( src.angle_sigma2_ ),
	hpairset_( src.hpairset_ )
{}


/// @brief default destructor
NatbiasHelixPairPotential::~NatbiasHelixPairPotential()
{}

/// @brief set HelixPairingSet
void
NatbiasHelixPairPotential::set_hpairset( HelixPairingSetOP const hpairset )
{
	hpairset_ = hpairset;
}

/// @brief set parameters for distance score between mid points of helices
void
NatbiasHelixPairPotential::set_params_distance_pot( Real w, Real d, Real s )
{
	dist_wts_ = w;
	mid_dist_ = d;
	dist_sigma2_ = s;
}

/// @brief set parameters for angle score of helix pair
void
NatbiasHelixPairPotential::set_params_angle_pot( Real w, Real d, Real s )
{
	angle_wts_ = w;
	cross_angle_ = d;
	angle_sigma2_ = s;
}

/// @brief show parameters
void
NatbiasHelixPairPotential::show_params() const
{
	TR << "HH_dist, wts, dist, sigma2 " << dist_wts_ << " " << mid_dist_ << " " << dist_sigma2_ << std::endl;
	TR << "HH_angle, wts, dist, sigma2 " << angle_wts_ << " " << cross_angle_ << " " << angle_sigma2_ << std::endl;
}

/// @brief
void
NatbiasHelixPairPotential::show( Pose const & pose ) const
{
	if( ! hpairset_ ) {
		TR << "No helix pairings to be calculated. ";
		return;
	}

	Size nhpairs = hpairset_->helix_pairings().size();
	SS_Info2_OP ssinfo = new SS_Info2( pose );
	utility::vector1< Real > hh_scores( nhpairs, 0.0 );

	Real hh_score( 0.0 );
	score( ssinfo, hh_score );

	Size num( 0 );
	HelixPairings const & hpairs = hpairset_->helix_pairings();
	TR << "name distance cross_angle align_angle score " << std::endl;
	for( HelixPairings::const_iterator it=hpairs.begin(), ite=hpairs.end() ; it != ite; ++it ) {
		num++;
		HelixPairing const & hpair( **it );
		TR << hpair << " " << hh_scores_[ num ] << std::endl;
	}

} // show


/// @brief
void
NatbiasHelixPairPotential::score(	SS_Info2_COP const ss_info, Real & hh_score ) const
{
	using protocols::fldsgn::topology::Helix;
	using protocols::fldsgn::topology::Helices;
	using protocols::fldsgn::topology::Strands;

	if( ! hpairset_ ) {
		return;
	}

	// calc geometry
	hpairset_->calc_geometry( ss_info );

	// helices for checking bend of them
	Helices const & helices( ss_info->helices() );

	// take helix pairings
	HelixPairings const & hpairs = hpairset_->helix_pairings();

	Real asigma = 2*angle_sigma2_;
	Real dsigma = 2*dist_sigma2_;

	Size num( 0 );

	hh_scores_.resize( hpairset_->helix_pairings().size() );
	for( HelixPairings::const_iterator it=hpairs.begin(), ite=hpairs.end(); it != ite; ++it ) {

		num++;
		hh_scores_[ num ] = 0.0;
		HelixPairing const & hpair( **it );

		// check bending of helix and give high penalty if helix bend is big
		if ( helices[ hpair.h1() ]->bend() > bend_angle_ ) hh_scores_[ num ] += 10.0;
		if ( helices[ hpair.h2() ]->bend() > bend_angle_ ) hh_scores_[ num ] += 10.0;

		if ( hpair.dist() <= mid_dist_ ) {

			// Edist = -wd if ( d < d0 )
			hh_scores_[ num ] += -dist_wts_;

			if ( hpair.cross_angle() <= cross_angle_ ) {
				// Eangle = -wa if( a < a0 )
				hh_scores_[ num ] += -angle_wts_;
			} else {
				//  Eangle = -wa * exp( -( c-c0 )**2/(2*sa) ) if ( a > a0 )
				Real r = numeric::square( hpair.cross_angle() - cross_angle_ )/asigma;
				Real s = - angle_wts_*exp( -r );
				hh_scores_[ num ] += s;
			}

		} else {
			// Edist = -wd * exp( -( d-d0 )**2/(2*sd) ) if ( d > d0 )
			Real r = numeric::square( hpair.dist() - mid_dist_ )/dsigma;
			Real s = - dist_wts_*exp( -r );
			hh_scores_[ num ] += s;
		}

		hh_score += hh_scores_[ num ];

		TR.Debug << hpair << " " << hh_scores_[ num ] << std::endl;

	} // for ( HelixPairings )

}  // score


} // ns sspot
} // ns potentials
}	// ns fldsgn
}	// ns protocols
