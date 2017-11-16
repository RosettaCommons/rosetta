// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fldsgn/potentials/sspot/NatbiasHelicesSheetPotential.cc
/// @brief native-biased centroid score for helices on sheets
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit header
#include <protocols/fldsgn/potentials/sspot/NatbiasHelicesSheetPotential.hh>

// Package headers
#include <protocols/fldsgn/topology/HelixPairing.hh>
#include <protocols/fldsgn/topology/HSSTriplet.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <basic/Tracer.hh>

// numeric
#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>

#include <utility/vector1.hh>


static basic::Tracer TR( "protocols.fldsgn.potentials.sspot.NatbiasHelicesSheetPotential", basic::t_info );

namespace protocols {
namespace fldsgn {
namespace potentials {
namespace sspot {


/// @Brief default constructor
NatbiasHelicesSheetPotential::NatbiasHelicesSheetPotential():
	hss3set_( /* NULL */ ),
	hpairset_( /* NULL */ )
{
	set_params();
}

/// @brief value constructor
NatbiasHelicesSheetPotential::NatbiasHelicesSheetPotential( HSSTripletSetOP const hss3set ):
	hss3set_( hss3set ),
	hpairset_( /* NULL */ )
{
	set_params();
}

/// @brief value constructor
NatbiasHelicesSheetPotential::NatbiasHelicesSheetPotential( HSSTripletSetOP const hss3set, HelixPairingSetOP const hpairset ):
	hss3set_( hss3set ),
	hpairset_( hpairset )
{
	set_params();
}

/// @brief value constructor
NatbiasHelicesSheetPotential::NatbiasHelicesSheetPotential( NatbiasHelicesSheetPotential const & src ):
	ReferenceCount(),
	hss3set_( src.hss3set_ ),
	hpairset_( src.hpairset_ ),
	hs_dist_wts_( src.hs_dist_wts_ ),
	hs_dist_( src.hs_dist_ ),
	hs_dist_sigma2_( src.hs_dist_sigma2_ ),
	hsheet_dist_repulsive_( src.hsheet_dist_repulsive_ ),
	hs_angle_wts_( src.hs_angle_wts_ ),
	hs_angle_( src.hs_angle_ ),
	hs_angle_sigma2_( src.hs_angle_sigma2_ ),
	hh_align_angle_wts_( src.hh_align_angle_wts_ ),
	hh_align_angle_( src.hh_align_angle_ ),
	hh_align_angle_sigma2_( src.hh_align_angle_sigma2_ )
{}

/// @brief default destructor
NatbiasHelicesSheetPotential::~NatbiasHelicesSheetPotential()
{}

/// @brief set parameters
void
NatbiasHelicesSheetPotential::set_params()
{
	// hs_dist_wts_ + hs_angle_wts_ supposed to be 1.0
	hs_dist_wts_ =  0.5;
	hs_dist_ = 13.0;
	hs_dist_sigma2_ = 6.5;

	hs_angle_wts_ = 0.5;
	hs_angle_ = -12.5;
	hs_angle_sigma2_ = 10.0;

	hsheet_dist_repulsive_ = 6.5;

	hh_align_angle_wts_ = 1.0;
	hh_align_angle_ = 15.0;
	hh_align_angle_sigma2_ = 15.0;
}

/// @brif set HSSTrirpletSet
void
NatbiasHelicesSheetPotential::hss_triplet_set( HSSTripletSetOP const hss3set )
{
	hss3set_ = hss3set;
}

/// @brief set HelixPairingSet
void
NatbiasHelicesSheetPotential::hpairset( HelixPairingSetOP const hpairset )
{
	hpairset_ = hpairset;
}

/// @brief set dist parameters for helix-sheet interaction
void
NatbiasHelicesSheetPotential::set_atrdist_params_helix_strands(
	Real const hs_dist_wts,
	Real const hs_dist,
	Real const hs_dist_sigma2 )
{
	hs_dist_wts_ = hs_dist_wts;
	hs_dist_ = hs_dist;
	hs_dist_sigma2_ = hs_dist_sigma2;
}

/// @brief set dist parameters for helix-sheet interaction
void
NatbiasHelicesSheetPotential::set_repldist_params_helix_sheet( Real const hsheet_dist_repulsive )
{
	hsheet_dist_repulsive_ = hsheet_dist_repulsive;
}

/// @brief set angle parameters for helix-sheet interaction
void
NatbiasHelicesSheetPotential::set_angle_params_helix_sheet(
	Real const hs_angle_wts,
	Real const hs_angle,
	Real const hs_angle_sigma2 )
{
	hs_angle_wts_ = hs_angle_wts;
	hs_angle_ = hs_angle;
	hs_angle_sigma2_ = hs_angle_sigma2;
}

/// @brief set angle parameters for helices on sheet
void
NatbiasHelicesSheetPotential::set_angle_params_helices_on_sheet(
	Real const hh_align_angle_wts,
	Real const hh_align_angle,
	Real const hh_align_angle_sigma2 )
{
	hh_align_angle_wts_ = hh_align_angle_wts;
	hh_align_angle_ = hh_align_angle;
	hh_align_angle_sigma2_ = hh_align_angle_sigma2;
}

/// @brief
void
NatbiasHelicesSheetPotential::show_params() const
{

	TR << "Distance pot of helix & strands: wts, atr_dist, sigma2 : "
		<< hs_dist_wts_ << " " << hs_dist_ << " " << hs_dist_sigma2_ << std::endl;

	TR << "Distance pot of helix & sheet: " << hsheet_dist_repulsive_ << std::endl;

	TR << "Angle pot of helix & strands: wts, angle, sigma2 : "
		<< hs_angle_wts_ << " " << hs_angle_ << " " << hs_angle_sigma2_ << std::endl;

	TR << "Angle pot of helix-helix projected onto sheet: wts, angle, sigma2 : "
		<< hh_align_angle_wts_ << " " << hh_align_angle_ << " "<< hh_align_angle_sigma2_ << std::endl;
}

/// @brief calc score
void
NatbiasHelicesSheetPotential::score( SS_Info2_COP const ss_info, Real & hh_score, Real & hs_score ) const
{
	using protocols::fldsgn::topology::Helix;
	using protocols::fldsgn::topology::Helices;
	using protocols::fldsgn::topology::Strand;
	using protocols::fldsgn::topology::Strands;
	using protocols::fldsgn::topology::HelixPairing;
	using protocols::fldsgn::topology::HelixPairings;
	using protocols::fldsgn::topology::HSSTriplets;

	// parameters for score calculation
	Real hs_dist_wts = 0.5*hs_dist_wts_;
	Real hs_dist_sigma2 = 2*hs_dist_sigma2_;
	Real hs_angle_sigma2 = 2*hs_angle_sigma2_;
	Real judge_hs_close = -0.5;

	// set helices and strands from ss_info
	Helices const & helices( ss_info->helices() );
	Strands const & strands( ss_info->strands() );

	// a helix is close to a sheet ?
	// If true, the helix-pair score on sheet for the helix will be calculated
	utility::vector1< bool > is_hs_close( ss_info->helices().size(), false );

	Size num( 0 );
	hs_scores_.resize( hss3set_->size() );

	// interaction between helix and sheet
	for ( HSSConstIterator it=hss3set_->begin(), ite=hss3set_->end(); it !=ite; ++it ) {

		num++;
		hs_scores_[ num ] = 0.0;

		HSSTripletOP hssop( *it );

		runtime_assert( helices.size() >= hssop->helix() );
		runtime_assert( strands.size() >= hssop->strand1() );
		runtime_assert( strands.size() >= hssop->strand2() );

		hssop->calc_geometry( ss_info );

		Real hsheet_dist = hssop->hsheet_dist();
		Real hs1_dist = hssop->hs1_dist();
		Real hs2_dist = hssop->hs2_dist();

		// repulsive between helix and sheet
		if ( hsheet_dist < hsheet_dist_repulsive_ ) {

			hs_scores_[ num ] += 10.0;
			TR.Debug << "hseet_dist=" << hsheet_dist << std::endl;

		} else {

			Real dist_score1( 0.0 ), dist_score2( 0.0 ), angle_score( 0.0 );

			// distance energy between helix and 1st strand
			if ( hs1_dist <= hs_dist_ ) {
				dist_score1 = -1.0;
			} else {
				Real r = numeric::square( hs1_dist - hs_dist_ )/hs_dist_sigma2;
				dist_score1 = -exp( -r );
			}

			// distance energy between helix and 2nd strand
			if ( hs2_dist <= hs_dist_ ) {
				dist_score2 = -1.0;
			} else {
				Real r = numeric::square( hs2_dist - hs_dist_ )/hs_dist_sigma2;
				dist_score2 = -exp( -r );
			}

			TR.Debug << "HS_dist_score: "
				<< hssop->helix() << "-" << hssop->strand1() << "," << hssop->strand2() << " "
				<< hsheet_dist << " " << hs1_dist << " " << hs2_dist << " "
				<< dist_score1 << " " << dist_score2 << std::endl;

			// angle between helix and sheet
			if ( dist_score1 <= judge_hs_close && dist_score2 <= judge_hs_close ) {

				is_hs_close[ hssop->helix() ] = true;

				Real hs_angle = hssop->hs_angle();

				///  90 < hs_angle < 180 has to have penalty !! this functionality need to be implemented.
				if ( hs_angle >= hs_angle_ ) {
					angle_score = -1.0;
				} else {
					Real r = numeric::square( hs_angle - hs_angle_ )/hs_angle_sigma2;
					angle_score = -exp( -r );
				}

				TR.Debug << "HS_angle_score: "
					<< hssop->helix() << "-" << hssop->strand1() << "," << hssop->strand2() << " "
					<< hs_angle << " " << angle_score << " " << std::endl;
			}

			// add to hs_score
			hs_scores_[ num ] = hs_dist_wts * ( dist_score1 + dist_score2 ) + hs_angle_wts_ * angle_score;
		}

		hs_score += hs_scores_[ num ];

	}

	Real hh_align_angle_sigma2 = 2*hh_align_angle_sigma2_;

	// interaction of helix pair on sheet
	if ( hpairset_ ) {

		hh_scores_.resize( hpairset_->size() );

		Size num( 0 );
		HelixPairings const & hpairs = hpairset_->helix_pairings();
		for ( HelixPairings::const_iterator it=hpairs.begin(), ite=hpairs.end(); it != ite; ++it ) {

			num++;
			hh_scores_[ num ] = 0.0;

			HelixPairing const & hpair( **it );

			runtime_assert( helices.size() >= hpair.h1() && helices.size() >= hpair.h2() );

			// only if both helices are close to sheet, helix pairing score is calculated
			if ( !is_hs_close[ hpair.h1() ] || !is_hs_close[ hpair.h2() ] ) continue;

			Helix const & h1 = *helices[ hpair.h1() ];
			Helix const & h2 = *helices[ hpair.h2() ];

			HSSTripletOP const hssop1 = get_hssop( hpair.h1() );
			HSSTripletOP const hssop2 = get_hssop( hpair.h2() );

			// make sure that HSSTripletOP is defined
			runtime_assert( hssop1 && hssop2 );

			// get helix-pair orientation
			Real flip( 1.0 );

			Size const s1 = hssop1->strand1();
			Size const s2 = hssop2->strand2();

			Real ag = numeric::conversions::degrees( angle_of( strands[ s1 ]->orient(), strands[ s2 ]->orient() ) );
			if ( ag > 90.0 ) flip = -1.0;
			Vector const v1 = strands[ s2 ]->mid_pos() - strands[ s1 ]->mid_pos();
			Vector const v2 = strands[ s1 ]->orient() + flip*strands[ s2 ]->orient();
			Vector const v3 = h1.orient().project_parallel( v1 ) + h1.orient().project_parallel( v2 );
			Vector const v4 = h2.orient().project_parallel( v1 ) + h2.orient().project_parallel( v2 );

			if ( hpair.orient() == 'A' ) flip = -1.0;
			Real angle = numeric::conversions::degrees( angle_of( v3, flip*v4 ) );

			Real score;
			if ( angle <= hh_align_angle_ ) {
				score = -hh_align_angle_wts_;
			} else {
				Real r = numeric::square( angle - hh_align_angle_ )/hh_align_angle_sigma2;
				score = -hh_align_angle_wts_*exp( -r );
			}
			hh_scores_[ num ] = score;
			hh_score += hh_scores_[ num ];

			TR.Debug << "HH_align_score : " << hpair.h1() << "-" << hpair.h2() << " " << angle << " " << score << std::endl;

		} // for( HelixPairings )
	} // if( hpairset_ )


}  // score

NatbiasHelicesSheetPotential::HSSTripletOP
NatbiasHelicesSheetPotential::get_hssop( Size const helix_id ) const
{
	using protocols::fldsgn::topology::HSSTriplets;
	HSSTriplets const triplets = hss3set_->hss_triplets( helix_id );
	if ( triplets.empty() ) return HSSTripletOP();

	// TL: the code above depends on there only being one triplet for each helix
	//     however, there is no reason why a user might now want to specify more than one.
	//     If a user wants to use this potential with helices in more than one triplet, the code
	//     above needs to be modified.  As a stopgap, we throw an error to prevent unwanted/unexpected
	//     behavior
	if ( triplets.size() > 1 ) {
		std::stringstream msg;
		msg << "Helix " << helix_id << " is present in more than one HSS triplet. This is not yet supported by "
			<< "NatbiasHelicesSheetPotential -- please update protocols/fldsgn/sspot/NatbiasHelicesSheetPotential.cc"
			<< std::endl;
		utility_exit_with_message( msg.str() );
	}

	return *triplets.begin();
}


} // ns sspot
} // ns potentials
} // ns fldsgn
} // ns protocols

