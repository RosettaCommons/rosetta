// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit header
#include <protocols/fldsgn/potentials/sspot/SSPairPotential.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/fldsgn/topology/BB_Pos.hh>
#include <protocols/fldsgn/topology/DimerPairing.hh>
#include <protocols/fldsgn/potentials/sspot/util.hh>

#include <basic/database/open.hh>
#include <basic/Tracer.hh>

// Numeric
#include <numeric/conversions.hh>

// utility
#include <utility/io/izstream.hh>

// C++ headers
#include <cmath>
#include <iostream>

//#include <ObjexxFCL/FArray1Da.hh>
#include <ObjexxFCL/format.hh>

#include <utility/vector1.hh>
#include <ObjexxFCL/FArray3D.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>
#include <core/scoring/SecondaryStructurePotential.hh>


static THREAD_LOCAL basic::Tracer TR( "protocols.fldsgn.potentials.sspot.SSPairPotential", basic::t_info );

namespace protocols {
namespace fldsgn {
namespace potentials {
namespace sspot {

/// @brief default constructor
SSPairPotential::SSPairPotential():
	strand_dist_cutoff_( 6.5 ),
	dimer_seqsep_cutoff_( 11 ),  // Size or 6 ?
	lowstrand_( 1 ),
	phithetascore_( 2, 3, 36, 36 ),
	dotscore_( 6 ),
	rsigma_dot_( 12, 18, 2, 2, rsigma_dot_initializer )
{
	load_phi_theta_bins();
	load_dotscore_bins();
}

/// @brief destructor
SSPairPotential::~SSPairPotential()
{}

/// @brief return score of phitheta
core::Real
SSPairPotential::calc_phithetascore( Size const strand_seqsep, Real const phi, Real const theta ) const
{
	Size istrand_seqsep;
	if ( strand_seqsep >= 2 && strand_seqsep <= 10 ) {
		istrand_seqsep = 2;
	} else {
		if ( strand_seqsep > 10 ) {
			istrand_seqsep = 3;
		} else {
			istrand_seqsep = 1;
		}
	}
	Size iphi = static_cast< Size >( 1 + ( phi + 180.0 )/10 );
	if ( iphi > 36 ) {
		iphi = 36;
	} else if ( iphi < 1 ) {
		iphi = 1;
	}
	Size itheta = static_cast< Size >( 1 + ( theta/5 ) );
	if ( itheta > 36 ) {
		itheta = 36;
	} else if ( itheta < 1 ) {
		itheta = 1;
	}
	return phithetascore_( 2, int( istrand_seqsep ), int( iphi ), int( itheta ) );
}

/// @brief return dot score
core::Real
SSPairPotential::calc_dotscore( Real const dpall ) const
{
	int idot;
	if ( dpall > 0.33 && dpall <= 0.66 ) {
		idot = 2;
	} else if ( dpall > 0.66 && dpall <= 1.00 ) {
		idot = 3;
	} else if ( dpall > 1.00 && dpall <= 1.33 ) {
		idot = 4;
	} else if ( dpall > 1.33 && dpall <= 1.60 ) {
		idot = 5;
	} else if ( dpall > 1.60 && dpall <= 1.80 ) {
		idot = 6;
	} else if ( dpall > 1.80 && dpall <= 2.00 ) {
		idot = 5;
	} else {
		idot = 1;
	}
	return dotscore_( idot );
}

/// @brief return rsigma score
core::Real
SSPairPotential::calc_rsigmascore( Real sig, Real dist, Size const sign1, Size const sign2 ) const
{
	if ( sig > 179.0 ) {
		sig = 179.0;
	} else if ( sig < 0.0 ) {
		sig = 0.0;
	}
	Size isig = static_cast< Size >( sig / 10 ) + 1;
	if ( dist > 6.4 ) {
		dist = 6.4;
	} else if ( dist < 3.5 ) {
		dist = 3.5;
	}
	Size idist = static_cast< Size >( ( dist - 3.5 ) / 0.25 ) + 1;

	// FIX THIS !!!!!!!!!!!!!!
	// The definition of dimer signs (sign1,sign2) appears inverted (1 should be 2, vice versa).
	Real tempscore_rsigma = rsigma_dot_( int( idist ), int( isig ), 3 - sign1, 3 - sign2 );

	// Modify sigma potential to no longer give an rsigma bonus to strands in wrong register.
	if ( sign1 == 1 && sign2 == 1 && sig < 110. && sig > 70.  ) tempscore_rsigma = 0.0;
	if ( sign1 == 1 && sign2 == 2 && (sig < 75. || sig > 95.) ) tempscore_rsigma = 0.0;
	if ( sign1 == 2 && sign2 == 1 && (sig < 90. || sig > 110.) ) tempscore_rsigma = 0.0;
	if ( sign1 == 2 && sign2 == 2 && sig < 120. && sig > 80.  ) tempscore_rsigma = 0.0;

	return tempscore_rsigma;
}

/// @brief calculate sum of dot product of the co vectors of strand dimers ss1 and ss2
/// @brief with the vector connecting the midpoints of the dimer vectors (vdist)
/// @brief also determine return the sign of the dot products for each dimer
/// @brief to determine which direction the CO groups point
void
SSPairPotential::pair_dp(
	Size const & ss1,
	Size const & ss2,
	BB_Pos const & bb_pos,
	Real & dp,
	Vector const & mid_vector,
	Size & sign1,
	Size & sign2
) const
{
	// length of C=O bond
	static Real const dist_co_inv = { 1.0f / 1.231015f };

	// unit vector of mid_vector
	Vector const u_midvec( mid_vector.normalized_or_zero() );


	Real dp1( 0.0 );
	Real sdp1( 0.0 );
	for ( Size i=ss1; i<=ss1+1; ++i ) {
		Vector temp;
		if ( i == ss1+1 ) {
			temp = dist_co_inv * ( bb_pos.C(i) - bb_pos.O(i) );
		} else {
			temp = dist_co_inv * ( bb_pos.O(i) - bb_pos.C(i) );
		}
		Real const tempdot = temp.dot( u_midvec );
		dp1 += std::abs(tempdot);
		sdp1 += tempdot;
	}
	dp1 *= 0.5;


	Real dp2( 0.0 );
	Real sdp2( 0.0 );
	for ( Size i=ss2; i<=ss2+1; ++i ) {
		Vector temp;
		if ( i == ss2+1 ) {
			temp = dist_co_inv * ( bb_pos.C(i) - bb_pos.O(i) );
		} else {
			temp = dist_co_inv * ( bb_pos.O(i) - bb_pos.C(i) );
		}
		Real const tempdot = temp.dot( u_midvec );
		dp2 += std::abs(tempdot);
		sdp2 += tempdot;
	}
	dp2 *= 0.5;

	dp = dp1 + dp2;

	//js These signs tell whether the first c=o bond vector of a dimer points
	//js at the other dimer.  sign1 = 1 means that the first c=o bond of dimer1
	//js points at dimer2.  sign2 = 1 means that dimer2 points at dimer1.  When
	//js sign1 or sign2 equals 2, that dimer points away from the other dimer.
	sign1 = ( sdp1 > 0.0 ? 2 : 1 );
	sign2 = ( sdp2 < 0.0 ? 2 : 1 );
}


/// @brief
void
SSPairPotential::score(
	Pose const & pose,
	SS_Info2 const & ss_info,
	DimerPairings & dimer_pairs,
	Real & ss_score ) const
{
	using protocols::fldsgn::potentials::sspot::spherical;
	using protocols::fldsgn::potentials::sspot::get_foldtree_seqsep;
	using protocols::fldsgn::topology::StrandOP;
	using protocols::fldsgn::topology::Strands;

	ss_score = 0.0;
	Real rsigma_score = 0.0;
	dimer_pairs.clear();

	Real ssdist_12_ = -1.6408;

	Strands const & strands( ss_info.strands() );
	BB_Pos const & bb_pos( ss_info.bb_pos() );

	EnergyGraph const & energy_graph( pose.energies().energy_graph() );

	for ( Size istrand=1; istrand<=strands.size(); istrand++ ) {

		StrandOP const strand( strands[ istrand ] );
		if ( strand->length() < 2 ) continue;
		for ( Size ss1=strand->begin(); ss1<strand->end(); ss1++ ) {

			for ( core::graph::Graph::EdgeListConstIter
					iru  = energy_graph.get_node( ss1 )->const_upper_edge_list_begin(),
					irue = energy_graph.get_node( ss1 )->const_upper_edge_list_end();
					iru != irue; ++iru ) {

				Size const ss2( (*iru)->get_second_node_ind() );
				Size jstrand = ss_info.strand_id( ss2 );
				if ( pose.residue_type( ss2 ).is_upper_terminus() ) continue;
				if ( ss_info.strand_id( ss2+1 ) == 0 || jstrand == 0 || istrand == jstrand ) continue;

				Size const strand_seqsep = get_foldtree_seqsep( pose, strands[ jstrand ]->begin()-1, strands[ istrand ]->end()+1 ) + 1;
				Size const dimer_seqsep  = get_foldtree_seqsep( pose, ss2, ss1 );

				Vector const & pt1( bb_pos.N( ss1   ) );
				Vector const & pt2( bb_pos.C( ss1+1 ) );
				Vector const & pt3( bb_pos.N( ss2   ) );
				Vector const & pt4( bb_pos.C( ss2+1 ) );

				// midpoint coordinates of two dimers
				Vector cen1 = Real( 0.5 )*( pt1 + pt2 );
				Vector cen2 = Real( 0.5 )*( pt3 + pt4 );

				// vector between midpoints
				Vector mid_vector = cen2 -cen1;
				Real dist_dimers = mid_vector.length();

				if ( dist_dimers < strand_dist_cutoff_ ) {

					// calc phi and theta score between dimers
					Real phi, theta;
					spherical( pt2, pt4, phi, theta, cen1, cen2, mid_vector );
					Real phithetascore = calc_phithetascore( strand_seqsep, phi, theta );

					if ( phithetascore < 0.0 ) {

						// calc dot product and dist score
						// car if the sequence distance between pairs is too small, don't add this term to the total
						// car LOCAL STRANDS DON'T GET SO MUCH OF A SCORE BONUS, add terms to pair only if pair separation > dimer_seqsep_cutoff_
						Real dotscore( 0.0 );
						Real distscore( 0.0 );
						Real dimer_pair_score( 0.0 );
						Real dpall;
						Size sign1, sign2;
						pair_dp( ss1, ss2, bb_pos, dpall, mid_vector, sign1, sign2 );
						dotscore = calc_dotscore( dpall );

						// rsigma score
						Real const u21z = ( pt2 - cen1 ).normalized_or_zero().dot( mid_vector.normalized_or_zero() );
						Real sig = numeric::conversions::degrees( numeric::arccos( u21z ) ); //std::acos( sin_cos_range( u21z ) );
						Real tempscore_rsigma = calc_rsigmascore( sig, dist_dimers, sign1, sign2 );
						if ( tempscore_rsigma > 0.0 ) {
							//car add in all unfavorable rsigma scores (these pairs may not contribute
							//car to ss_score if they have a favorable dimer_pair_score)
							//car note there are no positive scores in structure_db.cc as of 7/03
							rsigma_score += tempscore_rsigma;
						} else {
							//car and favorable ones if phitheta and dot score favorable & dist<6.5
							//car note that rsigma is not subject to the favorable-dimer-interactions-must
							//car be-consistent-pairwise-strands rule that is applied to the ss_score below
							if ( dotscore < 0.0 ) rsigma_score += tempscore_rsigma;
						}


						if ( dimer_seqsep >= dimer_seqsep_cutoff_ ) {
							// dist score
							if ( lowstrand_ > 0.5 ) {  // what is lowstrand ??????
								distscore = ssdist_12_;  // I don't think this make asense
							}
							dimer_pair_score = phithetascore + distscore + dotscore;
						} else {
							//// rhiju Allow penalty for locally paired strands
							//Real penalty_localpair( 0.0 );
							//Real seqsep = std::abs( ss2 - ss1 );
							//if ( seqsep <= 11 ) penalty_localpair = localstrandpair_penalty_ * SS_penalty( seqsep );
							// sum all scores
							//Real dimer_pair_score = phithetascore + distscore + dotscore + penalty_localpair;
						}

						if ( dimer_pair_score < 0.0 ) {
							dimer_pairs.push_back( topology::DimerPairingOP( new DimerPairing( ss1, ss2, dist_dimers, phi, theta,
								sig, dpall, sign1, sign2, dimer_pair_score ) ) );
						} else {
							ss_score += dimer_pair_score;
						}

					} // if phithetascore < 0.0
				} // if ( dist_dimers < strand_dist_cutoff_ )
			} //
		} // for ( ss1 )
	} // for( istrand )

	dimer_pairs.finalize( ss_info );
	for ( DimerPairings::iterator it= dimer_pairs.begin(), ite= dimer_pairs.end(); it != ite; ++it ) {
		DimerPairing const & pairing( **it );
		if ( !pairing.valid() ) continue;
		ss_score += pairing.score();

		int const res1 ( pairing.res1() );
		int const res2 ( pairing.res2() );
		Size const strand1( ss_info.strand_id( res1 ) );
		Size const strand2( ss_info.strand_id( res2 ) );
		int const edge1( strands[ strand1 ]->end()+1 );
		int const edge2( strands[ strand2 ]->begin()-1 );

		if ( std::abs( res2 - res1 ) > 11 && std::abs( edge2-edge1 ) > 8 ) {
			ss_score -= 0.2;
		}
	}

	ss_score *= 0.498 * 0.75;
	ss_score += rsigma_score * 0.1;

} // score


/// @brief load phi/theta bins for use in secondary structure scoring
void
SSPairPotential::load_phi_theta_bins( String const & ss_filename )
{
	using ObjexxFCL::format::skip;
	typedef ObjexxFCL::FArray3D< Real > FArray3D_real;
	FArray3D_real pts_SS( 36, 36, 3 );

	FArray1D_int iptsn( 36 );
	for ( int itheta = 1; itheta <= 36; ++itheta ) {
		iptsn(itheta) = 100;
	}

	// FIXME: need equivalent to open_data_file() function here
	utility::io::izstream SS_stream;
	basic::database::open( SS_stream, ss_filename );
	for ( int isep = 1; isep <= 3; ++isep ) {
		for ( int itheta = 1; itheta <= 36; ++itheta ) {
			for ( int iph = 1; iph <= 36; ++iph ) {
				SS_stream >> pts_SS(itheta,iph,isep) >> skip;
			}
		}
		if ( isep == 1 ) SS_stream.seek_beg();
	}
	SS_stream.close();
	SS_stream.clear();

	for ( int isep = 1; isep <= 3; ++isep ) {
		Real tot = 0.0;
		Real totn = 0.0;
		for ( int iph = 1; iph <= 36; ++iph ) {
			for ( int itheta = 1; itheta <= 36; ++itheta ) {
				pts_SS( itheta, iph, isep ) += iptsn( itheta )*0.000001f;  //  SMALL COUNTS CORRECTION
				tot += pts_SS( itheta, iph, isep );
				totn += iptsn( itheta );
			}
		}
		for ( int iph = 1; iph <= 36; ++iph ) {
			for ( int itheta = 1; itheta <= 36; ++itheta ) {
				phithetascore_( 2, isep, iph, itheta ) = -std::log( pts_SS( itheta, iph, isep )/tot) + std::log(iptsn( itheta )/totn );
			}
		}
	}

} // load_phi_theta_bins

/// @brief
void
SSPairPotential::load_dotscore_bins()
{
	// triangle-2 random numbers
	//     data idsn/56,167,278,278,167,56/
	// sort of triangle-4 random numbers
	FArray1D_int idsn( 6 );
	idsn( 1 ) = 5596;
	idsn( 2 ) = 16581;
	idsn( 3 ) = 27823;
	idsn( 4 ) = 27823;
	idsn( 5 ) = 16581;
	idsn( 6 ) = 5596;

	FArray1D_int ids( 6 );
	ids( 1 ) = 1;
	ids( 2 ) = 48;
	ids( 3 ) = 368;
	ids( 4 ) = 2378;
	ids( 5 ) = 7141;
	ids( 6 ) = 8904;

	Real tot = 0.0;
	Real totn = 0.0;
	for ( int idot = 1; idot <= 6; ++idot ) {
		tot += ids(idot);
		totn += idsn(idot);
	}
	for ( int idot = 1; idot <= 6; ++idot ) {
		if ( ids(idot) != 0 ) {
			dotscore_(idot) = -std::log(ids(idot)/tot) + std::log(idsn(idot)/totn);
		} else {
			dotscore_(idot) = 0.0;
		}
	}
}

/// @brief
void
SSPairPotential::rsigma_dot_initializer(
	FArray4D_real & rsigma_dot
)
{
	core::scoring::SecondaryStructurePotential ssp;
	ssp.rsigma_dot_initializer( rsigma_dot );
}
/* // section 12
//js --------------------------
//js new rsigma stats that take into account whether the first
//js c=o bond vector points towards away from the other dimer
rsigma_dot(  1,  1,  1,  2 ) = -1.038100;
rsigma_dot(  1,  2,  1,  2 ) = -1.038100;
rsigma_dot(  1,  3,  1,  2 ) = -1.038100;
rsigma_dot(  1,  4,  1,  2 ) = -1.038100;
rsigma_dot(  1,  5,  1,  2 ) = -1.038100;
rsigma_dot(  1,  6,  1,  2 ) = -1.038100;
rsigma_dot(  1,  7,  1,  2 ) = -1.038100;
rsigma_dot(  1,  8,  1,  2 ) = -2.984000;
rsigma_dot(  1,  9,  1,  2 ) = -3.746110;
rsigma_dot(  1, 10,  1,  2 ) = -3.746110;
rsigma_dot(  1, 11,  1,  2 ) = -2.647500;
rsigma_dot(  1, 12,  1,  2 ) = -1.038100;
rsigma_dot(  1, 13,  1,  2 ) = -1.038100;
rsigma_dot(  1, 14,  1,  2 ) = -1.038100;
rsigma_dot(  1, 15,  1,  2 ) = -1.038100;
rsigma_dot(  1, 16,  1,  2 ) = -1.038100;
rsigma_dot(  1, 17,  1,  2 ) = -1.038100;
rsigma_dot(  1, 18,  1,  2 ) = -1.038100;
rsigma_dot(  2,  1,  1,  2 ) = -0.973500;
rsigma_dot(  2,  2,  1,  2 ) = -0.973500;
rsigma_dot(  2,  3,  1,  2 ) = -0.973500;
rsigma_dot(  2,  4,  1,  2 ) = -0.973500;
rsigma_dot(  2,  5,  1,  2 ) = -0.973500;
rsigma_dot(  2,  6,  1,  2 ) = -0.973500;
rsigma_dot(  2,  7,  1,  2 ) = -2.072100;
rsigma_dot(  2,  8,  1,  2 ) = -0.973500;
rsigma_dot(  2,  9,  1,  2 ) = -4.943820;
rsigma_dot(  2, 10,  1,  2 ) = -5.178220;
rsigma_dot(  2, 11,  1,  2 ) = -4.823670;
rsigma_dot(  2, 12,  1,  2 ) = -2.919400;
rsigma_dot(  2, 13,  1,  2 ) = -2.583000;
rsigma_dot(  2, 14,  1,  2 ) = -0.973500;
rsigma_dot(  2, 15,  1,  2 ) = -0.973500;
rsigma_dot(  2, 16,  1,  2 ) = -0.973500;
rsigma_dot(  2, 17,  1,  2 ) = -0.973500;
rsigma_dot(  2, 18,  1,  2 ) = -0.973500;
rsigma_dot(  3,  1,  1,  2 ) = -0.912900;
rsigma_dot(  3,  2,  1,  2 ) = -0.912900;
rsigma_dot(  3,  3,  1,  2 ) = -0.912900;
rsigma_dot(  3,  4,  1,  2 ) = -0.912900;
rsigma_dot(  3,  5,  1,  2 ) = -0.912900;
rsigma_dot(  3,  6,  1,  2 ) = -0.912900;
rsigma_dot(  3,  7,  1,  2 ) = -2.011500;
rsigma_dot(  3,  8,  1,  2 ) = -4.048400;
rsigma_dot(  3,  9,  1,  2 ) = -5.981810;
rsigma_dot(  3, 10,  1,  2 ) = -7.006470;
rsigma_dot(  3, 11,  1,  2 ) = -6.696730;
rsigma_dot(  3, 12,  1,  2 ) = -4.131780;
rsigma_dot(  3, 13,  1,  2 ) = -0.912900;
rsigma_dot(  3, 14,  1,  2 ) = -0.912900;
rsigma_dot(  3, 15,  1,  2 ) = -0.912900;
rsigma_dot(  3, 16,  1,  2 ) = -0.912900;
rsigma_dot(  3, 17,  1,  2 ) = -0.912900;
rsigma_dot(  3, 18,  1,  2 ) = -0.912900;
rsigma_dot(  4,  1,  1,  2 ) = -0.855700;
rsigma_dot(  4,  2,  1,  2 ) = -0.855700;
rsigma_dot(  4,  3,  1,  2 ) = -0.855700;
rsigma_dot(  4,  4,  1,  2 ) = -0.855700;
rsigma_dot(  4,  5,  1,  2 ) = -0.855700;
rsigma_dot(  4,  6,  1,  2 ) = -2.465200;
rsigma_dot(  4,  7,  1,  2 ) = -3.420700;
rsigma_dot(  4,  8,  1,  2 ) = -3.900270;
rsigma_dot(  4,  9,  1,  2 ) = -6.315330;
rsigma_dot(  4, 10,  1,  2 ) = -8.116270;
rsigma_dot(  4, 11,  1,  2 ) = -8.329380;
rsigma_dot(  4, 12,  1,  2 ) = -5.298390;
rsigma_dot(  4, 13,  1,  2 ) = -1.954400;
rsigma_dot(  4, 14,  1,  2 ) = -0.855700;
rsigma_dot(  4, 15,  1,  2 ) = -0.855700;
rsigma_dot(  4, 16,  1,  2 ) = -0.855700;
rsigma_dot(  4, 17,  1,  2 ) = -0.855700;
rsigma_dot(  4, 18,  1,  2 ) = -0.855700;
rsigma_dot(  5,  1,  1,  2 ) = -0.801700;
rsigma_dot(  5,  2,  1,  2 ) = -0.801700;
rsigma_dot(  5,  3,  1,  2 ) = -0.801700;
rsigma_dot(  5,  4,  1,  2 ) = -0.801700;
rsigma_dot(  5,  5,  1,  2 ) = -0.801700;
rsigma_dot(  5,  6,  1,  2 ) = -1.900300;
rsigma_dot(  5,  7,  1,  2 ) = -3.634900;
rsigma_dot(  5,  8,  1,  2 ) = -5.119160;
rsigma_dot(  5,  9,  1,  2 ) = -6.381410;
rsigma_dot(  5, 10,  1,  2 ) = -9.110860;
rsigma_dot(  5, 11,  1,  2 ) = -9.615860;
rsigma_dot(  5, 12,  1,  2 ) = -6.566870;
rsigma_dot(  5, 13,  1,  2 ) = -2.411100;
rsigma_dot(  5, 14,  1,  2 ) = -1.900300;
rsigma_dot(  5, 15,  1,  2 ) = -0.801700;
rsigma_dot(  5, 16,  1,  2 ) = -0.801700;
rsigma_dot(  5, 17,  1,  2 ) = -0.801700;
rsigma_dot(  5, 18,  1,  2 ) = -0.801700;
rsigma_dot(  6,  1,  1,  2 ) = -0.750400;
rsigma_dot(  6,  2,  1,  2 ) = -0.750400;
rsigma_dot(  6,  3,  1,  2 ) = -0.750400;
rsigma_dot(  6,  4,  1,  2 ) = -0.750400;
rsigma_dot(  6,  5,  1,  2 ) = -0.750400;
rsigma_dot(  6,  6,  1,  2 ) = -2.359800;
rsigma_dot(  6,  7,  1,  2 ) = -4.511580;
rsigma_dot(  6,  8,  1,  2 ) = -5.325090;
rsigma_dot(  6,  9,  1,  2 ) = -6.768980;
rsigma_dot(  6, 10,  1,  2 ) = -8.613260;
rsigma_dot(  6, 11,  1,  2 ) = -9.405770;
rsigma_dot(  6, 12,  1,  2 ) = -7.032650;
rsigma_dot(  6, 13,  1,  2 ) = -3.885880;
rsigma_dot(  6, 14,  1,  2 ) = -3.315300;
rsigma_dot(  6, 15,  1,  2 ) = -0.750400;
rsigma_dot(  6, 16,  1,  2 ) = -0.750400;
rsigma_dot(  6, 17,  1,  2 ) = -0.750400;
rsigma_dot(  6, 18,  1,  2 ) = -0.750400;
rsigma_dot(  7,  1,  1,  2 ) = -0.701600;
rsigma_dot(  7,  2,  1,  2 ) = -0.701600;
rsigma_dot(  7,  3,  1,  2 ) = -0.701600;
rsigma_dot(  7,  4,  1,  2 ) = -0.701600;
rsigma_dot(  7,  5,  1,  2 ) = -0.701600;
rsigma_dot(  7,  6,  1,  2 ) = -2.311000;
rsigma_dot(  7,  7,  1,  2 ) = -4.365150;
rsigma_dot(  7,  8,  1,  2 ) = -5.446520;
rsigma_dot(  7,  9,  1,  2 ) = -6.795160;
rsigma_dot(  7, 10,  1,  2 ) = -7.145720;
rsigma_dot(  7, 11,  1,  2 ) = -7.527050;
rsigma_dot(  7, 12,  1,  2 ) = -6.885740;
rsigma_dot(  7, 13,  1,  2 ) = -4.508250;
rsigma_dot(  7, 14,  1,  2 ) = -1.800200;
rsigma_dot(  7, 15,  1,  2 ) = -2.647500;
rsigma_dot(  7, 16,  1,  2 ) = -0.701600;
rsigma_dot(  7, 17,  1,  2 ) = -0.701600;
rsigma_dot(  7, 18,  1,  2 ) = -0.701600;
rsigma_dot(  8,  1,  1,  2 ) = -0.655100;
rsigma_dot(  8,  2,  1,  2 ) = -0.655100;
rsigma_dot(  8,  3,  1,  2 ) = -0.655100;
rsigma_dot(  8,  4,  1,  2 ) = -0.655100;
rsigma_dot(  8,  5,  1,  2 ) = -0.655100;
rsigma_dot(  8,  6,  1,  2 ) = -3.699600;
rsigma_dot(  8,  7,  1,  2 ) = -4.089060;
rsigma_dot(  8,  8,  1,  2 ) = -5.417250;
rsigma_dot(  8,  9,  1,  2 ) = -6.293430;
rsigma_dot(  8, 10,  1,  2 ) = -6.663890;
rsigma_dot(  8, 11,  1,  2 ) = -6.788470;
rsigma_dot(  8, 12,  1,  2 ) = -6.368800;
rsigma_dot(  8, 13,  1,  2 ) = -4.732610;
rsigma_dot(  8, 14,  1,  2 ) = -1.753700;
rsigma_dot(  8, 15,  1,  2 ) = -0.655100;
rsigma_dot(  8, 16,  1,  2 ) = -0.655100;
rsigma_dot(  8, 17,  1,  2 ) = -0.655100;
rsigma_dot(  8, 18,  1,  2 ) = -0.655100;
rsigma_dot(  9,  1,  1,  2 ) = -0.610600;
rsigma_dot(  9,  2,  1,  2 ) = -0.610600;
rsigma_dot(  9,  3,  1,  2 ) = -0.610600;
rsigma_dot(  9,  4,  1,  2 ) = -0.610600;
rsigma_dot(  9,  5,  1,  2 ) = -1.709200;
rsigma_dot(  9,  6,  1,  2 ) = -3.906460;
rsigma_dot(  9,  7,  1,  2 ) = -3.906460;
rsigma_dot(  9,  8,  1,  2 ) = -5.355550;
rsigma_dot(  9,  9,  1,  2 ) = -6.070210;
rsigma_dot(  9, 10,  1,  2 ) = -5.971910;
rsigma_dot(  9, 11,  1,  2 ) = -6.070210;
rsigma_dot(  9, 12,  1,  2 ) = -5.852370;
rsigma_dot(  9, 13,  1,  2 ) = -4.901080;
rsigma_dot(  9, 14,  1,  2 ) = -2.807800;
rsigma_dot(  9, 15,  1,  2 ) = -2.220100;
rsigma_dot(  9, 16,  1,  2 ) = -0.610600;
rsigma_dot(  9, 17,  1,  2 ) = -1.709200;
rsigma_dot(  9, 18,  1,  2 ) = -0.610600;
rsigma_dot( 10,  1,  1,  2 ) = -0.568100;
rsigma_dot( 10,  2,  1,  2 ) = -0.568100;
rsigma_dot( 10,  3,  1,  2 ) = -0.568100;
rsigma_dot( 10,  4,  1,  2 ) = -2.177500;
rsigma_dot( 10,  5,  1,  2 ) = -2.765300;
rsigma_dot( 10,  6,  1,  2 ) = -4.002050;
rsigma_dot( 10,  7,  1,  2 ) = -4.802170;
rsigma_dot( 10,  8,  1,  2 ) = -5.598500;
rsigma_dot( 10,  9,  1,  2 ) = -5.809810;
rsigma_dot( 10, 10,  1,  2 ) = -5.809810;
rsigma_dot( 10, 11,  1,  2 ) = -5.732850;
rsigma_dot( 10, 12,  1,  2 ) = -5.572010;
rsigma_dot( 10, 13,  1,  2 ) = -4.858520;
rsigma_dot( 10, 14,  1,  2 ) = -3.935360;
rsigma_dot( 10, 15,  1,  2 ) = -1.666700;
rsigma_dot( 10, 16,  1,  2 ) = -1.666700;
rsigma_dot( 10, 17,  1,  2 ) = -0.568100;
rsigma_dot( 10, 18,  1,  2 ) = -1.666700;
rsigma_dot( 11,  1,  1,  2 ) = -0.527200;
rsigma_dot( 11,  2,  1,  2 ) = -0.527200;
rsigma_dot( 11,  3,  1,  2 ) = -0.527200;
rsigma_dot( 11,  4,  1,  2 ) = -1.625900;
rsigma_dot( 11,  5,  1,  2 ) = -3.746110;
rsigma_dot( 11,  6,  1,  2 ) = -4.190800;
rsigma_dot( 11,  7,  1,  2 ) = -4.761350;
rsigma_dot( 11,  8,  1,  2 ) = -5.897880;
rsigma_dot( 11,  9,  1,  2 ) = -5.850250;
rsigma_dot( 11, 10,  1,  2 ) = -5.789930;
rsigma_dot( 11, 11,  1,  2 ) = -5.371430;
rsigma_dot( 11, 12,  1,  2 ) = -5.608640;
rsigma_dot( 11, 13,  1,  2 ) = -5.306360;
rsigma_dot( 11, 14,  1,  2 ) = -4.288440;
rsigma_dot( 11, 15,  1,  2 ) = -3.092200;
rsigma_dot( 11, 16,  1,  2 ) = -0.527200;
rsigma_dot( 11, 17,  1,  2 ) = -0.527200;
rsigma_dot( 11, 18,  1,  2 ) = -0.527200;
rsigma_dot( 12,  1,  1,  2 ) = -0.488000;
rsigma_dot( 12,  2,  1,  2 ) = -0.488000;
rsigma_dot( 12,  3,  1,  2 ) = -0.488000;
rsigma_dot( 12,  4,  1,  2 ) = -3.623500;
rsigma_dot( 12,  5,  1,  2 ) = -4.722120;
rsigma_dot( 12,  6,  1,  2 ) = -4.692710;
rsigma_dot( 12,  7,  1,  2 ) = -5.197550;
rsigma_dot( 12,  8,  1,  2 ) = -5.544260;
rsigma_dot( 12,  9,  1,  2 ) = -5.740290;
rsigma_dot( 12, 10,  1,  2 ) = -5.478450;
rsigma_dot( 12, 11,  1,  2 ) = -5.347830;
rsigma_dot( 12, 12,  1,  2 ) = -5.363220;
rsigma_dot( 12, 13,  1,  2 ) = -5.103140;
rsigma_dot( 12, 14,  1,  2 ) = -4.722120;
rsigma_dot( 12, 15,  1,  2 ) = -3.984530;
rsigma_dot( 12, 16,  1,  2 ) = -3.196100;
rsigma_dot( 12, 17,  1,  2 ) = -0.488000;
rsigma_dot( 12, 18,  1,  2 ) = -0.488000;

// section 11
rsigma_dot(  1,  1,  1,  1 ) = -0.552300;
rsigma_dot(  1,  2,  1,  1 ) = -0.552300;
rsigma_dot(  1,  3,  1,  1 ) = -0.552300;
rsigma_dot(  1,  4,  1,  1 ) = -0.552300;
rsigma_dot(  1,  5,  1,  1 ) = -0.552300;
rsigma_dot(  1,  6,  1,  1 ) = -0.552300;
rsigma_dot(  1,  7,  1,  1 ) = -0.552300;
rsigma_dot(  1,  8,  1,  1 ) = -0.552300;
rsigma_dot(  1,  9,  1,  1 ) = -0.552300;
rsigma_dot(  1, 10,  1,  1 ) = -1.651000;
rsigma_dot(  1, 11,  1,  1 ) = -0.552300;
rsigma_dot(  1, 12,  1,  1 ) = -0.552300;
rsigma_dot(  1, 13,  1,  1 ) = -0.552300;
rsigma_dot(  1, 14,  1,  1 ) = -0.552300;
rsigma_dot(  1, 15,  1,  1 ) = -0.552300;
rsigma_dot(  1, 16,  1,  1 ) = -0.552300;
rsigma_dot(  1, 17,  1,  1 ) = -0.552300;
rsigma_dot(  1, 18,  1,  1 ) = -0.552300;
rsigma_dot(  2,  1,  1,  1 ) = -0.487800;
rsigma_dot(  2,  2,  1,  1 ) = -0.487800;
rsigma_dot(  2,  3,  1,  1 ) = -0.487800;
rsigma_dot(  2,  4,  1,  1 ) = -0.487800;
rsigma_dot(  2,  5,  1,  1 ) = -0.487800;
rsigma_dot(  2,  6,  1,  1 ) = -0.487800;
rsigma_dot(  2,  7,  1,  1 ) = -0.487800;
rsigma_dot(  2,  8,  1,  1 ) = -0.487800;
rsigma_dot(  2,  9,  1,  1 ) = -0.487800;
rsigma_dot(  2, 10,  1,  1 ) = -1.586400;
rsigma_dot(  2, 11,  1,  1 ) = -2.097200;
rsigma_dot(  2, 12,  1,  1 ) = -1.586400;
rsigma_dot(  2, 13,  1,  1 ) = -0.487800;
rsigma_dot(  2, 14,  1,  1 ) = -0.487800;
rsigma_dot(  2, 15,  1,  1 ) = -0.487800;
rsigma_dot(  2, 16,  1,  1 ) = -0.487800;
rsigma_dot(  2, 17,  1,  1 ) = -0.487800;
rsigma_dot(  2, 18,  1,  1 ) = -0.487800;
rsigma_dot(  3,  1,  1,  1 ) = -0.427200;
rsigma_dot(  3,  2,  1,  1 ) = -0.427200;
rsigma_dot(  3,  3,  1,  1 ) = -0.427200;
rsigma_dot(  3,  4,  1,  1 ) = -0.427200;
rsigma_dot(  3,  5,  1,  1 ) = -0.427200;
rsigma_dot(  3,  6,  1,  1 ) = -1.525800;
rsigma_dot(  3,  7,  1,  1 ) = -2.036600;
rsigma_dot(  3,  8,  1,  1 ) = -2.036600;
rsigma_dot(  3,  9,  1,  1 ) = -0.427200;
rsigma_dot(  3, 10,  1,  1 ) = -2.373100;
rsigma_dot(  3, 11,  1,  1 ) = -1.525800;
rsigma_dot(  3, 12,  1,  1 ) = -2.992100;
rsigma_dot(  3, 13,  1,  1 ) = -2.373100;
rsigma_dot(  3, 14,  1,  1 ) = -0.427200;
rsigma_dot(  3, 15,  1,  1 ) = -0.427200;
rsigma_dot(  3, 16,  1,  1 ) = -0.427200;
rsigma_dot(  3, 17,  1,  1 ) = -0.427200;
rsigma_dot(  3, 18,  1,  1 ) = -0.427200;
rsigma_dot(  4,  1,  1,  1 ) = -0.370000;
rsigma_dot(  4,  2,  1,  1 ) = -0.370000;
rsigma_dot(  4,  3,  1,  1 ) = -0.370000;
rsigma_dot(  4,  4,  1,  1 ) = -0.370000;
rsigma_dot(  4,  5,  1,  1 ) = -0.370000;
rsigma_dot(  4,  6,  1,  1 ) = -0.370000;
rsigma_dot(  4,  7,  1,  1 ) = -2.315900;
rsigma_dot(  4,  8,  1,  1 ) = -2.567200;
rsigma_dot(  4,  9,  1,  1 ) = -2.767900;
rsigma_dot(  4, 10,  1,  1 ) = -3.588900;
rsigma_dot(  4, 11,  1,  1 ) = -3.414500;
rsigma_dot(  4, 12,  1,  1 ) = -4.033580;
rsigma_dot(  4, 13,  1,  1 ) = -3.804010;
rsigma_dot(  4, 14,  1,  1 ) = -0.370000;
rsigma_dot(  4, 15,  1,  1 ) = -0.370000;
rsigma_dot(  4, 16,  1,  1 ) = -0.370000;
rsigma_dot(  4, 17,  1,  1 ) = -0.370000;
rsigma_dot(  4, 18,  1,  1 ) = -0.370000;
rsigma_dot(  5,  1,  1,  1 ) = -0.316000;
rsigma_dot(  5,  2,  1,  1 ) = -0.316000;
rsigma_dot(  5,  3,  1,  1 ) = -0.316000;
rsigma_dot(  5,  4,  1,  1 ) = -0.316000;
rsigma_dot(  5,  5,  1,  1 ) = -0.316000;
rsigma_dot(  5,  6,  1,  1 ) = -1.925400;
rsigma_dot(  5,  7,  1,  1 ) = -3.683200;
rsigma_dot(  5,  8,  1,  1 ) = -4.122610;
rsigma_dot(  5,  9,  1,  1 ) = -3.812460;
rsigma_dot(  5, 10,  1,  1 ) = -3.683200;
rsigma_dot(  5, 11,  1,  1 ) = -3.360500;
rsigma_dot(  5, 12,  1,  1 ) = -4.969910;
rsigma_dot(  5, 13,  1,  1 ) = -4.869830;
rsigma_dot(  5, 14,  1,  1 ) = -3.149200;
rsigma_dot(  5, 15,  1,  1 ) = -0.316000;
rsigma_dot(  5, 16,  1,  1 ) = -0.316000;
rsigma_dot(  5, 17,  1,  1 ) = -0.316000;
rsigma_dot(  5, 18,  1,  1 ) = -0.316000;
rsigma_dot(  6,  1,  1,  1 ) = -0.264700;
rsigma_dot(  6,  2,  1,  1 ) = -0.264700;
rsigma_dot(  6,  3,  1,  1 ) = -0.264700;
rsigma_dot(  6,  4,  1,  1 ) = -0.264700;
rsigma_dot(  6,  5,  1,  1 ) = -0.264700;
rsigma_dot(  6,  6,  1,  1 ) = -2.829600;
rsigma_dot(  6,  7,  1,  1 ) = -5.680760;
rsigma_dot(  6,  8,  1,  1 ) = -5.076840;
rsigma_dot(  6,  9,  1,  1 ) = -4.271990;
rsigma_dot(  6, 10,  1,  1 ) = -4.234950;
rsigma_dot(  6, 11,  1,  1 ) = -3.761170;
rsigma_dot(  6, 12,  1,  1 ) = -5.567960;
rsigma_dot(  6, 13,  1,  1 ) = -6.108200;
rsigma_dot(  6, 14,  1,  1 ) = -4.956010;
rsigma_dot(  6, 15,  1,  1 ) = -1.874100;
rsigma_dot(  6, 16,  1,  1 ) = -0.264700;
rsigma_dot(  6, 17,  1,  1 ) = -0.264700;
rsigma_dot(  6, 18,  1,  1 ) = -0.264700;
rsigma_dot(  7,  1,  1,  1 ) = -0.215900;
rsigma_dot(  7,  2,  1,  1 ) = -0.215900;
rsigma_dot(  7,  3,  1,  1 ) = -0.215900;
rsigma_dot(  7,  4,  1,  1 ) = -0.215900;
rsigma_dot(  7,  5,  1,  1 ) = -0.215900;
rsigma_dot(  7,  6,  1,  1 ) = -5.011660;
rsigma_dot(  7,  7,  1,  1 ) = -7.313420;
rsigma_dot(  7,  8,  1,  1 ) = -6.580620;
rsigma_dot(  7,  9,  1,  1 ) = -4.726730;
rsigma_dot(  7, 10,  1,  1 ) = -4.748470;
rsigma_dot(  7, 11,  1,  1 ) = -2.923900;
rsigma_dot(  7, 12,  1,  1 ) = -5.700670;
rsigma_dot(  7, 13,  1,  1 ) = -6.799280;
rsigma_dot(  7, 14,  1,  1 ) = -6.059410;
rsigma_dot(  7, 15,  1,  1 ) = -2.161800;
rsigma_dot(  7, 16,  1,  1 ) = -0.215900;
rsigma_dot(  7, 17,  1,  1 ) = -0.215900;
rsigma_dot(  7, 18,  1,  1 ) = -0.215900;
rsigma_dot(  8,  1,  1,  1 ) = -0.169300;
rsigma_dot(  8,  2,  1,  1 ) = -0.169300;
rsigma_dot(  8,  3,  1,  1 ) = -0.169300;
rsigma_dot(  8,  4,  1,  1 ) = -0.169300;
rsigma_dot(  8,  5,  1,  1 ) = -0.169300;
rsigma_dot(  8,  6,  1,  1 ) = -6.132930;
rsigma_dot(  8,  7,  1,  1 ) = -8.735140;
rsigma_dot(  8,  8,  1,  1 ) = -7.789560;
rsigma_dot(  8,  9,  1,  1 ) = -5.199790;
rsigma_dot(  8, 10,  1,  1 ) = -4.981530;
rsigma_dot(  8, 11,  1,  1 ) = -3.882920;
rsigma_dot(  8, 12,  1,  1 ) = -5.059700;
rsigma_dot(  8, 13,  1,  1 ) = -7.329420;
rsigma_dot(  8, 14,  1,  1 ) = -7.217740;
rsigma_dot(  8, 15,  1,  1 ) = -4.723230;
rsigma_dot(  8, 16,  1,  1 ) = -0.169300;
rsigma_dot(  8, 17,  1,  1 ) = -0.169300;
rsigma_dot(  8, 18,  1,  1 ) = -0.169300;
rsigma_dot(  9,  1,  1,  1 ) = -0.124900;
rsigma_dot(  9,  2,  1,  1 ) = -0.124900;
rsigma_dot(  9,  3,  1,  1 ) = -0.124900;
rsigma_dot(  9,  4,  1,  1 ) = -1.223500;
rsigma_dot(  9,  5,  1,  1 ) = -3.069300;
rsigma_dot(  9,  6,  1,  1 ) = -6.756900;
rsigma_dot(  9,  7,  1,  1 ) = -8.997660;
rsigma_dot(  9,  8,  1,  1 ) = -7.292710;
rsigma_dot(  9,  9,  1,  1 ) = -5.704630;
rsigma_dot(  9, 10,  1,  1 ) = -4.442390;
rsigma_dot(  9, 11,  1,  1 ) = -4.202430;
rsigma_dot(  9, 12,  1,  1 ) = -4.567550;
rsigma_dot(  9, 13,  1,  1 ) = -7.550850;
rsigma_dot(  9, 14,  1,  1 ) = -8.291970;
rsigma_dot(  9, 15,  1,  1 ) = -5.617960;
rsigma_dot(  9, 16,  1,  1 ) = -1.223500;
rsigma_dot(  9, 17,  1,  1 ) = -1.223500;
rsigma_dot(  9, 18,  1,  1 ) = -0.124900;
rsigma_dot( 10,  1,  1,  1 ) = -0.082300;
rsigma_dot( 10,  2,  1,  1 ) = -0.082300;
rsigma_dot( 10,  3,  1,  1 ) = -0.082300;
rsigma_dot( 10,  4,  1,  1 ) = -1.180900;
rsigma_dot( 10,  5,  1,  1 ) = -3.932490;
rsigma_dot( 10,  6,  1,  1 ) = -6.858840;
rsigma_dot( 10,  7,  1,  1 ) = -8.186740;
rsigma_dot( 10,  8,  1,  1 ) = -6.228670;
rsigma_dot( 10,  9,  1,  1 ) = -5.720690;
rsigma_dot( 10, 10,  1,  1 ) = -4.548250;
rsigma_dot( 10, 11,  1,  1 ) = -4.451790;
rsigma_dot( 10, 12,  1,  1 ) = -4.548250;
rsigma_dot( 10, 13,  1,  1 ) = -7.242410;
rsigma_dot( 10, 14,  1,  1 ) = -8.619920;
rsigma_dot( 10, 15,  1,  1 ) = -6.349540;
rsigma_dot( 10, 16,  1,  1 ) = -3.026800;
rsigma_dot( 10, 17,  1,  1 ) = -1.180900;
rsigma_dot( 10, 18,  1,  1 ) = -1.180900;
rsigma_dot( 11,  1,  1,  1 ) = -0.041500;
rsigma_dot( 11,  2,  1,  1 ) = -0.041500;
rsigma_dot( 11,  3,  1,  1 ) = -0.041500;
rsigma_dot( 11,  4,  1,  1 ) = -1.987400;
rsigma_dot( 11,  5,  1,  1 ) = -4.275620;
rsigma_dot( 11,  6,  1,  1 ) = -6.582550;
rsigma_dot( 11,  7,  1,  1 ) = -6.915710;
rsigma_dot( 11,  8,  1,  1 ) = -6.166200;
rsigma_dot( 11,  9,  1,  1 ) = -5.483930;
rsigma_dot( 11, 10,  1,  1 ) = -4.595390;
rsigma_dot( 11, 11,  1,  1 ) = -4.616230;
rsigma_dot( 11, 12,  1,  1 ) = -4.931860;
rsigma_dot( 11, 13,  1,  1 ) = -6.246070;
rsigma_dot( 11, 14,  1,  1 ) = -8.296560;
rsigma_dot( 11, 15,  1,  1 ) = -6.734840;
rsigma_dot( 11, 16,  1,  1 ) = -3.337400;
rsigma_dot( 11, 17,  1,  1 ) = -1.140100;
rsigma_dot( 11, 18,  1,  1 ) = -0.041500;
rsigma_dot( 12,  1,  1,  1 ) = -0.002300;
rsigma_dot( 12,  2,  1,  1 ) = -0.002300;
rsigma_dot( 12,  3,  1,  1 ) = -0.002300;
rsigma_dot( 12,  4,  1,  1 ) = -2.199500;
rsigma_dot( 12,  5,  1,  1 ) = -4.729680;
rsigma_dot( 12,  6,  1,  1 ) = -6.131350;
rsigma_dot( 12,  7,  1,  1 ) = -5.729140;
rsigma_dot( 12,  8,  1,  1 ) = -5.604410;
rsigma_dot( 12,  9,  1,  1 ) = -5.211780;
rsigma_dot( 12, 10,  1,  1 ) = -4.892640;
rsigma_dot( 12, 11,  1,  1 ) = -4.992730;
rsigma_dot( 12, 12,  1,  1 ) = -5.275290;
rsigma_dot( 12, 13,  1,  1 ) = -5.945090;
rsigma_dot( 12, 14,  1,  1 ) = -7.716970;
rsigma_dot( 12, 15,  1,  1 ) = -6.693140;
rsigma_dot( 12, 16,  1,  1 ) = -3.972590;
rsigma_dot( 12, 17,  1,  1 ) = -1.948200;
rsigma_dot( 12, 18,  1,  1 ) = -0.002300;

// section 22
rsigma_dot(  1,  1,  2,  2 ) = -0.550000;
rsigma_dot(  1,  2,  2,  2 ) = -0.550000;
rsigma_dot(  1,  3,  2,  2 ) = -0.550000;
rsigma_dot(  1,  4,  2,  2 ) = -0.550000;
rsigma_dot(  1,  5,  2,  2 ) = -0.550000;
rsigma_dot(  1,  6,  2,  2 ) = -0.550000;
rsigma_dot(  1,  7,  2,  2 ) = -0.550000;
rsigma_dot(  1,  8,  2,  2 ) = -2.159400;
rsigma_dot(  1,  9,  2,  2 ) = -0.550000;
rsigma_dot(  1, 10,  2,  2 ) = -2.159400;
rsigma_dot(  1, 11,  2,  2 ) = -0.550000;
rsigma_dot(  1, 12,  2,  2 ) = -0.550000;
rsigma_dot(  1, 13,  2,  2 ) = -0.550000;
rsigma_dot(  1, 14,  2,  2 ) = -0.550000;
rsigma_dot(  1, 15,  2,  2 ) = -0.550000;
rsigma_dot(  1, 16,  2,  2 ) = -0.550000;
rsigma_dot(  1, 17,  2,  2 ) = -0.550000;
rsigma_dot(  1, 18,  2,  2 ) = -0.550000;
rsigma_dot(  2,  1,  2,  2 ) = -0.485500;
rsigma_dot(  2,  2,  2,  2 ) = -0.485500;
rsigma_dot(  2,  3,  2,  2 ) = -0.485500;
rsigma_dot(  2,  4,  2,  2 ) = -0.485500;
rsigma_dot(  2,  5,  2,  2 ) = -0.485500;
rsigma_dot(  2,  6,  2,  2 ) = -1.584100;
rsigma_dot(  2,  7,  2,  2 ) = -3.050400;
rsigma_dot(  2,  8,  2,  2 ) = -2.094900;
rsigma_dot(  2,  9,  2,  2 ) = -2.431400;
rsigma_dot(  2, 10,  2,  2 ) = -1.584100;
rsigma_dot(  2, 11,  2,  2 ) = -0.485500;
rsigma_dot(  2, 12,  2,  2 ) = -0.485500;
rsigma_dot(  2, 13,  2,  2 ) = -0.485500;
rsigma_dot(  2, 14,  2,  2 ) = -0.485500;
rsigma_dot(  2, 15,  2,  2 ) = -0.485500;
rsigma_dot(  2, 16,  2,  2 ) = -0.485500;
rsigma_dot(  2, 17,  2,  2 ) = -0.485500;
rsigma_dot(  2, 18,  2,  2 ) = -0.485500;
rsigma_dot(  3,  1,  2,  2 ) = -0.424800;
rsigma_dot(  3,  2,  2,  2 ) = -0.424800;
rsigma_dot(  3,  3,  2,  2 ) = -0.424800;
rsigma_dot(  3,  4,  2,  2 ) = -0.424800;
rsigma_dot(  3,  5,  2,  2 ) = -0.424800;
rsigma_dot(  3,  6,  2,  2 ) = -3.258000;
rsigma_dot(  3,  7,  2,  2 ) = -4.231500;
rsigma_dot(  3,  8,  2,  2 ) = -3.980180;
rsigma_dot(  3,  9,  2,  2 ) = -2.370700;
rsigma_dot(  3, 10,  2,  2 ) = -1.523400;
rsigma_dot(  3, 11,  2,  2 ) = -2.822700;
rsigma_dot(  3, 12,  2,  2 ) = -0.424800;
rsigma_dot(  3, 13,  2,  2 ) = -0.424800;
rsigma_dot(  3, 14,  2,  2 ) = -0.424800;
rsigma_dot(  3, 15,  2,  2 ) = -0.424800;
rsigma_dot(  3, 16,  2,  2 ) = -0.424800;
rsigma_dot(  3, 17,  2,  2 ) = -0.424800;
rsigma_dot(  3, 18,  2,  2 ) = -0.424800;
rsigma_dot(  4,  1,  2,  2 ) = -0.367700;
rsigma_dot(  4,  2,  2,  2 ) = -0.367700;
rsigma_dot(  4,  3,  2,  2 ) = -0.367700;
rsigma_dot(  4,  4,  2,  2 ) = -0.367700;
rsigma_dot(  4,  5,  2,  2 ) = -3.312100;
rsigma_dot(  4,  6,  2,  2 ) = -3.412200;
rsigma_dot(  4,  7,  2,  2 ) = -5.509340;
rsigma_dot(  4,  8,  2,  2 ) = -4.128880;
rsigma_dot(  4,  9,  2,  2 ) = -3.075700;
rsigma_dot(  4, 10,  2,  2 ) = -2.765600;
rsigma_dot(  4, 11,  2,  2 ) = -2.765600;
rsigma_dot(  4, 12,  2,  2 ) = -2.765600;
rsigma_dot(  4, 13,  2,  2 ) = -0.367700;
rsigma_dot(  4, 14,  2,  2 ) = -0.367700;
rsigma_dot(  4, 15,  2,  2 ) = -0.367700;
rsigma_dot(  4, 16,  2,  2 ) = -0.367700;
rsigma_dot(  4, 17,  2,  2 ) = -0.367700;
rsigma_dot(  4, 18,  2,  2 ) = -0.367700;
rsigma_dot(  5,  1,  2,  2 ) = -0.313600;
rsigma_dot(  5,  2,  2,  2 ) = -0.313600;
rsigma_dot(  5,  3,  2,  2 ) = -0.313600;
rsigma_dot(  5,  4,  2,  2 ) = -1.412200;
rsigma_dot(  5,  5,  2,  2 ) = -3.021700;
rsigma_dot(  5,  6,  2,  2 ) = -4.779520;
rsigma_dot(  5,  7,  2,  2 ) = -6.825360;
rsigma_dot(  5,  8,  2,  2 ) = -4.657420;
rsigma_dot(  5,  9,  2,  2 ) = -3.358100;
rsigma_dot(  5, 10,  2,  2 ) = -3.609400;
rsigma_dot(  5, 11,  2,  2 ) = -4.245440;
rsigma_dot(  5, 12,  2,  2 ) = -3.021700;
rsigma_dot(  5, 13,  2,  2 ) = -1.412200;
rsigma_dot(  5, 14,  2,  2 ) = -0.313600;
rsigma_dot(  5, 15,  2,  2 ) = -0.313600;
rsigma_dot(  5, 16,  2,  2 ) = -0.313600;
rsigma_dot(  5, 17,  2,  2 ) = -0.313600;
rsigma_dot(  5, 18,  2,  2 ) = -0.313600;
rsigma_dot(  6,  1,  2,  2 ) = -0.262300;
rsigma_dot(  6,  2,  2,  2 ) = -0.262300;
rsigma_dot(  6,  3,  2,  2 ) = -0.262300;
rsigma_dot(  6,  4,  2,  2 ) = -1.871800;
rsigma_dot(  6,  5,  2,  2 ) = -2.208200;
rsigma_dot(  6,  6,  2,  2 ) = -6.862190;
rsigma_dot(  6,  7,  2,  2 ) = -7.463490;
rsigma_dot(  6,  8,  2,  2 ) = -4.837030;
rsigma_dot(  6,  9,  2,  2 ) = -4.436700;
rsigma_dot(  6, 10,  2,  2 ) = -4.916280;
rsigma_dot(  6, 11,  2,  2 ) = -5.687270;
rsigma_dot(  6, 12,  2,  2 ) = -4.405450;
rsigma_dot(  6, 13,  2,  2 ) = -2.208200;
rsigma_dot(  6, 14,  2,  2 ) = -0.262300;
rsigma_dot(  6, 15,  2,  2 ) = -1.360900;
rsigma_dot(  6, 16,  2,  2 ) = -0.262300;
rsigma_dot(  6, 17,  2,  2 ) = -0.262300;
rsigma_dot(  6, 18,  2,  2 ) = -0.262300;
rsigma_dot(  7,  1,  2,  2 ) = -0.213500;
rsigma_dot(  7,  2,  2,  2 ) = -0.213500;
rsigma_dot(  7,  3,  2,  2 ) = -0.213500;
rsigma_dot(  7,  4,  2,  2 ) = -0.213500;
rsigma_dot(  7,  5,  2,  2 ) = -4.256580;
rsigma_dot(  7,  6,  2,  2 ) = -8.696750;
rsigma_dot(  7,  7,  2,  2 ) = -7.275720;
rsigma_dot(  7,  8,  2,  2 ) = -4.746130;
rsigma_dot(  7,  9,  2,  2 ) = -4.145350;
rsigma_dot(  7, 10,  2,  2 ) = -4.923060;
rsigma_dot(  7, 11,  2,  2 ) = -6.725270;
rsigma_dot(  7, 12,  2,  2 ) = -5.620700;
rsigma_dot(  7, 13,  2,  2 ) = -2.611400;
rsigma_dot(  7, 14,  2,  2 ) = -0.213500;
rsigma_dot(  7, 15,  2,  2 ) = -0.213500;
rsigma_dot(  7, 16,  2,  2 ) = -0.213500;
rsigma_dot(  7, 17,  2,  2 ) = -0.213500;
rsigma_dot(  7, 18,  2,  2 ) = -0.213500;
rsigma_dot(  8,  1,  2,  2 ) = -0.167000;
rsigma_dot(  8,  2,  2,  2 ) = -0.167000;
rsigma_dot(  8,  3,  2,  2 ) = -0.167000;
rsigma_dot(  8,  4,  2,  2 ) = -1.265600;
rsigma_dot(  8,  5,  2,  2 ) = -6.383610;
rsigma_dot(  8,  6,  2,  2 ) = -9.098560;
rsigma_dot(  8,  7,  2,  2 ) = -6.418910;
rsigma_dot(  8,  8,  2,  2 ) = -4.561460;
rsigma_dot(  8,  9,  2,  2 ) = -4.174340;
rsigma_dot(  8, 10,  2,  2 ) = -4.979190;
rsigma_dot(  8, 11,  2,  2 ) = -7.462060;
rsigma_dot(  8, 12,  2,  2 ) = -7.337890;
rsigma_dot(  8, 13,  2,  2 ) = -3.722400;
rsigma_dot(  8, 14,  2,  2 ) = -0.167000;
rsigma_dot(  8, 15,  2,  2 ) = -0.167000;
rsigma_dot(  8, 16,  2,  2 ) = -0.167000;
rsigma_dot(  8, 17,  2,  2 ) = -0.167000;
rsigma_dot(  8, 18,  2,  2 ) = -0.167000;
rsigma_dot(  9,  1,  2,  2 ) = -0.122600;
rsigma_dot(  9,  2,  2,  2 ) = -0.122600;
rsigma_dot(  9,  3,  2,  2 ) = -0.122600;
rsigma_dot(  9,  4,  2,  2 ) = -1.221200;
rsigma_dot(  9,  5,  2,  2 ) = -7.827370;
rsigma_dot(  9,  6,  2,  2 ) = -8.067400;
rsigma_dot(  9,  7,  2,  2 ) = -5.809530;
rsigma_dot(  9,  8,  2,  2 ) = -4.966740;
rsigma_dot(  9,  9,  2,  2 ) = -4.327250;
rsigma_dot(  9, 10,  2,  2 ) = -4.884730;
rsigma_dot(  9, 11,  2,  2 ) = -7.072410;
rsigma_dot(  9, 12,  2,  2 ) = -8.543900;
rsigma_dot(  9, 13,  2,  2 ) = -5.252450;
rsigma_dot(  9, 14,  2,  2 ) = -1.732000;
rsigma_dot(  9, 15,  2,  2 ) = -0.122600;
rsigma_dot(  9, 16,  2,  2 ) = -0.122600;
rsigma_dot(  9, 17,  2,  2 ) = -0.122600;
rsigma_dot(  9, 18,  2,  2 ) = -0.122600;
rsigma_dot( 10,  1,  2,  2 ) = -0.080000;
rsigma_dot( 10,  2,  2,  2 ) = -0.080000;
rsigma_dot( 10,  3,  2,  2 ) = -1.689400;
rsigma_dot( 10,  4,  2,  2 ) = -2.913200;
rsigma_dot( 10,  5,  2,  2 ) = -7.579420;
rsigma_dot( 10,  6,  2,  2 ) = -6.907620;
rsigma_dot( 10,  7,  2,  2 ) = -5.422330;
rsigma_dot( 10,  8,  2,  2 ) = -5.161400;
rsigma_dot( 10,  9,  2,  2 ) = -4.590850;
rsigma_dot( 10, 10,  2,  2 ) = -4.970340;
rsigma_dot( 10, 11,  2,  2 ) = -6.331900;
rsigma_dot( 10, 12,  2,  2 ) = -8.504850;
rsigma_dot( 10, 13,  2,  2 ) = -6.881280;
rsigma_dot( 10, 14,  2,  2 ) = -2.644900;
rsigma_dot( 10, 15,  2,  2 ) = -0.080000;
rsigma_dot( 10, 16,  2,  2 ) = -1.178600;
rsigma_dot( 10, 17,  2,  2 ) = -0.080000;
rsigma_dot( 10, 18,  2,  2 ) = -0.080000;
rsigma_dot( 11,  1,  2,  2 ) = -0.039200;
rsigma_dot( 11,  2,  2,  2 ) = -1.137800;
rsigma_dot( 11,  3,  2,  2 ) = -2.437100;
rsigma_dot( 11,  4,  2,  2 ) = -4.356660;
rsigma_dot( 11,  5,  2,  2 ) = -6.159470;
rsigma_dot( 11,  6,  2,  2 ) = -6.109910;
rsigma_dot( 11,  7,  2,  2 ) = -5.082600;
rsigma_dot( 11,  8,  2,  2 ) = -5.203960;
rsigma_dot( 11,  9,  2,  2 ) = -5.029610;
rsigma_dot( 11, 10,  2,  2 ) = -5.180840;
rsigma_dot( 11, 11,  2,  2 ) = -5.633880;
rsigma_dot( 11, 12,  2,  2 ) = -7.850740;
rsigma_dot( 11, 13,  2,  2 ) = -8.236990;
rsigma_dot( 11, 14,  2,  2 ) = -3.258000;
rsigma_dot( 11, 15,  2,  2 ) = -2.604100;
rsigma_dot( 11, 16,  2,  2 ) = -1.648600;
rsigma_dot( 11, 17,  2,  2 ) = -0.039200;
rsigma_dot( 11, 18,  2,  2 ) = -0.039200;
rsigma_dot( 12,  1,  2,  2 ) = 0.000000;
rsigma_dot( 12,  2,  2,  2 ) = -2.197200;
rsigma_dot( 12,  3,  2,  2 ) = -2.708000;
rsigma_dot( 12,  4,  2,  2 ) = -5.093700;
rsigma_dot( 12,  5,  2,  2 ) = -5.802070;
rsigma_dot( 12,  6,  2,  2 ) = -5.587200;
rsigma_dot( 12,  7,  2,  2 ) = -5.198450;
rsigma_dot( 12,  8,  2,  2 ) = -5.459540;
rsigma_dot( 12,  9,  2,  2 ) = -5.068860;
rsigma_dot( 12, 10,  2,  2 ) = -5.272950;
rsigma_dot( 12, 11,  2,  2 ) = -5.370590;
rsigma_dot( 12, 12,  2,  2 ) = -6.597100;
rsigma_dot( 12, 13,  2,  2 ) = -7.631870;
rsigma_dot( 12, 14,  2,  2 ) = -5.043380;
rsigma_dot( 12, 15,  2,  2 ) = -3.218800;
rsigma_dot( 12, 16,  2,  2 ) = -2.833200;
rsigma_dot( 12, 17,  2,  2 ) = -2.197200;
rsigma_dot( 12, 18,  2,  2 ) = 0.000000;

// section 21
rsigma_dot(  1,  1,  2,  1 ) = -0.974100;
rsigma_dot(  1,  2,  2,  1 ) = -0.974100;
rsigma_dot(  1,  3,  2,  1 ) = -0.974100;
rsigma_dot(  1,  4,  2,  1 ) = -0.974100;
rsigma_dot(  1,  5,  2,  1 ) = -0.974100;
rsigma_dot(  1,  6,  2,  1 ) = -0.974100;
rsigma_dot(  1,  7,  2,  1 ) = -0.974100;
rsigma_dot(  1,  8,  2,  1 ) = -2.920000;
rsigma_dot(  1,  9,  2,  1 ) = -3.171400;
rsigma_dot(  1, 10,  2,  1 ) = -0.974100;
rsigma_dot(  1, 11,  2,  1 ) = -2.583600;
rsigma_dot(  1, 12,  2,  1 ) = -0.974100;
rsigma_dot(  1, 13,  2,  1 ) = -0.974100;
rsigma_dot(  1, 14,  2,  1 ) = -0.974100;
rsigma_dot(  1, 15,  2,  1 ) = -0.974100;
rsigma_dot(  1, 16,  2,  1 ) = -0.974100;
rsigma_dot(  1, 17,  2,  1 ) = -0.974100;
rsigma_dot(  1, 18,  2,  1 ) = -0.974100;
rsigma_dot(  2,  1,  2,  1 ) = -0.909600;
rsigma_dot(  2,  2,  2,  1 ) = -0.909600;
rsigma_dot(  2,  3,  2,  1 ) = -0.909600;
rsigma_dot(  2,  4,  2,  1 ) = -0.909600;
rsigma_dot(  2,  5,  2,  1 ) = -0.909600;
rsigma_dot(  2,  6,  2,  1 ) = -0.909600;
rsigma_dot(  2,  7,  2,  1 ) = -0.909600;
rsigma_dot(  2,  8,  2,  1 ) = -4.343580;
rsigma_dot(  2,  9,  2,  1 ) = -6.015540;
rsigma_dot(  2, 10,  2,  1 ) = -4.205430;
rsigma_dot(  2, 11,  2,  1 ) = -2.008200;
rsigma_dot(  2, 12,  2,  1 ) = -2.519000;
rsigma_dot(  2, 13,  2,  1 ) = -0.909600;
rsigma_dot(  2, 14,  2,  1 ) = -0.909600;
rsigma_dot(  2, 15,  2,  1 ) = -0.909600;
rsigma_dot(  2, 16,  2,  1 ) = -0.909600;
rsigma_dot(  2, 17,  2,  1 ) = -0.909600;
rsigma_dot(  2, 18,  2,  1 ) = -0.909600;
rsigma_dot(  3,  1,  2,  1 ) = -0.849000;
rsigma_dot(  3,  2,  2,  1 ) = -0.849000;
rsigma_dot(  3,  3,  2,  1 ) = -0.849000;
rsigma_dot(  3,  4,  2,  1 ) = -0.849000;
rsigma_dot(  3,  5,  2,  1 ) = -0.849000;
rsigma_dot(  3,  6,  2,  1 ) = -0.849000;
rsigma_dot(  3,  7,  2,  1 ) = -3.046200;
rsigma_dot(  3,  8,  2,  1 ) = -5.990630;
rsigma_dot(  3,  9,  2,  1 ) = -7.627760;
rsigma_dot(  3, 10,  2,  1 ) = -6.013760;
rsigma_dot(  3, 11,  2,  1 ) = -2.794900;
rsigma_dot(  3, 12,  2,  1 ) = -2.458400;
rsigma_dot(  3, 13,  2,  1 ) = -0.849000;
rsigma_dot(  3, 14,  2,  1 ) = -0.849000;
rsigma_dot(  3, 15,  2,  1 ) = -0.849000;
rsigma_dot(  3, 16,  2,  1 ) = -0.849000;
rsigma_dot(  3, 17,  2,  1 ) = -0.849000;
rsigma_dot(  3, 18,  2,  1 ) = -0.849000;
rsigma_dot(  4,  1,  2,  1 ) = -0.791800;
rsigma_dot(  4,  2,  2,  1 ) = -0.791800;
rsigma_dot(  4,  3,  2,  1 ) = -0.791800;
rsigma_dot(  4,  4,  2,  1 ) = -0.791800;
rsigma_dot(  4,  5,  2,  1 ) = -0.791800;
rsigma_dot(  4,  6,  2,  1 ) = -2.401300;
rsigma_dot(  4,  7,  2,  1 ) = -3.736250;
rsigma_dot(  4,  8,  2,  1 ) = -6.512120;
rsigma_dot(  4,  9,  2,  1 ) = -8.856760;
rsigma_dot(  4, 10,  2,  1 ) = -7.235940;
rsigma_dot(  4, 11,  2,  1 ) = -4.010690;
rsigma_dot(  4, 12,  2,  1 ) = -2.401300;
rsigma_dot(  4, 13,  2,  1 ) = -1.890400;
rsigma_dot(  4, 14,  2,  1 ) = -0.791800;
rsigma_dot(  4, 15,  2,  1 ) = -0.791800;
rsigma_dot(  4, 16,  2,  1 ) = -0.791800;
rsigma_dot(  4, 17,  2,  1 ) = -0.791800;
rsigma_dot(  4, 18,  2,  1 ) = -0.791800;
rsigma_dot(  5,  1,  2,  1 ) = -0.737700;
rsigma_dot(  5,  2,  2,  1 ) = -0.737700;
rsigma_dot(  5,  3,  2,  1 ) = -0.737700;
rsigma_dot(  5,  4,  2,  1 ) = -0.737700;
rsigma_dot(  5,  5,  2,  1 ) = -0.737700;
rsigma_dot(  5,  6,  2,  1 ) = -2.683700;
rsigma_dot(  5,  7,  2,  1 ) = -4.171730;
rsigma_dot(  5,  8,  2,  1 ) = -7.436010;
rsigma_dot(  5,  9,  2,  1 ) = -9.773380;
rsigma_dot(  5, 10,  2,  1 ) = -8.510500;
rsigma_dot(  5, 11,  2,  1 ) = -4.669570;
rsigma_dot(  5, 12,  2,  1 ) = -3.445800;
rsigma_dot(  5, 13,  2,  1 ) = -1.836400;
rsigma_dot(  5, 14,  2,  1 ) = -0.737700;
rsigma_dot(  5, 15,  2,  1 ) = -0.737700;
rsigma_dot(  5, 16,  2,  1 ) = -0.737700;
rsigma_dot(  5, 17,  2,  1 ) = -0.737700;
rsigma_dot(  5, 18,  2,  1 ) = -0.737700;
rsigma_dot(  6,  1,  2,  1 ) = -0.686500;
rsigma_dot(  6,  2,  2,  1 ) = -0.686500;
rsigma_dot(  6,  3,  2,  1 ) = -0.686500;
rsigma_dot(  6,  4,  2,  1 ) = -0.686500;
rsigma_dot(  6,  5,  2,  1 ) = -0.686500;
rsigma_dot(  6,  6,  2,  1 ) = -0.686500;
rsigma_dot(  6,  7,  2,  1 ) = -4.297370;
rsigma_dot(  6,  8,  2,  1 ) = -7.670240;
rsigma_dot(  6,  9,  2,  1 ) = -9.274480;
rsigma_dot(  6, 10,  2,  1 ) = -8.780530;
rsigma_dot(  6, 11,  2,  1 ) = -5.395980;
rsigma_dot(  6, 12,  2,  1 ) = -2.632400;
rsigma_dot(  6, 13,  2,  1 ) = -3.084300;
rsigma_dot(  6, 14,  2,  1 ) = -0.686500;
rsigma_dot(  6, 15,  2,  1 ) = -0.686500;
rsigma_dot(  6, 16,  2,  1 ) = -0.686500;
rsigma_dot(  6, 17,  2,  1 ) = -0.686500;
rsigma_dot(  6, 18,  2,  1 ) = -0.686500;
rsigma_dot(  7,  1,  2,  1 ) = -0.637700;
rsigma_dot(  7,  2,  2,  1 ) = -0.637700;
rsigma_dot(  7,  3,  2,  1 ) = -0.637700;
rsigma_dot(  7,  4,  2,  1 ) = -0.637700;
rsigma_dot(  7,  5,  2,  1 ) = -1.736300;
rsigma_dot(  7,  6,  2,  1 ) = -3.682200;
rsigma_dot(  7,  7,  2,  1 ) = -4.645000;
rsigma_dot(  7,  8,  2,  1 ) = -6.967380;
rsigma_dot(  7,  9,  2,  1 ) = -7.766160;
rsigma_dot(  7, 10,  2,  1 ) = -7.608390;
rsigma_dot(  7, 11,  2,  1 ) = -5.416790;
rsigma_dot(  7, 12,  2,  1 ) = -3.773160;
rsigma_dot(  7, 13,  2,  1 ) = -3.345700;
rsigma_dot(  7, 14,  2,  1 ) = -2.834900;
rsigma_dot(  7, 15,  2,  1 ) = -2.247100;
rsigma_dot(  7, 16,  2,  1 ) = -0.637700;
rsigma_dot(  7, 17,  2,  1 ) = -0.637700;
rsigma_dot(  7, 18,  2,  1 ) = -0.637700;
rsigma_dot(  8,  1,  2,  1 ) = -0.591100;
rsigma_dot(  8,  2,  2,  1 ) = -0.591100;
rsigma_dot(  8,  3,  2,  1 ) = -0.591100;
rsigma_dot(  8,  4,  2,  1 ) = -0.591100;
rsigma_dot(  8,  5,  2,  1 ) = -2.788400;
rsigma_dot(  8,  6,  2,  1 ) = -4.087650;
rsigma_dot(  8,  7,  2,  1 ) = -4.702020;
rsigma_dot(  8,  8,  2,  1 ) = -6.084200;
rsigma_dot(  8,  9,  2,  1 ) = -6.877140;
rsigma_dot(  8, 10,  2,  1 ) = -6.827510;
rsigma_dot(  8, 11,  2,  1 ) = -5.980210;
rsigma_dot(  8, 12,  2,  1 ) = -4.087650;
rsigma_dot(  8, 13,  2,  1 ) = -3.635700;
rsigma_dot(  8, 14,  2,  1 ) = -4.202060;
rsigma_dot(  8, 15,  2,  1 ) = -2.537100;
rsigma_dot(  8, 16,  2,  1 ) = -0.591100;
rsigma_dot(  8, 17,  2,  1 ) = -0.591100;
rsigma_dot(  8, 18,  2,  1 ) = -0.591100;
rsigma_dot(  9,  1,  2,  1 ) = -0.546700;
rsigma_dot(  9,  2,  2,  1 ) = -0.546700;
rsigma_dot(  9,  3,  2,  1 ) = -0.546700;
rsigma_dot(  9,  4,  2,  1 ) = -0.546700;
rsigma_dot(  9,  5,  2,  1 ) = -3.254700;
rsigma_dot(  9,  6,  2,  1 ) = -4.102040;
rsigma_dot(  9,  7,  2,  1 ) = -4.916140;
rsigma_dot(  9,  8,  2,  1 ) = -5.788440;
rsigma_dot(  9,  9,  2,  1 ) = -5.777800;
rsigma_dot(  9, 10,  2,  1 ) = -6.126420;
rsigma_dot(  9, 11,  2,  1 ) = -5.997730;
rsigma_dot(  9, 12,  2,  1 ) = -4.721080;
rsigma_dot(  9, 13,  2,  1 ) = -4.102040;
rsigma_dot(  9, 14,  2,  1 ) = -4.157610;
rsigma_dot(  9, 15,  2,  1 ) = -3.111600;
rsigma_dot(  9, 16,  2,  1 ) = -1.645300;
rsigma_dot(  9, 17,  2,  1 ) = -0.546700;
rsigma_dot(  9, 18,  2,  1 ) = -0.546700;
rsigma_dot( 10,  1,  2,  1 ) = -0.504100;
rsigma_dot( 10,  2,  2,  1 ) = -0.504100;
rsigma_dot( 10,  3,  2,  1 ) = -0.504100;
rsigma_dot( 10,  4,  2,  1 ) = -0.504100;
rsigma_dot( 10,  5,  2,  1 ) = -3.548700;
rsigma_dot( 10,  6,  2,  1 ) = -4.873580;
rsigma_dot( 10,  7,  2,  1 ) = -5.176960;
rsigma_dot( 10,  8,  2,  1 ) = -5.756400;
rsigma_dot( 10,  9,  2,  1 ) = -5.573030;
rsigma_dot( 10, 10,  2,  1 ) = -5.702630;
rsigma_dot( 10, 11,  2,  1 ) = -5.668920;
rsigma_dot( 10, 12,  2,  1 ) = -5.078840;
rsigma_dot( 10, 13,  2,  1 ) = -4.678520;
rsigma_dot( 10, 14,  2,  1 ) = -4.310790;
rsigma_dot( 10, 15,  2,  1 ) = -4.217700;
rsigma_dot( 10, 16,  2,  1 ) = -2.113600;
rsigma_dot( 10, 17,  2,  1 ) = -0.504100;
rsigma_dot( 10, 18,  2,  1 ) = -0.504100;
rsigma_dot( 11,  1,  2,  1 ) = -0.463300;
rsigma_dot( 11,  2,  2,  1 ) = -0.463300;
rsigma_dot( 11,  3,  2,  1 ) = -0.463300;
rsigma_dot( 11,  4,  2,  1 ) = -2.409200;
rsigma_dot( 11,  5,  2,  1 ) = -4.313460;
rsigma_dot( 11,  6,  2,  1 ) = -5.225480;
rsigma_dot( 11,  7,  2,  1 ) = -5.493750;
rsigma_dot( 11,  8,  2,  1 ) = -5.353660;
rsigma_dot( 11,  9,  2,  1 ) = -5.154660;
rsigma_dot( 11, 10,  2,  1 ) = -5.397780;
rsigma_dot( 11, 11,  2,  1 ) = -5.726000;
rsigma_dot( 11, 12,  2,  1 ) = -5.480590;
rsigma_dot( 11, 13,  2,  1 ) = -4.018660;
rsigma_dot( 11, 14,  2,  1 ) = -4.313460;
rsigma_dot( 11, 15,  2,  1 ) = -4.574180;
rsigma_dot( 11, 16,  2,  1 ) = -2.861200;
rsigma_dot( 11, 17,  2,  1 ) = -0.463300;
rsigma_dot( 11, 18,  2,  1 ) = -1.561900;
rsigma_dot( 12,  1,  2,  1 ) = -0.424100;
rsigma_dot( 12,  2,  2,  1 ) = -0.424100;
rsigma_dot( 12,  3,  2,  1 ) = -2.822000;
rsigma_dot( 12,  4,  2,  1 ) = -3.132100;
rsigma_dot( 12,  5,  2,  1 ) = -5.115440;
rsigma_dot( 12,  6,  2,  1 ) = -5.697090;
rsigma_dot( 12,  7,  2,  1 ) = -5.747100;
rsigma_dot( 12,  8,  2,  1 ) = -5.386930;
rsigma_dot( 12,  9,  2,  1 ) = -5.019210;
rsigma_dot( 12, 10,  2,  1 ) = -5.186260;
rsigma_dot( 12, 11,  2,  1 ) = -5.775950;
rsigma_dot( 12, 12,  2,  1 ) = -5.600240;
rsigma_dot( 12, 13,  2,  1 ) = -4.394380;
rsigma_dot( 12, 14,  2,  1 ) = -3.257300;
rsigma_dot( 12, 15,  2,  1 ) = -3.920600;
rsigma_dot( 12, 16,  2,  1 ) = -3.643000;
rsigma_dot( 12, 17,  2,  1 ) = -1.522700;
rsigma_dot( 12, 18,  2,  1 ) = -0.424100;
}*/


} // ns sspot
} // ns potentials
} // ns fldsgn
} // ns protocols
