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
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit header
#include <protocols/fldsgn/potentials/sspot/HSPairPotential.hh>

// Package headers
#include <protocols/fldsgn/potentials/sspot/util.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/fldsgn/topology/BB_Pos.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>

// utility
#include <utility/io/izstream.hh>

// C++ headers
#include <cmath>
#include <iostream>

#include <ObjexxFCL/format.hh>

#include <utility/vector1.hh>
#include <ObjexxFCL/FArray3D.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>


static THREAD_LOCAL basic::Tracer TR( "protocols.fldsgn.potentials.sspot.HSPairPotential", basic::t_info );

namespace protocols {
namespace fldsgn {
namespace potentials {
namespace sspot {

/// @brief default constructor
HSPairPotential::HSPairPotential():
	dist_cutoff_( 12.0 ),
	phithetascore_( 2, 3, 36, 36 )
{
	load_phi_theta_bins();
}

/// @brief default destructor
HSPairPotential::~HSPairPotential()
{}

/// @brief
HSPairPotential::Real
HSPairPotential::calc_phithetascore( Size const strand_seqsep, Real const phi, Real const theta ) const
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
	return phithetascore_( 1, int( istrand_seqsep ), int( iphi ), int( itheta ) );
}


/// @brief function reads in two points in sequence and returns two points in space,
/// @brief the endpoints of the axis through an alpha-helix
void
HSPairPotential::helix_end(
	Size const & pos1,
	BB_Pos const & bb_pos,
	Vector & p1,
	Vector & p2
) const
{
	static Real const eleven_inv = 1.0 / 11.0;
	Size const s1 = pos1;
	Size const s2 = pos1+1;
	Size const s3 = pos1+2;
	Size const s4 = pos1+3;

	Vector const Epos_sum( (                  bb_pos.CA( s1 ) + bb_pos.C( s1 ) ) +
		( bb_pos.N( s2 ) + bb_pos.CA( s2 ) + bb_pos.C( s2 ) ) +
		( bb_pos.N( s3 ) + bb_pos.CA( s3 ) + bb_pos.C( s3 ) ) +
		( bb_pos.N( s4 ) + bb_pos.CA( s4 )                  ) );

	p1 = ( Epos_sum + bb_pos.N( s1 ) ) * eleven_inv;
	p2 = ( Epos_sum + bb_pos.C( s4 ) ) * eleven_inv;
}


/// @brief
void
HSPairPotential::score(
	Pose const & pose,
	SS_Info2 const & ss_info,
	Real & hs_score
) const
{
	using protocols::fldsgn::potentials::sspot::get_foldtree_seqsep;
	using protocols::fldsgn::potentials::sspot::spherical;
	using protocols::fldsgn::topology::Helices;
	using protocols::fldsgn::topology::Strands;
	using protocols::fldsgn::topology::StrandOP;
	using protocols::fldsgn::topology::HelixOP;

	hs_score = 0.0;

	Helices const & helices( ss_info.helices() );
	Strands const & strands( ss_info.strands() );
	BB_Pos const & bb_pos( ss_info.bb_pos() );

	EnergyGraph const & energy_graph( pose.energies().energy_graph() );

	for ( Size ihelix=1; ihelix<=helices.size(); ihelix++ ) {

		HelixOP const helix( helices[ ihelix ] );

		if ( helix->length() <= 3 ) continue;

		for ( Size ss1=helix->begin(); ss1<=helix->end()-3; ss1++ ) {

			Vector pt1, pt2;
			helix_end( ss1, bb_pos, pt1, pt2 );

			for ( core::graph::Graph::EdgeListConstIter
					iru  = energy_graph.get_node( ss1+1 )->const_edge_list_begin(),
					irue = energy_graph.get_node( ss1+1 )->const_edge_list_end();
					iru != irue; ++iru ) {

				Size ss2( (*iru)->get_second_node_ind() );
				//Edges always have first node < second node. Just in case we picked the wrong one:
				if ( ss1+1 == ss2 ) ss2 = (*iru)->get_first_node_ind();

				Size jstrand = ss_info.strand_id( ss2 );
				if ( jstrand == 0 || ss_info.strand_id( ss2+1 ) == 0 ) continue;
				if ( pose.residue_type( ss2 ).is_upper_terminus() ) continue;

				Vector const & pt3( bb_pos.N( ss2   ) );
				Vector const & pt4( bb_pos.C( ss2+1 ) );

				// midpoint coordinates of two dimers
				Vector cen1 = Real( 0.5 )*( pt1 + pt2 );
				Vector cen2 = Real( 0.5 )*( pt3 + pt4 );

				// vector between midpoints
				Vector mid_vector = cen2 -cen1;

				if ( mid_vector.length() <= dist_cutoff_ ) {
					Real phi, theta;
					spherical( pt2, pt4, phi, theta, cen1, cen2, mid_vector );

					Size const helix_end_1 = helix->begin() - 1;
					Size const helix_end_2 = helix->end() + 1;

					Size const strand_end_1 = strands[ jstrand ]->begin() - 1;
					Size const strand_end_2 = strands[ jstrand ]->end() + 1;

					Size seqsep = std::min( get_foldtree_seqsep( pose, helix_end_2, strand_end_1 ) + 1,
						get_foldtree_seqsep( pose, strand_end_2, helix_end_1 ) + 1 );

					hs_score += calc_phithetascore( seqsep, phi, theta );
				}
			} // loop over neighbors of ss1
		} // ss1
	} // ihelix

	// modify by proper weighting
	hs_score *= 0.090;
} // score


/// @brief load phi/theta bins for use in secondary structure scoring
void
HSPairPotential::load_phi_theta_bins( String const & hs_filename )
{
	using ObjexxFCL::format::skip;
	typedef ObjexxFCL::FArray3D< Real > FArray3D_real;
	FArray3D_real pts_HS( 36, 36, 3 );

	FArray1D_int iptsn( 36 );
	for ( int itheta = 1; itheta <= 36; ++itheta ) {
		iptsn(itheta) = 100;
	}

	utility::io::izstream HS_stream;
	basic::database::open( HS_stream, hs_filename );
	for ( int isep = 1; isep <= 3; ++isep ) {
		for ( int itheta = 1; itheta <= 36; ++itheta ) {
			for ( int iph = 1; iph <= 36; ++iph ) {
				HS_stream >> pts_HS( itheta, iph, isep ) >> skip;
			}
		}
	}
	HS_stream.close();
	HS_stream.clear();

	for ( int isep = 1; isep <= 3; ++isep ) {
		Real tot = 0.0;
		Real totn = 0.0;
		for ( int iph = 1; iph <= 36; ++iph ) {
			for ( int itheta = 1; itheta <= 36; ++itheta ) {
				//  SMALL COUNTS CORRECTION
				pts_HS( itheta, iph, isep ) += iptsn( itheta )*0.000001f;
				tot += pts_HS( itheta, iph, isep );
				totn += iptsn( itheta );
			}
		}
		for ( int iph = 1; iph <= 36; ++iph ) {
			for ( int itheta = 1; itheta <= 36; ++itheta ) {
				phithetascore_( 1, isep, iph, itheta ) = -std::log(pts_HS( itheta, iph, isep )/tot ) +
					std::log(iptsn( itheta )/totn );
			}
		}
	}
} // load_phi_theta_bins


} // ns sspot
} // ns potentials
} // ns fldsgn
} // ns protocols
