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

#include <core/scoring/dna/setup.hh>
#include <core/scoring/dna/base_geometry.hh>
#include <core/scoring/dna/BasePartner.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace dna {

static thread_local basic::Tracer TR( "core.scoring.dna.setup" );

void
set_base_partner(
	pose::Pose & pose
)
{
	//using core::pose::datacache::CacheableDataType::BASE_PARTNER;
	using basic::datacache::DataCache_CacheableData;
	utility::vector1< Size > partner;
	find_basepairs( pose, partner );
	pose.data().set( core::pose::datacache::CacheableDataType::BASE_PARTNER, DataCache_CacheableData::DataOP( new BasePartner( partner ) ) );
}

void
find_basepairs(
	pose::Pose const & pose,
	//utility::vector1< std::pair< int, int > > & pairs,
	utility::vector1< Size > & partner
)
{
	using conformation::Residue;
	using namespace chemical;

	Real const max_d( 4.0 );
	Size const nres( pose.total_residue() );

	//pairs.clear();
	partner.clear();
	partner.resize( nres, 0 );

	std::map< AA, AA > base_partner;
	base_partner[ na_ade ] = na_thy;
	base_partner[ na_thy ] = na_ade;
	base_partner[ na_gua ] = na_cyt;
	base_partner[ na_cyt ] = na_gua;

	std::map< AA, std::string > hbond_atom;
	hbond_atom[ na_ade ] = "N1";
	hbond_atom[ na_thy ] = "N3";
	hbond_atom[ na_gua ] = "N1";
	hbond_atom[ na_cyt ] = "N3";

	for ( Size i=1; i<= nres; ++i ) {
		Residue const & i_rsd( pose.residue(i) );
		AA const & i_aa( i_rsd.aa() );
		if ( i_rsd.is_DNA() ) {
			// hbond atom, base y-axis
			Vector const & i_xyz( i_rsd.xyz( hbond_atom[ i_aa ] ) );
			Vector const i_axis( get_y_axis( i_rsd, 1 /*strand*/ ) );

			Real best(1000.0);
			for ( Size j=1; j<= nres; ++j ) {
				if ( j==i ) continue;
				Residue const & j_rsd( pose.residue( j ) );
				AA const & j_aa( j_rsd.aa() );
				if ( j_aa == base_partner[ i_aa ] ) {
					Vector const & j_xyz( j_rsd.xyz( hbond_atom[ j_aa ] ) );
					Real d( i_xyz.distance( j_xyz ) );
					if ( d<max_d ) {
						Vector const j_axis( get_y_axis( j_rsd, 2 /*strand*/ ) );
						Vector const u( ( i_xyz - j_xyz ).normalized() );
						Real const dot1( dot( i_axis, j_axis ) );
						Real const dot2( dot( i_axis, u ) );
						Real const dot3( dot( j_axis, u ) );
						d -= dot1+dot2+dot3;
						if ( d<best && dot1 > 0.75 && dot2 > 0.75 && dot3 > 0.75 ) {
							best = d;
							partner[i] = j;
						}
					}
				}
			}
		}
	}

	for ( Size i=1; i<= nres; ++i ) {
		if ( !partner[i] ) continue;

		if ( partner[ partner[i] ] != i ) {
			partner[i] = 0;
			continue;
		}

		if ( i < partner[i] ) {
			//pairs.push_back( std::make_pair( i, partner[i] ) );
			TR(basic::t_debug) << "found basepair: " << i << ' ' << partner[i] << std::endl;
		}
	}
}


} // namespace dna
}} // scoring core
