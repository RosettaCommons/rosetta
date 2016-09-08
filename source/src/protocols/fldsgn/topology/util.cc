// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/fldsgn/topology/util.cc
/// @brief utilities for fldsgn/topology
/// @author Nobuyasu Koga ( nobuyasu@u.washington.edu )

// unit header
#include <protocols/fldsgn/topology/util.hh>
#include <protocols/fldsgn/topology/StrandPairing.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/forge/build/Interval.hh>

// Project Headers
#include <core/scoring/Energies.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/dssp/StrandPairing.hh>
#include <core/id/AtomID_Map.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/sasa.hh>

// utility headers
#include <utility/exit.hh>

// C++ headers
#include <utility/assert.hh>
#include <iostream>
#include <sstream>
#include <basic/Tracer.hh>
#include <map>

#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>
static THREAD_LOCAL basic::Tracer TR( "protocols.topology.util" );

using namespace core;
typedef std::string String;
typedef utility::vector1< Size > VecSize;

namespace protocols {
namespace fldsgn {
namespace topology {

/// @brief convert StrandParingSet of dssp to fldsgn::topology::StrandPairingSet
protocols::fldsgn::topology::StrandPairingSet
calc_strand_pairing_set(
	core::pose::Pose const & pose,
	protocols::fldsgn::topology::SS_Info2_COP const ssinfo,
	core::Size minimum_pair_length )
{
	using core::Size;

	core::scoring::dssp::Dssp dssp( pose );

	std::map< String, StrandPairingOP > newpairs;

	core::scoring::dssp::StrandPairingSet spairset;
	spairset = dssp.strand_pairing_set();

	for ( Size ispair=1; ispair<=spairset.size(); ispair++ ) {

		core::scoring::dssp::StrandPairing sp;
		sp = spairset.strand_pairing( ispair );
		Size begin1 ( sp.begin1() );
		Size end1 ( sp.end1() );

		for ( Size iaa=begin1; iaa<=end1; ++iaa ) {

			Size istrand ( ssinfo->strand_id( iaa ) );
			if ( istrand == 0 ) continue;
			Size ist_begin ( ssinfo->strand( istrand )->begin() );

			Size jaa ( sp.get_pair( iaa ) );
			if ( jaa == 0 ) continue;
			Size pleats ( sp.get_pleating( iaa ) );
			Size jstrand ( ssinfo->strand_id( jaa ) );
			if ( jstrand == 0 ) continue;

			Size jst_begin ( ssinfo->strand( jstrand )->begin() );
			Size jst_end ( ssinfo->strand( jstrand )->end() );
			Size jst_length ( jst_end - jst_begin );

			if ( istrand == jstrand ) continue;

			if ( istrand !=0 && jstrand !=0 && jaa !=0  ) {

				char orient;
				Real rgstr_shift;

				if ( sp.antiparallel() ) {
					orient = 'A';
					rgstr_shift = Real(iaa) - Real(ist_begin) - (Real(jst_length) - (Real(jaa) - Real(jst_begin)));
				} else {
					orient = 'P';
					rgstr_shift = Real(iaa) - Real(ist_begin) - (Real(jaa) - Real(jst_begin));
				}

				std::ostringstream spairname;
				spairname << istrand << "-" << jstrand << "." << orient ;
				std::map<String, StrandPairingOP>::iterator it( newpairs.find( spairname.str() ) );
				if ( it == newpairs.end() ) {
					StrandPairingOP strand_pair( new StrandPairing( istrand, jstrand, iaa, jaa, pleats, rgstr_shift, orient ) );
					newpairs.insert( std::map<String, StrandPairingOP>::value_type( spairname.str(), strand_pair ) );
				} else {
					(*it).second->elongate( iaa, jaa, pleats, pleats );
				}

			} // istrand !=0 && jstrand !=0 && jaa !=0
		} // iaa
	} // ispair

	StrandPairingSet spairset_new;
	std::map<String, StrandPairingOP>::iterator it( newpairs.begin() );
	while ( it != newpairs.end() ) {

		// skip if pair length < minimum_pair_length
		if ( (*it).second->size1() < minimum_pair_length || (*it).second->size2() < minimum_pair_length ) {
			++it;
			continue;
		}


		TR.Debug << "Defined strand pair : " << *((*it).second)
			<< ", 1st_strand begin: " << (*it).second->begin1() << ", end: " <<  (*it).second->end1()
			<< ", 2nd_strand begin: " << (*it).second->begin2() << ", end: " <<  (*it).second->end2() << std::endl;
		spairset_new.push_back( (*it).second );
		++it;
	}

	spairset_new.finalize();

	return spairset_new;

} // clac_strand_pairings


/// @brief calc delta sasa, when a molecule is splited to 2parts.
core::Real
calc_delta_sasa(
	core::pose::Pose const & pose,
	utility::vector1< protocols::forge::build::Interval > intervals,
	Real const pore_radius )
{

	/// calc surface areas of total residues

	// define atom_map for main-chain and CB
	core::id::AtomID_Map< bool > atom_map;
	core::pose::initialize_atomid_map( atom_map, pose, false );
	for ( Size ir = 1; ir <= pose.size(); ++ir ) {
		for ( Size j = 1; j<=5; ++j ) {
			core::id::AtomID atom( j, ir );
			atom_map.set( atom, true );
		}
	}

	utility::vector1< Real > rsd_sasa;
	core::id::AtomID_Map< Real > atom_sasa;
	core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, pore_radius, false, atom_map );

	/// calc surface areas of A, B regions which are dicretized by intervals.
	/// A is non-assigned regions by intervals and B is assigned regions by intervals

	core::id::AtomID_Map< bool > atom_map_A;
	core::id::AtomID_Map< bool > atom_map_B;
	utility::vector1< Real > rsd_sasa_A;
	utility::vector1< Real > rsd_sasa_B;
	core::pose::initialize_atomid_map( atom_map_A, pose, false );
	core::pose::initialize_atomid_map( atom_map_B, pose, false );

	utility::vector1< bool > position_A( pose.size(), true );
	for ( Size ii=1; ii<=intervals.size(); ii++ ) {
		for ( Size jj=intervals[ ii ].left; jj<=intervals[ ii ].right; ++jj ) {
			position_A[ jj ] = false;
		}
		for ( Size j=1; j<=5; ++j ) {
			core::id::AtomID atom1( j, intervals[ ii ].left );
			atom_map_A.set( atom1, true );

			core::id::AtomID atom2( j, intervals[ ii ].right );
			atom_map_A.set( atom2, true );

			core::id::AtomID atom3( j, intervals[ ii ].left-1 );
			atom_map_B.set( atom3, true );

			core::id::AtomID atom4( j, intervals[ ii ].right+1 );
			atom_map_B.set( atom4, true );
		}
	}

	for ( Size jj=1; jj<=pose.size(); jj++ ) {
		if ( position_A[ jj ] ) {
			for ( Size j=1; j<=5; ++j ) {
				core::id::AtomID atom( j, jj );
				atom_map_A.set( atom, true );
			}
		} else {
			for ( Size j=1; j<=5; ++j ) {
				core::id::AtomID atom( j, jj );
				atom_map_B.set( atom, true );
			}
		}
	}

	core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa_A, pore_radius, false, atom_map_A );
	core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa_B, pore_radius, false, atom_map_B );

	Real tot_all( 0.0 ), tot_A( 0.0 ), tot_B( 0.0 );
	for ( Size ii=2; ii<=pose.size()-1; ii++ ) {
		tot_all += rsd_sasa[ ii ];
		if ( position_A[ ii ] ) {
			tot_A += rsd_sasa_A[ ii ];
		} else {
			tot_B += rsd_sasa_B[ ii ];
		}
	}

	//   for( Size ii=1; ii<=pose.size(); ii++ ) {
	//  if( position_A[ ii ] ) {
	//   std::cout << ii << " " << rsd_sasa[ ii ] << " " << rsd_sasa_A[ ii ] << " " << rsd_sasa_B[ ii ] << "*" <<std::endl;
	//  } else {
	//   std::cout << ii << " " << rsd_sasa[ ii ] << " " << rsd_sasa_A[ ii ] << " " << rsd_sasa_B[ ii ] << std::endl;
	//  }
	// }
	// std::cout << tot_all << " " << tot_A << " " << tot_B << std::endl;

	return tot_A + tot_B - tot_all;

} // calc_delta_sasa


/// @brief check kink of helix, return number of loosen hydrogen
core::Size
check_kink_helix(
	core::pose::Pose const & pose,
	core::Size const begin,
	core::Size const end )
{
	using core::scoring::EnergiesCacheableDataType::HBOND_SET;
	using core::scoring::hbonds::HBondSet;

	core::pose::Pose copy_pose( pose );
	HBondSet const & hbond_set( static_cast< HBondSet const & > ( copy_pose.energies().data().get( HBOND_SET )) );
	//hbond_set.show( copy_pose );

	Size num_broken_hbond( 0 );
	for ( Size ii=begin; ii<=end; ++ii ) {
		if ( ! hbond_set.acc_bbg_in_bb_bb_hbond( ii ) ) {
			TR << "hbonds of " << ii << " is broken. " << std::endl;
			num_broken_hbond++;
		}
	}

	return num_broken_hbond;
}


/// @brief check kink of helix, return number of loosen hydrogen
utility::vector1< core::scoring::hbonds::HBond >
check_internal_hbonds(
	core::pose::Pose const & pose,
	core::Size const begin,
	core::Size const end )
{
	using core::scoring::EnergiesCacheableDataType::HBOND_SET;
	using core::scoring::hbonds::HBondSet;
	using core::scoring::hbonds::HBond;

	core::pose::Pose copy_pose( pose );
	HBondSet const & hbond_set( static_cast< HBondSet const & > ( copy_pose.energies().data().get( HBOND_SET )) );

	utility::vector1< HBond > hbonds;
	for ( core::Size i=1; i<=(core::Size)hbond_set.nhbonds(); ++i ) {
		Size don_pos = hbond_set.hbond((int)i).don_res();
		Size acc_pos = hbond_set.hbond((int)i).acc_res();
		if ( don_pos >= begin && don_pos <= end &&
				acc_pos >= begin && acc_pos <= end ) {
			hbonds.push_back( hbond_set.hbond((int)i) );
		}
	}

	return hbonds;

}


// /// @brief
// utility::vector1< Size >
// split_into_ss_chunk(
// core::pose::Pose const & pose,
// protocols::fldsgn::topology::SS_Info2_COP const ssinfo )
//{
//
// /// types of chunk
// /// 1. sheet with short loops ( <=6 residues )
// /// 2. -loop-helix-loop-helix-
// /// 3. long loop ( >6 residues ) between strand & helix, strand & strand, helix & helix
//
//
// for() {
//
// }
//
//
// return;
//
//}


} // namespace topology
} // namespace fldsgn
} // namespace protocols
