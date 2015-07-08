// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ./src/protocols/fldsgn/topology/Sheet.cc
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@u.washington.edu )

// unit headers
#include <protocols/fldsgn/topology/Sheet.hh>

// Project headers
#include <basic/Tracer.hh>
#include <core/graph/DisjointSets.hh>
#include <protocols/fldsgn/topology/StrandPairing.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>

#include <core/id/AtomID_Map.hh>
#include <core/scoring/sasa.hh>
#include <core/pose/util.hh>

// utility headers
#include <utility/exit.hh>

// C++ headers
#include <iostream>
#include <map>

#include <utility/vector1.hh>
#include <utility/options/IntegerVectorOption.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>

#include <numeric/NumericTraits.hh>

static thread_local basic::Tracer TR( "protocols.topology.Sheet" );

using namespace core;

namespace protocols {
namespace fldsgn {
namespace topology {

/// @details Auto-generated virtual destructor
SheetSet::~SheetSet() {}


/// @brief default constructor
Sheet::Sheet() :
	num_strands_( 0 ),
	is_barrel_( false ),
	is_geometry_initialized_( false )
{
	initialize();
}

/// @brief value constructor
Sheet::Sheet( VecSize const & order_strands, VecInt const & orient_strands, bool is_barrel ):
	num_strands_( order_strands.size() ),
	is_barrel_( is_barrel ),
	order_strands_( order_strands ),
	orient_strands_( orient_strands ),
	is_geometry_initialized_( false )
{
	runtime_assert( order_strands.size() == orient_strands.size() );
	initialize();
}

/// @brief copy constructor
Sheet::Sheet( Sheet const & s ):
	ReferenceCount(),
	num_strands_( s.num_strands_ ),
	is_barrel_( s.is_barrel_ ),
	order_strands_( s.order_strands_ ),
	orient_strands_( s.orient_strands_ ),
	strand_order_( s.strand_order_ ),
	sheet_plane_( s.sheet_plane_ ),
	sheet_center_( s.sheet_center_ ),
	ca_cb_orients_( s.ca_cb_orients_ ),
	is_geometry_initialized_( s.is_geometry_initialized_ )
{}

/// @brief initialize this class
void
Sheet::initialize()
{
	num_strands_ = order_strands_.size();
	for( Size i=1; i<=order_strands_.size(); i++ ){
		strand_order_.insert( std::map< Size, Size >::value_type( order_strands_[i], i ) );
	}
}

/// @brief default destructor
Sheet::~Sheet(){}


/// @brief return strand pairing
std::ostream & operator<<( std::ostream & out, const Sheet &s )
{
	utility::vector1<Size> const order( s.order_strands() );
	utility::vector1<Real> const orient( s.orient_strands() );

	for( Size i=1; i<= order.size(); i++ ){
		out << "#Strand " << i << ": " << order[ i ] << "," << orient[ i ] << std::endl;;
	}
	return out;
}

/// @brief whether the given strand number belongs to this sheet or not
bool
Sheet::is_member( Size const s )
{
	for( utility::vector1< Size >::const_iterator iter = order_strands_.begin(),
				 iter_end = order_strands_.end(); iter != iter_end; ++iter ) {
		if( s == *iter ){
			return true;
		}
	}
	return false;
}

/// @brief calc surface areas only with beta-sheet
utility::vector1< Real >
Sheet::calc_sasa_bothsides( Pose const & pose, SS_Info2_COP const ssinfo, Real pore_radius )
{
	runtime_assert( is_geometry_initialized_ );

	utility::vector1< Real > surfaces( 2, 0.0 );

	// define atom_map for main-chain and CB
	core::id::AtomID_Map< bool > atom_map;
	core::pose::initialize_atomid_map( atom_map, pose, false );
	protocols::fldsgn::topology::Strands const strands( ssinfo->strands() );
	for( Size ii=1; ii<=order_strands_.size(); ++ii ) {
		Size const st( order_strands_[ ii ] );

		Size begin( strands[ st ]->begin() );
		Size end( strands[ st ]->end() );

		if( begin == 1 ) begin = 2;
		if( end == pose.total_residue() ) end = pose.total_residue() - 1;

		for( Size ir=begin-1; ir<=end+1; ir++ ) {
			for ( Size j = 1; j<=5; ++j ) {
				core::id::AtomID atom( j, ir );
				atom_map.set( atom, true );
			}
		}
	}

	// calc sasa
	core::id::AtomID_Map< Real > atom_sasa;
	utility::vector1< Real > rsd_sasa;
	core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, pore_radius, false, atom_map );

	for ( Size iaa=2; iaa<=pose.total_residue()-1; iaa++ ) {

		if( ca_cb_orient( iaa ) == 1 ) {
			surfaces[ 1 ] += rsd_sasa[ iaa ];
		} else if( ca_cb_orient( iaa ) == -1 ) {
			surfaces[ 2 ] += rsd_sasa[ iaa ];
		}
		// std::cout << iaa << " " << ca_cb_orient(iaa) << " " << rsd_sasa[ iaa ] << std::endl;
	}
	// std::cout << surfaces[ 1 ] << " " << surfaces[ 2 ] << std::endl;

	return surfaces;
}

/// @brief calc geometry
void
Sheet::calc_geometry( SS_Info2_COP const ssinfo )
{
	runtime_assert( ssinfo->bbpos_is_set() );

	if( num_strands_ == 0 ) {
		return;
	}

	protocols::fldsgn::topology::Strands const strands( ssinfo->strands() );

	Vector const v1 = ( strands[ order_strands_[ num_strands_ ] ]->mid_pos() -
											strands[ order_strands_[ 1 ] ]->mid_pos() ).normalized();


	Vector v2( strands[ order_strands_[ 1 ] ]->orient() );
	for( Size ii=2; ii<=num_strands_; ii++ ) {
		Real sign( 1.0 );
		if( orient_strands_[ 1 ] != orient_strands_[ ii ] ) {
			sign = - 1.0;
		}
		v2 += sign*strands[ order_strands_[ ii ] ]->orient();
	}

	sheet_plane_ = v1.cross( v2.normalized() );

	sheet_center_.clear();
	for( Size ii=1; ii<=num_strands_; ii++ ) {
		sheet_center_ += strands[ ii ]->mid_pos();
	}
	sheet_center_ = sheet_center_/Real( num_strands_ );

	protocols::fldsgn::topology::BB_Pos const bb_pos( ssinfo->bb_pos() );
	ca_cb_orients_.resize( bb_pos.size() );
	for( Size ii=1; ii<=bb_pos.size(); ii++ ) {
		ca_cb_orients_[ ii ] = 0;
	}

	Real half_pi = numeric::NumericTraits<Real>::pi()/2.0;
	utility::vector1< Size > anchor;
	anchor.resize( strands.size(), 0 );
	for( Size ii=1; ii<=order_strands_.size(); ii++ ) {

		Size const st( order_strands_[ ii ] );
		Size const begin( strands[ st ]->begin() );
		Size const end( strands[ st ]->end() );
		anchor[ st ] = begin;

		Vector cacb = bb_pos.CB( begin ) - bb_pos.CA( begin );
		Real max_angle = fabs( half_pi - angle_of( sheet_plane_, cacb ) );

		for( Size jj=begin+1; jj<=end; jj++ ) {
			cacb = bb_pos.CB( jj ) - bb_pos.CA( jj );
			Real angle = fabs( half_pi - angle_of( sheet_plane_, cacb ) );
			if( max_angle < angle ) {
				anchor[ st ] = jj;
				max_angle = angle;
			}
		}
	}


	for( Size ii=1; ii<=order_strands_.size(); ++ii ) {

		Size const st( order_strands_[ ii ] );
		Size const begin( strands[ st ]->begin() );
		Size const end( strands[ st ]->end() );

		Vector const cacb = bb_pos.CB( anchor[ st ] ) - bb_pos.CA( anchor[ st ] );
		if( cacb.dot( sheet_plane_ ) > 0 ) {
			ca_cb_orients_[ anchor[ st ] ] = 1;
		} else {
			ca_cb_orients_[ anchor[ st ] ] = -1;
		}

		for( Size jj=anchor[ st ]; jj>=begin+1; --jj ) {
			Real dot = ( bb_pos.CB( jj ) - bb_pos.CA( jj ) ).dot( bb_pos.CB( jj-1 ) - bb_pos.CA( jj-1 ) );
			if( dot > 0 ) {
				ca_cb_orients_[ jj-1 ] = ca_cb_orients_[ jj ];
			} else {
				ca_cb_orients_[ jj-1 ] = -1 * ca_cb_orients_[ jj ];
			}
		}

		for( Size jj=anchor[ st ]; jj<=end-1; ++jj ) {
			Real dot = ( bb_pos.CB( jj ) - bb_pos.CA( jj ) ).dot( bb_pos.CB( jj+1 ) - bb_pos.CA( jj+1 ) );
			if( dot > 0 ) {
				ca_cb_orients_[ jj+1 ] = ca_cb_orients_[ jj ];
			} else {
				ca_cb_orients_[ jj+1 ] = -1 * ca_cb_orients_[ jj ];
			}
		}

	}

// for debug
// for( Size ii=1; ii<=order_strands_.size(); ++ii ) {
//	Size const st( order_strands_[ ii ] );
//	Size const begin( strands[ st ]->begin() );
//	Size const end( strands[ st ]->end() );
//	for( Size jj=begin; jj<=end; jj++ ) {
//		std::cout << jj << " " << ca_cb_orients_[ jj ] << std::endl;
//	}
//}

	is_geometry_initialized_ = true;
}

/// @brief
int
Sheet::which_side( Vector const vec ) const
{
	Real dot = sheet_plane_.dot( vec - sheet_center_ );
	if ( dot > 0 ) {
		return 1;
	} else {
		return -1;
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////
/// @brief default constructor
SheetSet::SheetSet()
{}

/// @brief value constructor
SheetSet::SheetSet( Sheets const & sheets ):
	sheets_( sheets )
{}

/// @brief value constructor
SheetSet::SheetSet( SS_Info2_COP const ssinfo, StrandPairingSetCOP const spairset )
{
	initialize( ssinfo, spairset );
}

/// @brief copy constructor
SheetSet::SheetSet( SheetSet const & s ) :
	ReferenceCount(),
	sheets_( s.sheets_ )
{}

/// @brief return strand pairing
std::ostream & operator<<( std::ostream & out, const SheetSet &s )
{
	out << "#### Sheet Info " << std::endl;
	for( Size i=1; i<=s.num_sheets() ; i++ ){
		SheetCOP sop ( s.sheet( i ) );
		if( sop->is_barrel() ){
			out << "## Barrel " << i << ":" << std::endl;
		}else{
			out << "## Sheet " << i << ":" << std::endl;
		}
		out << (*sop);
	}
	return out;
}

/// @brief
void
SheetSet::push_back( SheetOP const sop )
{
	sheets_.push_back( sop );
}

/// @brief
void SheetSet::clear()
{
	sheets_.clear();
}

/// @brief
SheetOP
SheetSet::sheet( Size const s ) const
{
	runtime_assert( s <= sheets_.size() );
	return sheets_[ s ];
}

/// @brief return a sheet number given a strand number
Size
SheetSet::which_sheet( Size const s ) const
{

	if( s == 0 ) return 0;
	return sheet_number_[ s ];
}

///	 @brief return strand pairing set
SheetSet::StrandPairingSet
SheetSet::spairset() const
{
	return spairset_;
}

/// @brief
void
SheetSet::set_sheet_number() const
{
	for( Size i=1; i<=sheets_.size(); i++ ){
		SheetCOP const sheet ( sheets_[ i ] );
		VecSize const & order ( sheet->order_strands() );
		for( VecSize::const_iterator it=order.begin(), ite=order.end(); it!=ite; ++it ){
			sheet_number_[ *it ] = i;
		}
	}
}

/// @brief
void
SheetSet::calc_geometry( SS_Info2_COP const ssinfo )
{
	for( Size i=1; i<=sheets_.size(); i++ ){
		sheets_[ i ]->calc_geometry( ssinfo );
	}
}

/// @brief
void
SheetSet::initialize( SS_Info2_COP const ssinfo, StrandPairingSetCOP const spairset_in )
{
	using core::graph::DisjointSets;
	using protocols::fldsgn::topology::StrandPairings;
	using protocols::fldsgn::topology::StrandPairingCOP;

	// initialize sheet_number_
	for( Size i=1; i<=ssinfo->strands().size(); i++ ) {
		sheet_number_.insert( std::map< Size, Size >::value_type( i, 0 ) );
	}

	// do nothing if strand pairs < 2
	if( spairset_in->num_strands() < 2 ) {
		return;
	}

	spairset_ = *( spairset_in->clone() );
	spairset_.make_strand_neighbor_two();

	runtime_assert( spairset_.finalized() );

	// strand pairings
	StrandPairings spairs = spairset_.strand_pairings();

	// number of strands included in spairset
	// Size num_strands = spairset.num_strands();

	// set neighbor strands and sheet_set
	DisjointSets sheet_set( ssinfo->strands().size() );
	for( StrandPairings::const_iterator iter = spairs.begin(); iter != spairs.end(); ++iter ) {
		StrandPairing const & sp( **iter );
		sheet_set.ds_union( sp.s1(), sp.s2() );
	}

	// calc order of strands
	std::map< Size, VecSize > sset =  sheet_set.sets();
	std::map< Size, VecSize >::iterator it = sset.begin(), end = sset.end();
	while( it != end ) {

		bool ibarrel (true);
		VecSize order_strands;

		VecSize list_strands = (*it).second;
		if( list_strands.size() <= 1 ) {
			++it;
			continue;
		}

		for( VecSize::const_iterator itt = list_strands.begin(); itt != list_strands.end(); ++itt ) {
			Size s( *itt );
			if( spairset_.neighbor_strands( s ).size() == 1 ) {
				order_strands.push_back( s );
				ibarrel = false;
			} else if( spairset_.neighbor_strands( s ).size() > 2 ) {
				TR.Error << "Error, num of neighbor strands > 2 " << std::endl;
				runtime_assert( false );
			}
		}

		if( ibarrel ){
			order_strands.push_back( list_strands.front() );
		}else{
			runtime_assert( order_strands.size() == 2 );
			sort ( order_strands.begin(), order_strands.end() );
			order_strands.erase( order_strands.end() - 1 );
		}

		VecSize const & neighbor ( spairset_.neighbor_strands( order_strands.front() ) );
		if( ibarrel ) {
			if( neighbor[1] < neighbor[2] ) {
				order_strands.insert( order_strands.begin() + 1, neighbor[1] );
			} else {
				order_strands.insert( order_strands.begin() + 1, neighbor[2] );
			}
		} else {
			order_strands.insert( order_strands.begin() + 1, neighbor[1] );
		}


		for( Size i=2; i<= list_strands.size()-1; i++ ){
			Size s1 = order_strands[ i ];
			VecSize const & neighbor( spairset_.neighbor_strands( s1 ) );

			for( Size j=1; j<=neighbor.size(); j++ ){
				if( neighbor[ j ] != order_strands[ i-1 ] ){
					order_strands.insert( order_strands.begin() + i, neighbor[j] );
				}
			}
		}

		// set orient of strands in sheet
		VecInt orient_strands; // Sergey: was declared as 'VecReal' seems to be a typo, converting to VecInt to fix wanrning
		orient_strands.push_back( 1 );
		for( Size i=1; i<=order_strands.size()-1; i++ ){
			Size s1( order_strands[ i ] );
			Size s2( order_strands[ i+1 ] );
			StrandPairingCOP const spairop( spairset_.strand_pairing( s1, s2 ) );

			if( spairop->orient() == 'P' ){
				orient_strands.push_back( orient_strands[ i ] );
			}else{
				orient_strands.push_back( orient_strands[ i ]*-1 );
			}
		}

		SheetOP sop( new Sheet( order_strands, orient_strands, ibarrel ) );
		sheets_.push_back( sop );

    ++it;

  } // while( it !=sset.end() )

	set_sheet_number();

} // SheetTopology::calc_sheetset


} // namespace topology
} // namespace fldsgn
} // namespace protocols

