// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// Copyright in the Rosetta software belongs to the developers and their institutions.
// For more information, see www.rosettacommons.org.

/// @file ./src/protocols/topology/BetaAlphaBetaMotif.cc
/// @brief  class for beta alpha beta motifs
/// @author Nobuyasu Koga ( nobuyasu@u.washington.edu )

// unit headers
#include <protocols/fldsgn/topology/BetaAlphaBetaMotif.hh>
#include <protocols/fldsgn/topology/StrandPairing.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/fldsgn/topology/HSSTriplet.hh>

// project headers
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <protocols/fldsgn/topology/Sheet.hh>

// utility headers
#include <utility/exit.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <cassert>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>

// AUTO-REMOVED #include <ObjexxFCL/format.hh>

#include <utility/vector1.hh>


static thread_local basic::Tracer TR( "protocols.fldsgn.topology.BetaAlphaBetaMotif" );

using namespace core;

namespace protocols {
namespace fldsgn {
namespace topology {

/// @brief default constructor
BetaAlphaBetaMotif::BetaAlphaBetaMotif():
	strand1_( 0 ),
	strand2_( 0 ),
	helix_( 0 ),
	cross_over_( 0 ),
	left_handed_( false ),
	hs_dist_( 0.0 ),
	hs_angle_( 0.0 ),
	hs1_dist_( 0.0 ),
	hs2_dist_( 0.0 ),
	hsheet_elev_angle_( 0.0 ),
	geometry_is_initialized_( false )
{
	helix_cycle_.push_back( 0 );
}

/// @brief value constructor
BetaAlphaBetaMotif::BetaAlphaBetaMotif(
  Size const & strand1,
	Size const & strand2,
	Size const & helix,
	Size const & cross_over ) :
	strand1_( strand1 ),
	strand2_( strand2 ),
	helix_( helix ),
	cross_over_( cross_over ),
	left_handed_( false ),
	hs_dist_( 0.0 ),
	hs_angle_( 0.0 ),
 	hs1_dist_( 0.0 ),
	hs2_dist_( 0.0 ),
	hsheet_elev_angle_( 0.0 ),
	geometry_is_initialized_( false )
{
	helix_cycle_.push_back( 0 );
}

/// @brief copy constructor
BetaAlphaBetaMotif::BetaAlphaBetaMotif( BetaAlphaBetaMotif const & s ):
	ReferenceCount(),
	strand1_( s.strand1_ ),
	strand2_( s.strand2_ ),
	helix_( s.helix_ ),
	cross_over_( s.cross_over_ ),
	left_handed_( s.left_handed_ ),
	sheet_plane_( s.sheet_plane_ ),
	sheet_pos_( s.sheet_pos_ ),
	hs_dist_( s.hs_dist_ ),
	hs_angle_( s.hs_angle_ ),
	hs1_dist_( s.hs1_dist_ ),
	hs2_dist_( s.hs2_dist_ ),
	hsheet_elev_angle_( s.hsheet_elev_angle_ ),
	helix_cycle_( s.helix_cycle_ ),
	geometry_is_initialized_( s.geometry_is_initialized_ )
{}

/// @brief destructor
BetaAlphaBetaMotif::~BetaAlphaBetaMotif(){}

/// @brief IO operator
std::ostream & operator<<(std::ostream & out, const BetaAlphaBetaMotif & s ) {

	if( s.is_lefthanded() ) {
		out << "# " << s.helix() << "," << s.strand1() << "-" << s.strand2() << "." << s.cross_over() << "." << "Left_handed";
	} else {
		out << "# " << s.helix() << "," << s.strand1() << "-" << s.strand2() << "." << s.cross_over();
	}
	out << " " << s.hsheet_dist() << " " << s.hs_angle() << " " << s.hsheet_elev_angle() << " " << s.hs1_dist() << " " << s.hs2_dist() << " "
		<< s.helix_cycle_as_string() << std::endl;

	return out;
}

/// @brief return name
std::string
BetaAlphaBetaMotif::name() const
{
	std::ostringstream name;
	name << helix() << "," << strand1() << "-" << strand2();
	return name.str();
}

/// @brief return helix cycle as string
std::string
BetaAlphaBetaMotif::helix_cycle_as_string() const
{
	std::ostringstream name;
	if( helix_cycle_[ 1 ] == 0 ) {
		name << "0";
	} else {
		if( helix_cycle_.size() == 1 ) {
			name << helix_cycle_[ 1 ];
		} else {
			if( helix_cycle_[ 1 ] < helix_cycle_[ 2 ] ) {
				name << helix_cycle_[ 1 ] << helix_cycle_[ 2 ];
			} else {
				name << helix_cycle_[ 2 ] << helix_cycle_[ 1 ];
			}
		}
	}
	return name.str();
}

/// @brief whether the CB->CA vector of the C-term residue of 1st strand is pointing inward or outward
/// 1: inward, 2: outward
core::Size
BetaAlphaBetaMotif::calc_inout( SS_Info2_COP const ssinfo, Size const resi ) const
{
	using core::Vector;
	using protocols::fldsgn::topology::BB_Pos;

	Real neighbor_dist( 14.0 );
	Size burial( 4 );

	BB_Pos bb_pos( ssinfo->bb_pos() );

	Size start = ssinfo->strand( strand1_ )->end() + 1;
	Size end = ssinfo->strand( strand2_ )->begin() - 1;

	Vector cbca = bb_pos.CB( resi ) - bb_pos.CA( resi );

	Size neighbor( 0 );
	for( Size ii=start; ii<=end; ii++ ) {

		Vector bb = bb_pos.CA( resi ) - bb_pos.CA( ii );

		Real orient = cbca.x()*( bb_pos.CA( ii ).x() - bb_pos.CA( resi ).x() ) +
			cbca.y()*( bb_pos.CA( ii ).y() - bb_pos.CA( resi ).y() ) +
			cbca.z()*( bb_pos.CA( ii ).z() - bb_pos.CA( resi ).z() );

		if( bb.length() <= neighbor_dist && orient > 0 ) {
			neighbor ++;
		}

		if( neighbor > burial ) {
			return 1;
		}
	}
	return 2;
}

/// @brief
bool compare( std::map< Size, Real >::const_iterator a, std::map< Size, Real >::const_iterator b )
{
	return ( (*a).second < (*b).second );
}



/// @brief calc helix cycle against sheet. Helix cycle is classified as 0, 1, 2, 3, 4,
/// which denote the position on helix, where the residue pointing to sheet plane.
///	0 means helix cycle is not calculated, or it's impossible to determine the helix cycle.
void
BetaAlphaBetaMotif::calc_helix_cycle( SS_Info2_COP const ssinfo )
{
	using protocols::fldsgn::topology::Helix;
	using protocols::fldsgn::topology::Strand;

	if( ! geometry_is_initialized_ ) {
		TR << "Geometry have to be initialized before calculating helix cycle. ";
		runtime_assert( false );
	}

	// get helix
	Helix  const hx ( *ssinfo->helix( helix_ ) );
	Strand const s1 ( *ssinfo->strand( strand1_ ) );
	Strand const s2 ( *ssinfo->strand( strand2_ ) );
	BB_Pos bb_pos( ssinfo->bb_pos() );

	Size h_end( 8 );
	if( hx.length() < 8 ) {
		h_end = hx.length();
	}


	bool use_ca( false );
	for( Size ii=1; ii<=h_end; ++ii ) {
		Size pos( hx.begin() + ii - 1 );
		if( bb_pos.CB( pos ).is_zero() ) use_ca = true;
	}

	std::map< Size, Real > pos_dist;
	for( Size ii=1; ii<=h_end; ++ii ) {

		Size pos( hx.begin() + ii - 1 );

		Vector hpos;
		if( use_ca ) {
			hpos = bb_pos.CA( pos );
		} else {
			hpos = bb_pos.CB( pos );
		}

		// get orient vector from the center of helical wheel to ca
		Vector const helix_center = hx.Nend_pos() + hx.Nend_orient().dot( hpos - hx.Nend_pos() )*hx.Nend_orient();
		Vector const center2hpos = ( hpos - helix_center ).normalized();

		// distance between ca and sheet
		Real const dist_hpos_sheet = sheet_plane_.dot( hpos - sheet_pos_ );

		if( dist_hpos_sheet < 0 ) {
			TR << " Distance between CA and sheet is negative value. Conformation is something funny. ";
			TR << " Cannot determine helix cycle. " << std::endl;
			TR << "Pos " << pos << ", " << dist_hpos_sheet << std::endl;
			return;
		}

		// angle between center2ca and sheet_plane
		Real cos = cos_of( center2hpos, -sheet_plane_ );
		if( cos < 0 ) continue;

		// distance from ca along the vector of center2ca
		Real dist = dist_hpos_sheet/cos;

		//std::cout << ii << " " << dist << std::endl;
		pos_dist.insert( std::pair< Size, Real >( ii, dist ) );

	}

	utility::vector1< std::map< Size, Real >::iterator > vec_ite;
	std::map< Size, Real >::iterator it = pos_dist.begin();
	while( it != pos_dist.end() ) {
		vec_ite.push_back( it );
		it++;
	}

	helix_cycle_.clear();
	std::sort( vec_ite.begin(), vec_ite.end(), compare );

	if( hx.length() > 8 ) {
		helix_cycle_.push_back( ( **vec_ite.begin() ).first );
		helix_cycle_.push_back( ( **( vec_ite.begin() + 1 ) ).first );
	} else {
		helix_cycle_.push_back( ( **vec_ite.begin() ).first );
	}

}


/// @brief
void
BetaAlphaBetaMotif::calc_geometry( SS_Info2_COP const ssinfo, SheetSetCOP const sheet_set )
{
	using protocols::fldsgn::topology::Helix;
	using protocols::fldsgn::topology::Strand;
	using protocols::fldsgn::topology::SheetCOP;
	using protocols::fldsgn::topology::HSSTriplet;
	using protocols::fldsgn::topology::StrandPairingSet;
	using protocols::fldsgn::topology::StrandPairingOP;

	if( ! ssinfo->bbpos_is_set() ) return;
	runtime_assert( ssinfo->bbpos_is_set() );

	Helix  const hx ( *ssinfo->helix( helix_ ) );
	Strand const s1 ( *ssinfo->strand( strand1_ ) );
	Strand const s2 ( *ssinfo->strand( strand2_ ) );
	BB_Pos bb_pos( ssinfo->bb_pos() );

	if( ! s1.is_geometry_initialized() ) return;
	if( ! s2.is_geometry_initialized() ) return;
	if( ! hx.is_geometry_initialized() ) return;

	HSSTriplet hss3( helix_, strand1_, strand2_ );
	hss3.calc_geometry( ssinfo );

	// get mid point of helix
	Vector const hmid = hx.mid_pos();

	//hs_dist_  = hss3.hsheet_dist();
	//hs_angle_ = hss3.hs_angle();
	left_handed_ = hss3.left_handed();

	// check strand1_ and strand2_ are in same sheet
	Size isheet = sheet_set->which_sheet( strand1_ );
	runtime_assert( isheet == sheet_set->which_sheet( strand1_ ) );
    SheetOP sheet( sheet_set->sheet( isheet ) );

	// get strand order of strand1 and strand2 in sheet
	Size s1_order = sheet->strand_order( strand1_ );
	Size s2_order = sheet->strand_order( strand2_ );
	runtime_assert( s1_order > 0 && s2_order > 0 );
	int sign( 1 );
	if( s1_order > s2_order ) sign = -1;

	// determine 1st reference strand for calculating geometry,
	// reference strands are center strands of sheets
	int ref1 = sign * ( cross_over_/2 + 1 ) + s1_order + sign*-1;

	runtime_assert( ref1 > 0 );

	int s1_orient = sheet->orient_strand( s1_order );

	// calc v1, v2, v3 to define sheet plane
	// v1 and v2 are orient vectors to define sheet plane,
	// v3 is a positional vector, which locates on sheet plane.
	Vector v1, v2, v3;

	int ref2 = ref1 + sign;
	runtime_assert( ref2 > 0 );

	Size ref1_stid = sheet->order_strand( ref1 );
	Size ref2_stid = sheet->order_strand( ref2 );

	// Strand const & ref_s1 = *ssinfo->strand( ref1_stid ); // Unused variable causes warning.
	// Strand const & ref_s2 = *ssinfo->strand( ref2_stid ); // Unused variable causes warning.

	int refs1_orient = sheet->orient_strand( ref1 );
	int refs2_orient = sheet->orient_strand( ref2 );

	Real ss_sign1( 1.0 ), ss_sign2( 1.0 );
	if( s1_orient != refs1_orient ) ss_sign1 = -1.0;
	if( s1_orient != refs2_orient ) ss_sign2 = -1.0;

	// get StrandPairingSet
	StrandPairingSet spairset = sheet_set->spairset();

	StrandPairingOP spair1 = spairset.strand_pairing( ref1_stid, ref2_stid );
	runtime_assert( " spair1->orient() == 'P' || spair1->orient() == 'A' " );
	runtime_assert( " spair1->size1() != 1 || spair1->size2() != 1 " );

	// get vectors to define sheet plane
	Size begin1, end1, begin2, end2;
	if( ref1_stid < ref2_stid ) {
		begin1 = spair1->begin1();
		end1   = spair1->end1();
		if( spair1->orient() == 'P' ) {
			begin2 = spair1->begin2();
			end2   = spair1->end2();
		} else {
			begin2 = spair1->end2();
			end2   = spair1->begin2();
		}

	} else {
		if( spair1->orient() == 'P' ) {
 			begin1 = spair1->begin2();
			end1   = spair1->end2();
		} else {
			begin1 = spair1->end2();
			end1   = spair1->begin2();
		}
		begin2 = spair1->begin1();
		end2   = spair1->end1();
	}

	Vector const s1p = ( bb_pos.C( end1 ) - bb_pos.N( begin1 ) ).normalized();
	Vector const s2p = ( bb_pos.C( end2 ) - bb_pos.N( begin2 ) ).normalized();
	Vector const s1p_mid = ( bb_pos.N( begin1 ) + bb_pos.C( end1 ) )/2.0;
	Vector const s2p_mid = ( bb_pos.N( begin2 ) + bb_pos.C( end2 ) )/2.0;

	if( cross_over_%2 == 0 ) {

		v1 = ( s2p_mid - s1p_mid ).normalized();
		v2 = ( ss_sign1*s1p + ss_sign2*s2p ).normalized();
		v3 = s1p_mid;

		hs1_dist_ = ( hmid - s1p_mid ).length();
		hs2_dist_ = ( hmid - s2p_mid ).length();

		//TR << begin1 << "-" << end1 << "," << begin2 << "-" << end2 << std::endl;
		//TR << "Ref " << ref1 << " " << ref2 << ", "
		//  << ref1_stid << " " << ref2_stid << ", "
		//   << refs1_orient << " " << refs2_orient << std::endl;

	} else {

		int ref3 = ref2 + sign;
		runtime_assert( ref3 > 0 );

		Size ref3_stid = sheet->order_strand( ref3 );
		// Strand const & ref_s3 = *ssinfo->strand( ref3_stid ); // Unused variable causes warning.
		int refs3_orient = sheet->orient_strand( ref3 );

		Real ss_sign3( 1.0 );
		if( s1_orient != refs3_orient ) ss_sign3 = -1.0;

		StrandPairingOP spair2 = spairset.strand_pairing( ref2_stid, ref3_stid );
		runtime_assert( " spair2->orient() == 'P' || spair2->orient() == 'A' " );
		runtime_assert( " spair2->size1() != 1 || spair2->size2() != 1 " );

		// get vectors to define sheet plane
		Size begin3, end3;
		if( ref2_stid < ref3_stid ) {
			if( spair2->orient() == 'P' ) {
				begin3 = spair2->begin2();
				end3   = spair2->end2();
			} else {
				begin3 = spair2->end2();
				end3   = spair2->begin2();
			}
		} else {
			begin3 = spair2->begin1();
			end3   = spair2->end1();
		}

		Vector const s3p = ( bb_pos.C( end3 ) - bb_pos.N( begin3 ) ).normalized();
		Vector const s3p_mid = ( bb_pos.N( begin3 ) + bb_pos.C( end3 ) )/2.0;

		hs1_dist_ = ( hmid - s1p_mid ).length();
		hs2_dist_ = ( hmid - s3p_mid ).length();

		v1 = ( s3p_mid - s1p_mid ).normalized();
		v2 = ( ss_sign1*s1p + ss_sign2*s2p + ss_sign3*s3p ).normalized();
		v3 = s1p_mid;

		//TR << begin1 << "-" << end1 << "," << begin2 << "-" << end2 << "," << begin3 << "-" << end3 << std::endl;
		//TR << "Ref " << ref1 << " " << ref2 << " " << ref3 << ", "
		//  << ref1_stid << " " << ref2_stid << " " << ref3_stid << ", "
		//   << s1_orient << " " << refs1_orient << " " << refs2_orient << " " << refs3_orient << std::endl;

	}

	sheet_plane_ = ( v1.cross( v2 ) ).normalized();
	sheet_pos_ = v3;
	// distance between helix and sheet
	hs_dist_ = sheet_plane_.dot( hmid - v3 );

	// angle between helix and strands
    Real dot2( -hx.orient().dot( sheet_plane_ ) );
	Real angle( numeric::conversions::degrees( angle_of( -hx.orient(), v2 ) ));
	if( dot2 > 0.0 ) {
		hsheet_elev_angle_ = angle;
	} else {
		hsheet_elev_angle_ = -angle;
	}

	// calc angle between strands and helix projected on sheet
	Vector const h1 = hx.Nend_pos() - sheet_plane_.dot( hx.Nend_pos() - v3 )*sheet_plane_;
	Vector const h2 = hx.Cend_pos() - sheet_plane_.dot( hx.Cend_pos() - v3 )*sheet_plane_;
	Vector const hx_on_sheet = ( h1 - h2 ).normalized();
	Real ori = cross( hx_on_sheet, v2 ).dot( sheet_plane_ );
	hs_angle_ = numeric::conversions::degrees( angle_of( v2, hx_on_sheet ) );
	if( ori < 0 ) hs_angle_ = -1*hs_angle_;

	geometry_is_initialized_ = true;

	calc_helix_cycle( ssinfo );

	// std::cout << *this;

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief default constructor
BetaAlphaBetaMotifSet::BetaAlphaBetaMotifSet() {}

/// @brief value constructor
BetaAlphaBetaMotifSet::BetaAlphaBetaMotifSet( BetaAlphaBetaMotifs const & bab_motifs ):
	bab_motifs_( bab_motifs )
{}


/// @brief value constructor
BetaAlphaBetaMotifSet::BetaAlphaBetaMotifSet( SS_Info2_COP const ssinfo, SheetSetCOP const sheet_set )
{
	set_babmotifs( ssinfo, sheet_set );
	calc_geometry( ssinfo, sheet_set );
}


/// @brief copy constructor
BetaAlphaBetaMotifSet::BetaAlphaBetaMotifSet( BetaAlphaBetaMotifSet const & s ):
	ReferenceCount(),
	bab_motifs_( s.bab_motifs_ )
{}

/// @brief destructor
BetaAlphaBetaMotifSet::~BetaAlphaBetaMotifSet(){}

/// @brief add BetaAlphaBetaMotif
void
BetaAlphaBetaMotifSet::push_back( BetaAlphaBetaMotifOP const bop )
{
	bab_motifs_.push_back( bop );
}

/// @brief
void
BetaAlphaBetaMotifSet::clear(){
	bab_motifs_.clear();
}

/// @brief
BetaAlphaBetaMotifs const &
BetaAlphaBetaMotifSet::bab_motifs() const
{
	return bab_motifs_;
}

/// @brief
BetaAlphaBetaMotifOP
BetaAlphaBetaMotifSet::bab_motif( Size const & i ) const
{
		runtime_assert( i <= bab_motifs_.size() );
		return bab_motifs_[ i ];
}

/// @brief
std::ostream & operator<<( std::ostream & out, const BetaAlphaBetaMotifSet & s )
{
	out << "### BetaAlphaBetaMotif Info " << std::endl;
	out << "# helix,strand1-strand2.num_crossover hsheet_dist hs_angl hsheet_elev_angl hs1_dist hs2_dist helix_cycle " << std::endl;
	for( BetaAlphaBetaMotifs::const_iterator iter = s.bab_motifs().begin(),
				 iter_end = s.bab_motifs().end(); iter != iter_end; ++iter ) {
		BetaAlphaBetaMotif const & bab( **iter );
		out << bab;
 	}
	return out;
}

/// @brief set bab motif
void
BetaAlphaBetaMotifSet::set_babmotifs( SS_Info2_COP const ssinfo, SheetSetCOP const sheet_set )
{
	using protocols::fldsgn::topology::Helix;
	using protocols::fldsgn::topology::Strand;
	using protocols::fldsgn::topology::Helices;
	using protocols::fldsgn::topology::Strands;

	Helices helices = ssinfo->helices();
	Strands strands = ssinfo->strands();

	if( strands.size() < 2 ) {
		return;
	}

	for( Size ist1=1; ist1<=strands.size()-1 ; ist1++ )  {

		Size ist2 = ist1+1;
		Strand const & s1 ( *strands[ ist1 ] );
		Strand const & s2 ( *strands[ ist2 ] );

		Size const sheet_num1 ( sheet_set->which_sheet( ist1 ) );
		Size const sheet_num2 ( sheet_set->which_sheet( ist2 ) );
		if( sheet_num1 == 0 || sheet_num2 == 0 ) continue;  // in the case, ist1 or ist2 does not belong to sheet
		if( sheet_num1 != sheet_num2 ) continue;

		Size const ord1 = sheet_set->sheet( sheet_num1 )->strand_order( ist1 );
		Size const ord2 = sheet_set->sheet( sheet_num2 )->strand_order( ist2 );

		utility::vector1< Real > orients( sheet_set->sheet( sheet_num1 )->orient_strands() );
		if( orients[ ord1 ] != orients[ ord2 ] ) continue;

		Size ih( 0 ), ihelix( 0 );
		Size helix_num( 0 );
		for( Helices::const_iterator jt=helices.begin(), jte=helices.end(); jt!=jte; ++jt ) {
			ih++;
			Helix const & hx( **jt );
			if( hx.begin() < s1.begin() ) continue;
			if( hx.begin() > s2.end() ) continue;
			if( hx.begin() > s1.end() && hx.end() < s2.begin() ) {
				ihelix = ih;
				helix_num++;
			}
		}//helices

		if( helix_num == 1 ) {
			runtime_assert( ihelix != 0 );
			Size cross_over = Size ( std::abs( Real(ord1) - Real(ord2) ) ) - 1;
			bab_motifs_.push_back( BetaAlphaBetaMotifOP( new BetaAlphaBetaMotif( ist1, ist2, ihelix, cross_over ) ) );
			TR.Debug << ist1 << " " << ist2 << " " << ihelix << " " << cross_over << std::endl;
		}

	}//strands


} // SheetTopology::set_bab_motifs


/// @brief
void
BetaAlphaBetaMotifSet::calc_geometry( SS_Info2_COP const ssinfo,  SheetSetCOP const sheet_set )
{
	if( ! ssinfo->bbpos_is_set() ) return;
	runtime_assert( ssinfo->bbpos_is_set() );
	for( BetaAlphaBetaMotifs::iterator it=bab_motifs_.begin(), ite=bab_motifs_.end(); it!=ite; ++it ) {
		BetaAlphaBetaMotifOP const babm( *it );
		babm->calc_geometry( ssinfo, sheet_set );
	}
}


} // namespace topology
} // namespace fldsgn
} // namespace protocols

