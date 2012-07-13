// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ./src/protocols/fldsgn/topology/SS_Info2.cc
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@u.washington.edu )

/// Unit headers
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/fldsgn/topology/util.hh>

/// Project headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <core/pose/Pose.hh>
#include <core/util/ABEGOManager.hh>
#include <basic/datacache/CacheableData.hh>

/// C++ Headers
#include <iostream>
#include <cassert>

// Numeric
#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>

#include <utility/vector1.hh>
#include <boost/lexical_cast.hpp>


static basic::Tracer TR( "protocols.fldsgn.topology.SS_Info2" );

using namespace core;

namespace protocols {
namespace fldsgn {
namespace topology {

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SS_Base
/// @brief default constructor
SS_Base::SS_Base():
	name_( 'D' ),
	begin_( 0 ),end_( 0 ),
	is_geometry_initialized_( false )
{
	Vector v;
	v.zero();
	orient_ = v;
	Nend_orient_ = v;
	Cend_orient_ = v;
	Nend_pos_ = v;
	Cend_pos_ = v;
	mid_pos_ = v;
}

/// @brief default constructor
SS_Base::SS_Base( char const & name ):
	name_( name ),
	begin_( 0 ),end_( 0 ),
	is_geometry_initialized_( false )
{
	Vector v;
	v.zero();
	orient_ = v;
	Nend_orient_ = v;
	Cend_orient_ = v;
	Nend_pos_ = v;
	Cend_pos_ = v;
	mid_pos_ = v;
}

/// @brief value constructor
SS_Base::SS_Base( char const & name, Size const & begin, Size const & end ):
	name_( name ),
	begin_( begin ), end_( end ),
	is_geometry_initialized_( false )
{
	Vector v;
	v.zero();
	orient_ = v;
	Nend_orient_ = v;
	Cend_orient_ = v;
	Nend_pos_ = v;
	Cend_pos_ = v;
	mid_pos_ = v;
}


/// @brief copy constructor
SS_Base::SS_Base(	SS_Base const & s ):
	ReferenceCount(),
	name_( s.name_ ),
	begin_ ( s.begin_ ),
	end_ ( s.end_ ),
	is_geometry_initialized_( s.is_geometry_initialized_ ),
	orient_ ( s.orient_ ),
	Nend_orient_( s.Nend_orient_ ),
	Cend_orient_( s.Cend_orient_ ),
	Nend_pos_( s.Nend_pos_ ),
	Cend_pos_( s.Cend_pos_ ),
	mid_pos_( s.mid_pos_ )
{}


/// @brief destructor
SS_Base::~SS_Base(){}

/// @brief check the given residue is a member of this SS_Base
bool
SS_Base::is_member( Size const s ) const
{
	if( s >= begin_ && s <= end_ ) {
		return true;		
	} else {
		return false;
	}
}
	

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief default constructor
Strand::Strand():
	SS_Base( 'E' )
{}


/// @brief value constructor
Strand::Strand( Size const & begin, Size const & end ):
	SS_Base( 'E', begin, end )
{}


/// @brief copy constructor
Strand::Strand(	Strand const & s ):
	SS_Base( s ),
	bulges_( s.bulges_ )
{}

	
/// @brief destructor
Strand::~Strand(){}


/// @brief
std::ostream & operator<<(std::ostream & out, const Strand & st )
{
	out << st.begin() << "-" << st.end() ;
	return out;
}

// @brief strand has bulge ?
bool
Strand::has_bulge() const
{
	if( bulges_.size() > 0 ) {
		return true;
	} else {
		return false;		
	}
}

// @brief
bool
Strand::set_bulge( Size const & s )
{
	runtime_assert( is_member( s ) );
	
	for( Size i=1; i<= bulges_.size(); i++ ) {
		if( bulges_[ i ] == s ) {
			return false;
		}
	}
	
	bulges_.push_back( s );
	return true;	
}
	
	
/// @brief set vector between Calpha atoms of edge residues
void
Strand::calc_geometry( BB_Pos const & bbpos )
{
	is_geometry_initialized( true );

	// orient is defined by begin and end Ca atoms of strand
	Nend_pos( (bbpos.N( begin() ) + bbpos.C( begin() ))/2.0 );
	Cend_pos( (bbpos.N( end() ) + bbpos.C( end() ))/2.0 );
	Vector v = Cend_pos() - Nend_pos();
	orient( v.normalized() );

	Size pos = ( end() - begin() )/2 + begin();
	if( ( end() - begin() )%2 == 0 ) {
		mid_pos( ( bbpos.N( pos ) + bbpos.C( pos ) )/2.0 );
	} else {
		mid_pos( ( bbpos.N( pos ) + bbpos.C( pos+1 ) )/2.0 );
	}

	Vector n = ( bbpos.C( begin() ) - bbpos.N( begin() ) )/2.0;
	Vector c = ( bbpos.C( end()   ) - bbpos.N( end()   ) )/2.0;
	Nend_orient( n.normalized() );
	Cend_orient( c.normalized() );
	
//	Nend_orient( ( mid_pos() - Nend_pos() ).normalized() );
//	Cend_orient( ( Cend_pos() - mid_pos() ).normalized() );
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief default constructor
Helix::Helix():
	SS_Base( 'H' ),
	bend_( 0 )
{}


/// @brief value constructor
Helix::Helix( Size const & begin, Size const & end ):
	SS_Base( 'H', begin, end ),
	bend_( 0 )
{}


/// @brief copy constructor
Helix::Helix( Helix const & s ):
	SS_Base( s ),
	bend_( s.bend_ )
{}


/// @brief destructor
Helix::~Helix(){}


/// @brief
std::ostream & operator<<(std::ostream & out, const Helix & hx )
{
	out << hx.begin() << "-" << hx.end() ;
	return out;
}

/// @brief
void
Helix::calc_geometry( BB_Pos const & bbpos )
{
	is_geometry_initialized( true );

	Size begin( (*this).begin() );
	Size end( (*this).end() );

	static Real const eleven_inv = 1.0 / 11.0;

	Size const s1 = int ( begin );
	Size const s2 = s1 + 1;
	Size const s3 = s2 + 1;
	Size const s4 = s3 + 1;

	Vector const p1 = (                 bbpos.CA( s1 ) + bbpos.C( s1 ) ) +
		                ( bbpos.N( s2 ) + bbpos.CA( s2 ) + bbpos.C( s2 ) ) +
									 	( bbpos.N( s3 ) + bbpos.CA( s3 ) + bbpos.C( s3 ) ) +
									 	( bbpos.N( s4 ) + bbpos.CA( s4 )                 ) ;

	Vector const Nend_pos1 = ( p1 + bbpos.N( s1 ) ) * eleven_inv;
	Vector const Nend_pos2 = ( p1 + bbpos.C( s4 ) ) * eleven_inv;
	Nend_orient( ( Nend_pos2 - Nend_pos1 ).normalized() );

	Size const n1 = int ( end - 3 );
	Size const n2 = n1 + 1;
	Size const n3 = n2 + 1;
	Size const n4 = n3 + 1;

	Vector const p2 = (                 bbpos.CA( n1 ) + bbpos.C( n1 ) ) +
										( bbpos.N( n2 ) + bbpos.CA( n2 ) + bbpos.C( n2 ) ) +
		                ( bbpos.N( n3 ) + bbpos.CA( n3 ) + bbpos.C( n3 ) ) +
									 	( bbpos.N( n4 ) + bbpos.CA( n4 )                 ) ;

	Vector const Cend_pos1 = ( p2 + bbpos.N( n1 ) ) * eleven_inv;
	Vector const Cend_pos2 = ( p2 + bbpos.C( n4 ) ) * eleven_inv;
	Cend_orient( ( Cend_pos2 - Cend_pos1 ).normalized() );

	// set vector of helix
	Vector v = Cend_pos2 - Nend_pos1;
	orient( v.normalized() );
	// set positonal vector of C- and N-terminal
	Nend_pos( Nend_pos1 );
	Cend_pos( Cend_pos2 );
	// set mid point
	mid_pos( ( Nend_pos1 + Cend_pos2 )/2.0 );
	// set bend of the helix
	bend_ = numeric::conversions::degrees( angle_of( Nend_orient(), Cend_orient() ) );

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief default constructor
Loop::Loop():
	SS_Base( 'L' ),
	type_("")
{}


/// @brief value constructor
Loop::Loop( Size const & begin, Size const & end, String const & type ):
	SS_Base( 'L', begin, end ),
	type_( type )
{}


/// @brief copy constructor
Loop::Loop( Loop const & s ):
	SS_Base( s ),
	type_( s.type_ )
{}


/// @brief destructor
Loop::~Loop(){}


/// @brief
std::ostream & operator<<(std::ostream & out, const Loop & lp )
{
	out << lp.begin() << "-" << lp.end() << lp.type();
	return out;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief default constructor
SS_Info2::SS_Info2():
	CacheableData(),
	bbpos_is_set_( false ),
	secstruct_( "" ),
	min_helix_length_( 0 ),
	min_strand_length_( 0 ),
	initialized_by_spairset_( false )
{}


/// @brief value constructor
SS_Info2::SS_Info2(	String const & secstruct ):
	CacheableData(),
	bbpos_is_set_( false ),
	min_helix_length_( 0 ),
	min_strand_length_( 0 ),
	initialized_by_spairset_( false )
{
	initialize( secstruct );
}


/// @brief value constructor
SS_Info2::SS_Info2(	Pose const & pose, String const & secstruct ):
	CacheableData(),
	min_helix_length_( 0 ),
	min_strand_length_( 0 ),
	initialized_by_spairset_( false )
{
	initialize( pose, secstruct );
}


/// @brief copy constructor
SS_Info2::SS_Info2( SS_Info2 const & s ):
	CacheableData(),
	bbpos_is_set_( s.bbpos_is_set_ ),
	min_helix_length_( s.min_helix_length_ ),
	min_strand_length_( s.min_strand_length_ ),
	secstruct_( s.secstruct_ ),
	abego_( s.abego_ ),
	bb_pos_( s.bb_pos_ ),
	strands_( s.strands_ ),
	strand_id_( s.strand_id_ ),
	helices_( s.helices_ ),
	helix_id_( s.helix_id_ ),
	loop_id_( s.loop_id_ ),
	ss_element_id_( s.ss_element_id_ ),
	ss_element_id2op_( s.ss_element_id2op_ ),
	HE_elements_( s.HE_elements_ ),
	hleid_to_ssid_( s.hleid_to_ssid_ ),
	ssid_to_hleid_( s.ssid_to_hleid_ ),
	chiral_sstriplet_( s.chiral_sstriplet_ ),
	initialized_by_spairset_( s.initialized_by_spairset_ )
{}


/// @brief destructor
SS_Info2::~SS_Info2(){}


/// @brief make clone
basic::datacache::CacheableDataOP
SS_Info2::clone() const
{
	return new SS_Info2( *this );
}


/// @brief
void
SS_Info2::clear_data()
{
	bbpos_is_set_ = false;
	secstruct_ = "";
	abego_.clear();
	bb_pos_.clear();
	strands_.clear();
	helices_.clear();
	strand_id_.clear();
	helix_id_.clear();
	loop_id_.clear();
	ss_element_id_.clear();
	ss_element_id2op_.clear();
	HE_elements_.clear(),
	hleid_to_ssid_.clear();
	ssid_to_hleid_.clear();
	chiral_sstriplet_.clear();
	initialized_by_spairset_ = false;
}
	
/// @brief resize vectors
void
SS_Info2::resize( Size const nres )
{
	bb_pos_.resize( int(nres) );
}


/// @brief set secondary structure
void
SS_Info2::secstruct( Size const ii, char const & sec  )
{	
	if ( sec != 'E' && sec != 'H' && sec != 'L' && sec != 'D' ) {
		TR.Error << "unrecognized secstruct char : " << sec << std::endl;
		utility_exit();
	}
	secstruct_.at( ii-1 ) = sec;
}
	
/// @brief set abego
void
SS_Info2::abego( Size const ii, String const & aa )
{
	abego_.at( ii-1 ) = aa;
}
	
/// @brief initialize parameters of this class
void
SS_Info2::initialize( Pose const & pose, String const & secstruct )
{
	// data all clear
	clear_data();
	
	bbpos_is_set_ = true;
	if( secstruct == "" ){
		secstruct_ = pose.secstruct();
	}else{
		secstruct_ = secstruct;
		runtime_assert( pose.total_residue() == secstruct_.length() );
	}
	resize( pose.total_residue() );
	
	// remove helix and strand from secstruct depending on their minmum lengths
	if( min_helix_length_ > 0 || min_strand_length_ > 0 ) {
		reduce_secstruct();
	}
	
	// identify secondary structures
	identify_ss( secstruct_ );
	
	// set bb_pos
	bb_pos_.take_coordinates_from_pose( pose );
	
	// set abego info
	abego_ = core::util::get_abego( pose, 1 /* abego level */ );
	
	// set orient of strands and helices
	set_SSorient();
	
	// This will be removed !!! 
	// set chirality for consecutive 3 secondary structure elements
	set_chiral_consecutive_sstriplets();
}


/// @brief initialize parameters of this class
void
SS_Info2::initialize( String const & secstruct )
{
	clear_data(); // data all clear
	Size nres( secstruct.size() );
	runtime_assert( nres > 0 );
	secstruct_ = secstruct;
	resize( nres );
	// remove helix and strand from secstruct depending on their lengths
	if( min_helix_length_ > 0 || min_strand_length_ > 0 ) {
		reduce_secstruct();
	}
	identify_ss( secstruct_ );
}

	
/// @brief set bulges using StrandPairingSet
void
SS_Info2::set_bulge( StrandPairingSetCOP const spairset )
{

	runtime_assert( spairset->finalized() );
	runtime_assert( strands().size() >= spairset->num_strands() );
  runtime_assert( bbpos_is_set_ );

	// get total residue  NEEDS TO BE MODIFIED FOR LIGAND
	Size nres( secstruct().length() );
	
	utility::vector1< bool > flag( secstruct_.length(), false );
	for( Size ii=1; ii<=spairset->strand_pairings().size(); ii++ ) {

		StrandPairingOP spop( spairset->strand_pairing( ii ) );
		for( Size jj=1; jj<=spop->residue_pair_vec1().size(); jj++ ) {
			runtime_assert( spop->residue_pair_vec1()[ jj ] <= nres );
			flag[	spop->residue_pair_vec1()[ jj ] ] = true;
		}
		
		for( Size jj=1; jj<=spop->residue_pair_vec2().size(); jj++ ) {
			runtime_assert( spop->residue_pair_vec2()[ jj ] <= nres );
			flag[	spop->residue_pair_vec2()[ jj ] ] = true;
		}
		
	}
	
	for( Size ii=1; ii<=spairset->strand_pairings().size(); ii++ ) {

		StrandPairingOP spop( spairset->strand_pairing( ii ) );
		
		StrandOP s1op = strand( spop->s1() );
		StrandOP s2op = strand( spop->s2() );		
		
		for( Size ll=1; ll<=2; ++ll ) {

			utility::vector1<Size> bulges;
			if( ll==1 ) {
				bulges = spop->bulges1();
			} else {
				bulges = spop->bulges2();
			}
						
			for( Size jj=1; jj<=bulges.size(); jj++ ) {
			
				Size res = bulges[ jj ];
				if( flag[ res ] ) continue;

				Size before = res - 1;
				Size after = res + 1;
				runtime_assert( before > 1 && after < nres );
				runtime_assert( ss_element_id( before ) == ss_element_id( res ) && ss_element_id( after ) == ss_element_id( res ) );

				Vector cbca   = ( bb_pos().CB( res )   - bb_pos().CA( res ) ).normalized();
				Vector cbca_b = ( bb_pos().CB( before ) - bb_pos().CA( before ) ).normalized();
				Vector cbca_a = ( bb_pos().CB( after ) - bb_pos().CA( after ) ).normalized();
			
				//if( cbca.dot( cbca_a ) < 0.0 && cbca.dot( cbca_b ) < 0.0 ) continue;
				Real angle1 = numeric::conversions::degrees( angle_of( cbca, cbca_b ) );
				Real angle2 = numeric::conversions::degrees( angle_of( cbca, cbca_a ) );
	
				if( angle1 > 120.0 && angle2 > 120.0 ) continue;
				
				bool bulge_on( false );
				if( ll==1 ) {
					bulge_on = s1op->set_bulge( res );
				} else {
					bulge_on = s2op->set_bulge( res );
				}
				if ( basic::options::option[ basic::options::OptionKeys::run::debug ]() ) {
					TR << spop->s1() << "-" << spop->s2() << "," << res << " " << angle1 << " " << angle2 << std::endl;				
				}
			}
		}
		
	}
	initialized_by_spairset_ = true;
}

/// @brief set bulges
utility::vector1<Size> const
SS_Info2::bulges( Size const id ) const
{
	runtime_assert( initialized_by_spairset_ );
	return strand( id )->bulges();
}
	
/// @brief make string of HE elements
std::string
SS_Info2::make_string_HE_elements() const
{
	String ss;
	std::ostringstream ssindex;
	for( Size ii=1; ii<=HE_elements_.size(); ++ii ) {	
		char s = ss_element_id2op_[ HE_elements_[ ii ] ]->name();
		ssindex << s;		
	}
	return ssindex.str();
}
	
/// @brief make string of HE elements
std::string
SS_Info2::make_string_chiralities() const
{
	String ss;
	std::ostringstream ssindex;
	for( Size ii=1; ii<=HE_elements_.size()-2; ++ii ) {	
		char s = chiral_sstriplet_[ HE_elements_[ ii ] ];
		ssindex << s;		
	}
	return ssindex.str();
}	
	
/// @brief
std::ostream & operator<<(std::ostream & out, const SS_Info2 & ssinfo )
{
	Size count( 0 );
	out << "#### SS_Info " << std::endl;
	out << "# SS order: " << ssinfo.make_string_HE_elements() << std::endl;
	out << "# SS chirality: " << ssinfo.make_string_chiralities() << std::endl;
	for ( Helices::const_iterator it=ssinfo.helices_.begin(), ite=ssinfo.helices_.end(); it != ite; ++it ) {
		count ++;
		Helix const & hx( **it );
		out << "# Helix  " << count << ": " << hx << std::endl;
	}
	count = 0;
	for ( Strands::const_iterator it=ssinfo.strands_.begin(), ite=ssinfo.strands_.end(); it != ite; ++it ) {
		count ++;
		Strand const & st( **it );
		out << "# Strand " << count << ": " << st << std::endl;
	}
	return out;
}


/// @brief set orientation vector of secondary structures given a pose
void
SS_Info2::set_SSorient( Pose const & pose )
{
	runtime_assert( pose.total_residue() == bb_pos_.size() );
	bb_pos_.take_coordinates_from_pose( pose );
	bbpos_is_set_ = true;
	set_SSorient();
}


/// @brief set orientation vector of secondary structures
void
SS_Info2::set_SSorient()
{
	set_SSorient_strand();
	set_SSorient_helix();
}

/// @brief set orientation vector of strand
void
SS_Info2::set_SSorient_strand()
{
	for( Strands::iterator iter = strands_.begin(),
		iter_end = strands_.end(); iter != iter_end; ++iter ) {
		StrandOP & st( *iter );
		if(  st->length() >= 2 ){
			st->calc_geometry( bb_pos_ );
		}
	} // Strands
}

/// @brief set orientation vector of helix
void
SS_Info2::set_SSorient_helix()
{
	for( Helices::iterator iter = helices_.begin(),
		iter_end = helices_.end(); iter != iter_end; ++iter ) {
		HelixOP & hx( *iter );
		if(  hx->length() >= 4 ){
			hx->calc_geometry( bb_pos_ );
		}
	} // Helices
}

/// @brief set chirality of consecutive secondary structure elements
void
SS_Info2::set_chiral_consecutive_sstriplets()
{
	using protocols::fldsgn::topology::calc_chirality_sstriplet;

	Size max_ssid = ss_element_id_.back();
	chiral_sstriplet_.resize( max_ssid );
	for( Size ii=1; ii<=max_ssid; ++ii ) {
		chiral_sstriplet_[ ii ] = 'C';		
	}

	if( HE_elements_.size() < 3 ) return;
	
	for( Size ii=1; ii<=HE_elements_.size()-2 ; ++ii ) {
		Size id1 = HE_elements_[ ii ];
		Size id2 = HE_elements_[ ii+1 ];
		Size id3 = HE_elements_[ ii+2 ];		
		char chiral = calc_chirality_sstriplet( ss_element_id2op( id1 ), ss_element_id2op( id2 ), ss_element_id2op( id3 ) );
		chiral_sstriplet_[ id1 ] = chiral;
	}
}

/// @brief get SS_elemen id from helix, strand, or loop index
Size 
SS_Info2::get_ssid_by_hleid( char const & ss, Size const & index )
{
	if ( ss != 'L' && ss != 'H' && ss != 'E' ) {
		TR << "get_hleid_ssid only get 'L', 'E', and 'H' " << ss << std::endl;
		utility_exit();
	} 
	
	std::ostringstream ssindex;
	ssindex << ss << "-" << index;
	
	if( ( ss == 'L' &&   loops_.size() < index ) ||
			( ss == 'H' && helices_.size() < index ) ||		
 			( ss == 'E' && strands_.size() < index ) )
	{		
		TR << "Bigger index number was specified. index=" << index << std::endl;
		utility_exit();
	}

	return hleid_to_ssid_[ ssindex.str() ];

}

/// @brief get SS_elemen id from helix, strand, or loop index
std::string
SS_Info2::get_hleid_by_ssid( Size const & index )
{
	if( index > ss_element_id_.back() ) {
		TR << "Bigger index number was specified. index=" << index << std::endl;
		utility_exit();
	}		
	return ssid_to_hleid_[ index ];
}	
	
/// @brief get chirality for consecutive 3 secondary structur elements given a first ss element described H, E, or L w/ its index
char
SS_Info2::get_chiral( char const & ss, Size const & index )
{
	Size id = get_ssid_by_hleid( ss, index );
	return get_chiral( id );		
}
	
	
/// @brief get chirality for consecutive 3 secondary structur elements given a first ss element id
char
SS_Info2::get_chiral( Size const id ) const
{
	
	runtime_assert( bbpos_is_set_ );
	
	if( ss_element_id_.back() - 2 < id ) {
		TR << "Bigger index number was specified." << index;
		utility_exit();
	}
	
	char chiral = chiral_sstriplet_[ id ];		
	if( chiral == 'C' ) {
		TR << "Id of ss element of which chirality was not defined was specified: the position can be loop or last ss elements. " << id << std::endl;
		utility_exit();		
	}
	
	return chiral;
}
	
/// @brief relate helix, strand, or loop index to SS element id
void
SS_Info2::set_hleid_ssid( char const & ss, Size const & index, Size const & id )
{
	runtime_assert( ss == 'L' || ss == 'H' || ss == 'E' );
	std::ostringstream ssindex;
	ssindex << ss << "-" << index;
	hleid_to_ssid_.insert( std::map< String, Size >::value_type( ssindex.str(), id ) );		
	ssid_to_hleid_.insert( std::map< Size, String >::value_type( id, ssindex.str() ) );		
}

/// @brief modify secstruct based on the length of ss elements, depending on min_helix_length_ and min_strand_length_
void
SS_Info2::reduce_secstruct()
{
	Size nres( secstruct().length() );
	
	Size begin( 1 );
	String prev( "X" );
	utility::vector1< bool > modify( nres, false );

	for( Size ii=1; ii<=nres; ++ii ) {
		
		String const & ss( secstruct().substr( ii-1, 1 ) );

		if( ss != prev && ii != 1 ) {
			if(( prev == "E" && (ii - begin) < min_strand_length_ ) || 
				 ( prev == "H" && (ii - begin) < min_helix_length_ ) ) {
				for( Size jj=begin; jj<=ii-1; ++jj ) {
					modify[ jj ] = true;
				}
			}			
			begin = ii;
		}
	
		prev = ss;
	}
	
	if(( prev == "E" && (nres - begin) < min_strand_length_ ) || 
		 ( prev == "H" && (nres - begin) < min_helix_length_ ) ) {
		for( Size jj=begin; jj<=nres; ++jj ) {
			modify[ jj ] = true;					
		}
	}	
		
	for( Size ii=1; ii<=nres; ++ii ) {
		if( modify[ ii ] ) {
			secstruct( ii, 'L' );
		}
	}
	
}
	
	
	
/// ALERT!! Need to be changed for multiple chains and w/ ligand
/// @brief identify the region of each secondary structure element
void
SS_Info2::identify_ss( String const & secstruct )
{
	bool flag_L( false );
	bool flag_E( false );
	bool flag_H( false );
	Size beginE, beginH, beginL;
	Size istrand( 0 ), ihelix( 0 ), iloop( 0 ), iss( 0 );
	
	Size nres( secstruct.length() );
	helix_id_.resize( nres );
	strand_id_.resize( nres );
	loop_id_.resize( nres );
	ss_element_id_.resize( nres );

	String prev( "X" );
	Size jss( 0 );
	for( Size i=1; i<=nres; ++i ) {
		String const & ss( secstruct.substr( i-1, 1 ) );

		strand_id_[ i ] = 0;
		helix_id_ [ i ] = 0;
		loop_id_  [ i ] = 0;

		if( ss =="E" ) {

			if ( flag_E == false ) {
				istrand ++;
				beginE = i;
			}
			flag_E = true;
			strand_id_[ i ] = istrand;

		} else if( ss == "H" ) {

			if ( flag_H == false ) {
				ihelix++;
				beginH = i;
			}
			flag_H = true;
			helix_id_[ i ] = ihelix;

		} else {

			if( flag_L == false ) {
				iloop++;
				beginL = i;
			}
			flag_L = true;
			loop_id_[ i ] = iloop;

		}

		if ( ss !="E" && flag_E == true ) {
			flag_E = false;
			if( ( i - beginE ) >= 2 ) {
				strands_.push_back( new Strand( beginE, i-1 ) );
			} else {
				strands_.push_back( new Strand( beginE, beginE ) );
			}			
			iss++;
			ss_element_id2op_[ iss ] = strands_.back();
			set_hleid_ssid( 'E', strands_.size(), iss );
			HE_elements_.push_back( iss );
		}

		if ( ss !="H" && flag_H == true ) {
			flag_H = false;
			if( ( i-beginH ) >= 2 ) {
				helices_.push_back( new Helix( beginH, i-1 ) );
			} else {
				helices_.push_back( new Helix( beginH, beginH ) );
			}			
			iss++;
			ss_element_id2op_[ iss ] = helices_.back();
			set_hleid_ssid( 'H', helices_.size(), iss );
			HE_elements_.push_back( iss );
		}

		if( ss !="L" && flag_L == true ) {
			flag_L = false;
			if( ( i-beginL ) >= 2 ) {
				loops_.push_back( new Loop( beginL, i-1 ) );
			} else {
				loops_.push_back( new Loop( beginL, beginL ) );
			}			
			iss++;
			ss_element_id2op_[ iss ] = loops_.back();
			set_hleid_ssid( 'L',loops_.size(), iss );				
		}

		if( prev != ss ) {
			jss ++;
			prev = ss;
		}

		ss_element_id_[ i ] = jss;
		
		//std::cout << i << " " << ss_element_id_[ i ] << std::endl;
	} // for( Size i )

	if( flag_E ) {
		if( ( secstruct.length() - beginE + 1 ) >= 2 ) {
			strands_.push_back( new Strand( beginE, secstruct.length() ) );
		} else {
			strands_.push_back( new Strand( beginE, beginE ) );
		}
		iss++;
		ss_element_id2op_[ iss ] = strands_.back();
		set_hleid_ssid( 'E', strands_.size(), iss );
		HE_elements_.push_back( iss );
	}
	if( flag_H ) {
		if( ( secstruct.length() - beginH + 1 ) >= 2 ) {
			helices_.push_back( new Helix( beginH, secstruct.length() ) );
		} else {
			helices_.push_back( new Helix( beginH, beginH ) );
		}
		iss++;
		ss_element_id2op_[ iss ] = helices_.back();
		set_hleid_ssid( 'H', helices_.size(), iss );
		HE_elements_.push_back( iss );
	}
	if( flag_L ) {
		if( ( secstruct.length() - beginL + 1 ) >= 2 ) {
			loops_.push_back( new Loop( beginL, secstruct.length() ) );
		} else {
			loops_.push_back( new Loop( beginL, beginL ) );
		}
		iss++;
		ss_element_id2op_[ iss ] = loops_.back();
		set_hleid_ssid( 'L', loops_.size(), iss );
	}
	
	if ( basic::options::option[ basic::options::OptionKeys::run::debug ]() ) {
  	TR << make_string_HE_elements() << std::endl;
	 	for( Size ii=1; ii<=ss_element_id_.back(); ++ii ) {
			TR << "ssid: " << ii << " " << ss_element_id2op_[ ii ]->name() << " " << ss_element_id2op_[ ii ] << std::endl;			
		}
		for( Size ii=1; ii<=helices_.size(); ++ii ) {
			TR << "Helix: " << ii << " " << get_ssid_by_hleid( 'H', ii )<< std::endl;
		}
		for( Size ii=1; ii<=strands_.size(); ++ii ) {
			TR << "Strand: " << ii << " " << get_ssid_by_hleid( 'E', ii )<< std::endl;
		}
		for( Size ii=1; ii<=loops_.size(); ++ii ) {
			TR << "Loop: " << ii << " " << get_ssid_by_hleid( 'L', ii )<< std::endl;
		}
	}

} // SS_Info2::identify_ss()

/// @brief redfine the secondary structure with abego
/// Non-consistent positions between dssp and abego turned into loop	
/// The strands that consists of only one residue will be deleted.	
void
SS_Info2::redefine_with_abego()
{
	if( secstruct() == "" || !bbpos_is_set() ) {
		utility_exit_with_message( "No initialization with pose before calling redefine_with_abego" );
	}

	Size nres( secstruct().length() );

	// detect bulge
	for( Size ii=1; ii<=nres; ii++ ) {
		if( secstruct( ii ) == 'E' ) {
			if( abego( ii ) != "B" && abego( ii ) != "Z" && abego( ii ) != "Y" &&
			    abego( ii ) != "S" && abego( ii ) != "P" ) {
				secstruct( ii, 'L' );
			}
		}
		if( secstruct( ii ) == 'H' ) {
			if( abego( ii ) != "A" && abego( ii ) != "M" && abego( ii ) != "N" ) {
				secstruct( ii, 'L' );
			}
		}
	}

	bool readE( false ), readH( false );
	Size lenE( 0 ), lenH( 0 );
	for( Size i=1; i<=nres; ++i ) {
		String const & ss( secstruct().substr( i-1, 1 ) );
		
		if( readE && ss != "E" ) {
			if( lenE == 1 ) {
				secstruct( i-1, 'L' );
			}
			readE = false;
			lenE = 0;
			continue;
		}
		if( ss =="E" ) {
			readE = true;
			lenE++;
		}
		
		if( readH && ss != "H" ) {
			if( lenH <= 3 ) {
				for( Size j=1; j<=lenH; j++ ) {
					secstruct( i-j, 'L' );
				}
			}
			readH = false;
			lenH = 0;
			continue;
		}
		if( ss =="H" ) {
			readH = true;
			lenH++;
		}
	}

//	bool flag( false );
//	for( Strands::iterator iter = strands_.begin(),
//		iter_end = strands_.end(); iter != iter_end; ++iter ) {
//		StrandOP & st( *iter );
		
//		if( st->end() > secstruct().length() - 2 ) continue;

//		Size ii( 1 );
		/// redefine C-terminal side of strand 
//		Size res( st->end() + ii )
//		while( abego( res ) == "B" || abego( res ) == "Z" || abego( res ) == "Y" ||
//			   abego( res ) == "S" || abego( res ) == "P" ) {
//			ii ++;
//			res = st->end() + ii;
//			if( secstruct().length() < res ) break;
//		}
//		if( ii > 2 ) {
//			flag = true;
//			for( Size jj=1; jj<=ii-2; ++jj ) {
//				secstruct( st->end() + jj, 'E' );
//			}
//		}
		
		/// redefine N-terminal side of strand
//		if( st->begin() < 3 ) continue;
		
//		ii = 1;
//		res = st->begin() - ii;
//		while( abego( res ) == "B" || abego( res ) == "Z" || abego( res ) == "Y" ||
//			   abego( res ) == "S" || abego( res ) == "P" ) {
//			ii ++;
//			res = st->begin() - ii;
//			if( res < 1 ) break;
//		}
		
//		if( ii > 2 ) {
//			flag = true;
//			for( Size jj=1; jj<=ii-2; ++jj ) {
//				secstruct( st->begin() - jj, 'E' );
//			}
//		}
//	} // Strands
	
	// clear only parts of data, since bbpos and abego have to be intact
	strands_.clear();
	helices_.clear();
	strand_id_.clear();
	helix_id_.clear();
	loop_id_.clear();
	ss_element_id_.clear();
	identify_ss( secstruct() );
	set_SSorient();
}

/// @brief set secondary structure into pose
void
SS_Info2::set_ss_into_pose( Pose & pose )
{
	runtime_assert( pose.total_residue() == secstruct().length() );
	for( Size i=1; i<= pose.total_residue(); i++ ){
		char ss = secstruct( i );
		pose.set_secstruct( i, boost::lexical_cast<char>( ss ) );
	}
}
	

} // namespace topology
} // namespace fldsgn
} // namespace protocols
