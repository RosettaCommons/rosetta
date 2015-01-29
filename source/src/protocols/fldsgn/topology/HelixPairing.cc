// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ./src/protocols/fldsgn/topology/HelixPairing.cc
/// @brief class for helix pairings
/// @author Nobuyasu Koga ( nobuyasu@u.washington.edu )

// unit headers
#include <protocols/fldsgn/topology/HelixPairing.hh>

// project headers
#include <protocols/fldsgn/topology/SS_Info2.hh>

// utility headers
#include <utility/string_util.hh>
#include <utility/exit.hh>

// C++ headers
#include <utility/assert.hh>
#include <iostream>
#include <sstream>
#include <boost/lexical_cast.hpp>
#include <basic/Tracer.hh>

#include <ObjexxFCL/format.hh>

// numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>

#include <utility/vector1.hh>


static thread_local basic::Tracer TR( "protocols.topology.HelixPairing" );

using namespace core;
typedef std::string String;

namespace protocols {
namespace fldsgn {
namespace topology {

/// @Brief default constructor
HelixPairing::HelixPairing():
	h1_( 0 ),
	h2_( 0 ),
	orient_( 'N' ),
	name_( "" ),
	dist_( 0.0 ),
	cross_angle_( 0.0 ),
	align_angle_( -99.0 ),
	loop_length_( 6 )
{}

/// @brief value constructor
HelixPairing::HelixPairing(
  Size const h1,
	Size const h2,
	char const o
):
	h1_( h1 ),
	h2_( h2 ),
	orient_( o ),
	dist_( 0.0 ),
	cross_angle_( 0.0 ),
	align_angle_( -99.0 ),
	loop_length_( 6 )
{
	runtime_assert( h1 < h2 );
	initialize();
}

/// @brief value constructor
HelixPairing::HelixPairing( String const & hp ):
	dist_( 0.0 ),
	cross_angle_( 0.0 ),
	align_angle_( -99.0 ),
	loop_length_( 6 )
{
	utility::vector1< String > parts( utility::string_split( hp, '.' ) );
	runtime_assert( parts.size() == 2 );

	utility::vector1< String > helices( utility::string_split( parts[1], '-' ) );
	h1_ = boost::lexical_cast<Size>( helices[1] );
	h2_ = boost::lexical_cast<Size>( helices[2] );
	runtime_assert( h1_ < h2_ );

	char para = parts[2][0];
	runtime_assert( para == 'P' || para == 'A' );
	orient_ = parts[2][0];

	name_ = hp;
}


void HelixPairing::initialize()
{
	using namespace boost;
	name_ = lexical_cast<String>(h1_) + '-' + lexical_cast<String>(h2_) + '.' + orient_ ;
}

/// @brief copy constructor
HelixPairing::HelixPairing( HelixPairing const & hp ) :
	ReferenceCount(),
	h1_( hp.h1_ ),
	h2_( hp.h2_ ),
	orient_( hp.orient_ ),
	name_( hp.name_ ),
	dist_( hp.dist_ ),
	cross_angle_( hp.cross_angle_ ),
	align_angle_( hp.align_angle_ ),
	loop_length_( hp.loop_length_ )
{}

/// @brief default destructor
HelixPairing::~HelixPairing(){}

/// @brief clone this object
HelixPairingOP HelixPairing::clone()
{
	return HelixPairingOP( new HelixPairing( *this ) );
}

/// @brief return name
std::ostream & operator<<(std::ostream & out, const HelixPairing &hp )
{
  using ObjexxFCL::format::F;
  using ObjexxFCL::format::A;
	out << A( 7, hp.name() ) << " "
			<< F( 8, 3, hp.dist() ) << F( 8, 3, hp.cross_angle() ) << F( 8, 3, hp.align_angle() );
	//out << A( 7, hp.name() ) << " "
	//<< F( 8, 3, hp.dist() ) << F( 8, 3, hp.cross_angle() );
	return out;
}

/// @brief
bool HelixPairing::is_parallel() const {
	if( orient_ == 'P' ) {
		return true;
	}else{
		return false;
	}
}

/// @brief
void
HelixPairing::calc_geometry( SS_Info2_COP const ss_info )
{

	using protocols::fldsgn::topology::Helix;
	using protocols::fldsgn::topology::Helices;
	using protocols::fldsgn::topology::Strands;

	Real flip( 1.0 );
	if( orient() == 'A' ) flip = -1.0;

	Helices const & helices( ss_info->helices() );
	Strands const & strands( ss_info->strands() );

	Helix const & hx1 = *helices[ h1() ];
	Helix const & hx2 = *helices[ h2() ];
	Vector const midvec = hx2.mid_pos() - hx1.mid_pos();

	dist_ = midvec.length();
	cross_angle_ = numeric::conversions::degrees( angle_of( hx1.orient(), flip*hx2.orient() ) );

	int r1 = int( hx1.begin() ) - int( loop_length_ );
	int r2 = int( hx2.begin() ) - int( loop_length_ );
	if( r1 < 1 || r2 < 1 ) return;

	Size const s1 = ss_info->strand_id( Size( r1 ) );
	Size const s2 = ss_info->strand_id( Size( r2 ) );
	if( s1 == 0 || s2 == 0 ) return;

	// check the strand pairing of s1 and s2 is parallel, otherwise skip this calculation
	Real dot = ss_info->strand( s1 )->orient().dot( ss_info->strand( s2 )->orient() );
	Real sign( 1.0 );
	if( dot < 0 ) sign = -1.0;

	Vector const v1 = ( strands[ s2 ]->mid_pos() - strands[ s1 ]->mid_pos() ).normalized();
	Vector const v2 = ( strands[ s1 ]->orient() + sign*strands[ s2 ]->orient() ).normalized();
	Vector sheet_plane = v1.cross( v2 );

	Vector const h1 = hx1.Nend_pos() - sheet_plane.dot( hx1.Nend_pos() - strands[ s2 ]->mid_pos() )*sheet_plane;
	Vector const h2 = hx1.Cend_pos() - sheet_plane.dot( hx1.Cend_pos() - strands[ s2 ]->mid_pos() )*sheet_plane;
	Vector const h3 = hx2.Nend_pos() - sheet_plane.dot( hx2.Nend_pos() - strands[ s2 ]->mid_pos() )*sheet_plane;
	Vector const h4 = hx2.Cend_pos() - sheet_plane.dot( hx2.Cend_pos() - strands[ s2 ]->mid_pos() )*sheet_plane;

	//Vector const v3 = hx1.orient().project_parallel( v1 ) + hx1.orient().project_parallel( v2 );
	//Vector const v4 = hx2.orient().project_parallel( v1 ) + hx2.orient().project_parallel( v2 );
	Vector const v3 = h2 - h1;
	Vector const v4 = h4 - h3;

	align_angle_ = numeric::conversions::degrees( angle_of( v3, flip*v4 ) );

	// TR.Debug << h1() << " " << h2() << " " << dist_ << " " << cross_angle_ << " " << align_angle_ << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief default constructor
HelixPairingSet::HelixPairingSet():
	hpairset_name_( "" ),
	num_helices_( 0 ),
	initialize_map_helix_pairings_( false )
{}

/// @brief value constructor
HelixPairingSet::HelixPairingSet( String const & helix_pairings ):
	hpairset_name_( helix_pairings ),
	num_helices_( 0 ),
	initialize_map_helix_pairings_( false )
{
	if( helix_pairings == "" ) {
		return;
	}

	utility::vector1< String > hpairs( utility::string_split( helix_pairings, ';' ) );
	for( utility::vector1< String >::const_iterator iter = hpairs.begin(); iter != hpairs.end() ; ++iter) {
		helix_pairings_.push_back( protocols::fldsgn::topology::HelixPairingOP( new HelixPairing( *iter ) ) );
	}
}

/// @brief value constructor
HelixPairingSet::HelixPairingSet( HelixPairings const & helix_pairings ):
	helix_pairings_( helix_pairings ),
	hpairset_name_( "" ),
	num_helices_( 0 ),
	initialize_map_helix_pairings_( false )
{
	// set hpairset_name_
	for ( HelixPairings::const_iterator it=helix_pairings_.begin(),	ite=helix_pairings_.end(); it != ite; ++it ) {
		HelixPairing const & hpair( **it );
		if( hpairset_name_ == "" ){
			hpairset_name_ = hpair.name();
		}else{
			hpairset_name_ += ';' + hpair.name();
		}
	}
}

/// @brief copy constructor
HelixPairingSet::HelixPairingSet( HelixPairingSet const & s ):
 	ReferenceCount(),
	helix_pairings_( s.helix_pairings_ ),
	hpairset_name_( s.hpairset_name_ ),
	num_helices_( s.num_helices_ ),
	initialize_map_helix_pairings_( s.initialize_map_helix_pairings_ ),
	map_helix_pairings_( s.map_helix_pairings_ )
{}

/// @brief destructor
HelixPairingSet::~HelixPairingSet(){}

/// @brief clone this object
HelixPairingSetOP
HelixPairingSet::clone() const
{
	return HelixPairingSetOP( new HelixPairingSet( *this ) );
}

/// @brief return helix pairing
std::ostream & operator<<( std::ostream & out, const HelixPairingSet &s )
{
 	out << "#### HelixPairingSet Info " << std::endl;
	out << "# " << s.name() << std::endl;
	out << "# name distance cross_angle align_angle " << std::endl;
	HelixPairings const & hpairs( s.helix_pairings() );
	for( HelixPairings::const_iterator iter = hpairs.begin(); iter != hpairs.end(); ++iter ) {
		out << "# " << (**iter) << std::endl;
  }
	return out;
}

/// @brief add HelixPairingOP to StandPairingSet
void
HelixPairingSet::push_back( HelixPairingOP const hop )
{
	helix_pairings_.push_back( hop );
	initialize_map_helix_pairings_ = false;
	if( hpairset_name_ == "" ){
		hpairset_name_ = hop->name();
	}else{
		hpairset_name_ += ';' + hop->name();
	}
}

/// @brief clear data of this HelixParingSet
void
HelixPairingSet::clear()
{
	helix_pairings_.clear();
	hpairset_name_  = "";
	initialize_map_helix_pairings_ = false;
	map_helix_pairings_.clear();
}

/// @brief return one of the stand_pairings give a number
HelixPairingOP
HelixPairingSet::helix_pairing( Size const s ) const
{
	runtime_assert( s <= helix_pairings_.size() );
	return helix_pairings_[ s ];
}

/// @brief
HelixPairingOP
HelixPairingSet::helix_pairing( Size const h1, Size const h2 )
{
	if( ! initialize_map_helix_pairings_ ){
		create_map_helix_pairings();
	}

	if( h1 <= num_helices_ && h2 <= num_helices_ ) {
		return map_helix_pairings_[ h1 ][ h2 ];
	} else {
		return 0;
	}
}

/// @brief return all helix pairings
HelixPairings const &
HelixPairingSet::helix_pairings() const
{
	return helix_pairings_;
}

/// @brief return the size of helix_pairings_
Size
HelixPairingSet::size() const
{
	return helix_pairings_.size();
}

/// @brief return hpairset_name_
String
HelixPairingSet::name() const
{
	return hpairset_name_;
}

/// @brief calculate geomtry of helix parings
void
HelixPairingSet::calc_geometry( SS_Info2_COP const ss_info )
{
	for( HelixPairings::const_iterator iter = helix_pairings_.begin(); iter != helix_pairings_.end(); ++iter ) {
		HelixPairingOP const hpair(*iter);
		hpair->calc_geometry( ss_info );
  }
}

/// @brief create 2D table of helix pairings
void HelixPairingSet::create_map_helix_pairings()
{
	initialize_map_helix_pairings_ = true;
	for ( HelixPairings::const_iterator it=helix_pairings_.begin(),	ite=helix_pairings_.end(); it != ite; ++it ) {
		HelixPairing hpair(**it);
		if(	hpair.h2() > num_helices_ ){
			num_helices_ = hpair.h2();
		}
	}
	map_helix_pairings_.resize( num_helices_ );
	for( Size i=1; i<=num_helices_; i++ ){
		map_helix_pairings_[i].resize( num_helices_ );
		for( Size j=1; j<=num_helices_; j++ ){
			map_helix_pairings_[i][j] = 0;
		}
	}

	for ( HelixPairings::const_iterator it=helix_pairings_.begin(),	ite=helix_pairings_.end(); it != ite; ++it ) {
		HelixPairingOP const hpair( *it );
		map_helix_pairings_[ hpair->h1() ][ hpair->h2() ] = hpair;
		map_helix_pairings_[ hpair->h2() ][ hpair->h1() ] = hpair;
	}
	//TR << "#helices = " << num_helices_ << std::endl;
}


} // namespace topology
} // namespace fldsgn
} // namespace protocols
