// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/fldsgn/topology/SS_Info2.cc
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@u.washington.edu )

/// Unit headers
#include <protocols/fldsgn/topology/SS_Info2.hh>

/// Project headers
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <basic/datacache/CacheableData.hh>

/// C++ Headers
#include <iostream>
#include <utility>
#include <utility/assert.hh>

// Numeric
#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>

#include <utility/vector1.hh>


using namespace core;

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace fldsgn {
namespace topology {

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SS_Base
/// @brief default constructor
SS_Base::SS_Base():
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
SS_Base::SS_Base( Size const & begin, Size const & end ):
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
SS_Base::SS_Base( SS_Base const & s ):
	ReferenceCount(),
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
SS_Base::~SS_Base()= default;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief default constructor
Strand::Strand():
	SS_Base()
{}


/// @brief value constructor
Strand::Strand( Size const & begin, Size const & end ):
	SS_Base( begin, end )
{}


/// @brief copy constructor
Strand::Strand( Strand const & /*s*/ ) = default;


/// @brief destructor
Strand::~Strand()= default;


/// @brief
std::ostream & operator<<(std::ostream & out, const Strand & st )
{
	out << st.begin() << "-" << st.end() ;
	return out;
}


/// @brief set vector between Calpha atoms of edge residues
void
Strand::calc_geometry( BB_Pos const & bbpos )
{
	is_geometry_initialized( true );

	// orient is defined by begin and end Ca atoms of strand
	Nend_pos( bbpos.N( begin() ) );
	Cend_pos( bbpos.C( end() ) );
	Vector v = Cend_pos() - Nend_pos();
	orient( v.normalized() );

	Size pos = ( end() - begin() )/2 + begin();
	if ( ( end() - begin() )%2 == 0 ) {
		mid_pos( ( bbpos.N( pos ) + bbpos.C( pos ) )/2.0 );
	} else {
		mid_pos( ( bbpos.N( pos ) + bbpos.C( pos+1 ) )/2.0 );
	}

	Nend_orient( ( mid_pos() - Nend_pos() ).normalized() );
	Cend_orient( ( Cend_pos() - mid_pos() ).normalized() );
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief default constructor
Helix::Helix():
	SS_Base(),
	bend_( 0 )
{}


/// @brief value constructor
Helix::Helix( Size const & begin, Size const & end ):
	SS_Base( begin, end ),
	bend_( 0 )
{}


/// @brief copy constructor
Helix::Helix( Helix const & /*s*/ ) = default;


/// @brief destructor
Helix::~Helix()= default;


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
	SS_Base(),
	type_("")
{}


/// @brief value constructor
Loop::Loop( Size const begin, Size const end, String const & type ):
	SS_Base( begin, end ),
	type_( type )
{}


/// @brief copy constructor
Loop::Loop( Loop const & /*s*/ ) = default;


/// @brief destructor
Loop::~Loop()= default;


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
	secstruct_( "" )
{}


/// @brief value constructor
SS_Info2::SS_Info2( String const & secstruct ):
	CacheableData()
{
	initialize( secstruct );
}


/// @brief value constructor
SS_Info2::SS_Info2( Pose const & pose, String const & secstruct ):
	CacheableData()
{
	initialize( pose, secstruct );
}


/// @brief copy constructor
SS_Info2::SS_Info2( SS_Info2 const & s ):
	CacheableData(),
	secstruct_( s.secstruct_ ),
	bb_pos_( s.bb_pos_ ),
	strands_( s.strands_ ),
	strand_id_( s.strand_id_ ),
	helices_( s.helices_ ),
	helix_id_( s.helix_id_ ),
	loop_id_( s.loop_id_ ),
	ss_element_id_( s.ss_element_id_ )
{}


/// @brief destructor
SS_Info2::~SS_Info2()= default;


/// @brief make clone
basic::datacache::CacheableDataOP
SS_Info2::clone() const
{
	return basic::datacache::CacheableDataOP( new SS_Info2( *this ) );
}


/// @brief
void
SS_Info2::clear_data()
{
	secstruct_ = "";
	bb_pos_.clear();
	strands_.clear();
	helices_.clear();
	strand_id_.clear();
	helix_id_.clear();
	loop_id_.clear();
	ss_element_id_.clear();
}


/// @brief resize vectors
void
SS_Info2::resize( Size const nres )
{
	bb_pos_.resize( int(nres) );
	helix_id_.resize( nres );
	strand_id_.resize( nres );
	loop_id_.resize( nres );
	ss_element_id_.resize( nres );
}


/// @brief initialize parameters of this class
void
SS_Info2::initialize( Pose const & pose, String const & secstruct )
{
	bbpos_is_set_ = true;

	// if the pose has ligands, they will be considered 'L'
	/*
	//flo sep'12
	//not clear what to do if the pose has ligands
	//ideally the class should be initalized with the ligand positions,
	//but obviously the ligands don't have any secondary structure
	//attempt at getting desired behavior: initialize with pose length,
	//but make sure that secondary structure is according to number
	//of protein res
	core::Size num_protein_res( pose.size() );
	for ( core::Size i = pose.size(); i != 0; i-- ) {
	if ( !pose.residue_type( i ).is_protein() ) num_protein_res--;
	}
	*/

	// data all clear
	clear_data();
	// set strands and helices
	if ( secstruct == "" ) {
		secstruct_ = pose.secstruct();
	} else {
		secstruct_ = secstruct;
	}
	runtime_assert( pose.size() == secstruct_.length() ); //flo sep'12 changed from pose.size() to num_protein_res
	for ( core::Size i=1; i<=pose.size(); ++i ) {
		if ( !pose.residue( i ).is_protein() ) {
			secstruct_[ i - 1 ] = 'L';
		}
	}
	resize( pose.size() );
	identify_ss( secstruct_ );
	// set bb_pos
	bb_pos_.take_coordinates_from_pose( pose );
	// set orient of strands and helices
	set_SSorient();
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
	identify_ss( secstruct_ );
}


/// @brief
std::ostream & operator<<(std::ostream & out, const SS_Info2 & ssinfo )
{
	Size count( 0 );
	out << "#### SS_Info " << std::endl;
	for ( auto const & helice : ssinfo.helices_ ) {
		count ++;
		Helix const & hx( *helice );
		out << "# Helix  " << count << ": " << hx << std::endl;
	}
	count = 0;
	for ( auto const & strand : ssinfo.strands_ ) {
		count ++;
		Strand const & st( *strand );
		out << "# Strand " << count << ": " << st << std::endl;
	}
	return out;
}


/// @brief set orientation vector of secondary structures given a pose
void
SS_Info2::set_SSorient( Pose const & pose )
{
	runtime_assert( pose.size() == bb_pos_.size() );
	bb_pos_.take_coordinates_from_pose( pose );
	bbpos_is_set_ = true;
	set_SSorient();
}


/// @brief set orientation vector of secondary structures given a pose which is defined in the constructor
void
SS_Info2::set_SSorient()
{
	for ( auto & st : strands_ ) {
		if (  st->length() >= 2 ) {
			st->calc_geometry( bb_pos_ );
		}
	} // Strands
	for ( auto & hx : helices_ ) {
		if (  hx->length() >= 4 ) {
			hx->calc_geometry( bb_pos_ );
		}
	} // Helices
}


/// @brief identify the region of each secondary structure element
void
SS_Info2::identify_ss( String const & secstruct )
{
	bool flag_L( false );
	bool flag_E( false );
	bool flag_H( false );
	Size beginE( 0 ), beginH( 0 ), beginL( 0 );
	Size istrand( 0 ), ihelix( 0 ), iloop( 0 ), iss( 0 );

	String prev( "" );
	for ( Size i=1; i<= secstruct.length(); ++i ) {
		String const & ss( secstruct.substr( i-1, 1 ) );

		strand_id_[ i ] = 0;
		helix_id_ [ i ] = 0;
		loop_id_  [ i ] = 0;

		if ( ss =="E" ) {

			if ( flag_E == false ) {
				istrand ++;
				beginE = i;
			}
			flag_E = true;
			strand_id_[ i ] = istrand;

		} else if ( ss == "H" ) {

			if ( flag_H == false ) {
				ihelix++;
				beginH = i;
			}
			flag_H = true;
			helix_id_[ i ] = ihelix;

		} else {

			if ( flag_L == false ) {
				iloop++;
				beginL = i;
			}
			flag_L = true;
			loop_id_[ i ] = iloop;

		}

		if ( ss !="E" && flag_E == true ) {
			flag_E = false;
			if ( ( i - beginE ) >= 2 ) {
				strands_.push_back( StrandOP( new Strand( beginE, i-1 ) ) );
			} else {
				strands_.push_back( StrandOP( new Strand( beginE, beginE ) ) );
			}
		}

		if ( ss !="H" && flag_H == true ) {
			flag_H = false;
			if ( ( i-beginH ) >= 2 ) {
				helices_.push_back( HelixOP( new Helix( beginH, i-1 ) ) );
			} else {
				helices_.push_back( HelixOP( new Helix( beginH, beginH ) ) );
			}
		}

		if ( ss !="L" && flag_L == true ) {
			flag_L = false;
			if ( ( i-beginL ) >= 2 ) {
				loops_.push_back( LoopOP( new Loop( beginL, i-1 ) ) );
			} else {
				loops_.push_back( LoopOP( new Loop( beginL, beginL ) ) );
			}
		}

		if ( prev != ss ) {
			iss ++;
			prev = ss;
		}
		ss_element_id_[ i ] = iss;

	} // for( Size i )

	if ( flag_E == true ) {
		if ( ( secstruct.length() - beginE + 1 ) >= 2 ) {
			strands_.push_back( StrandOP( new Strand( beginE, secstruct.length() ) ) );
		} else {
			strands_.push_back( StrandOP( new Strand( beginE, beginE ) ) );
		}
	}
	if ( flag_H == true ) {
		if ( ( secstruct.length() - beginH + 1 ) >= 2 ) {
			helices_.push_back( HelixOP( new Helix( beginH, secstruct.length() ) ) );
		} else {
			helices_.push_back( HelixOP( new Helix( beginH, beginH ) ) );
		}
	}
	if ( flag_L == true ) {
		if ( ( secstruct.length() - beginL + 1 ) >= 2 ) {
			loops_.push_back( LoopOP( new Loop( beginL, secstruct.length() ) ) );
		} else {
			loops_.push_back( LoopOP( new Loop( beginL, beginL ) ) );
		}
	}

} // SS_Info2::identify_ss()

} // namespace topology
} // namespace fldsgn
} // namespace protocols



#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::fldsgn::topology::SS_Base::save( Archive & arc ) const {
	arc( CEREAL_NVP( begin_ ) ); // Size
	arc( CEREAL_NVP( end_ ) ); // Size
	arc( CEREAL_NVP( is_geometry_initialized_ ) ); // _Bool
	arc( CEREAL_NVP( orient_ ) ); // Vector
	arc( CEREAL_NVP( Nend_orient_ ) ); // Vector
	arc( CEREAL_NVP( Cend_orient_ ) ); // Vector
	arc( CEREAL_NVP( Nend_pos_ ) ); // Vector
	arc( CEREAL_NVP( Cend_pos_ ) ); // Vector
	arc( CEREAL_NVP( mid_pos_ ) ); // Vector
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::fldsgn::topology::SS_Base::load( Archive & arc ) {
	arc( begin_ ); // Size
	arc( end_ ); // Size
	arc( is_geometry_initialized_ ); // _Bool
	arc( orient_ ); // Vector
	arc( Nend_orient_ ); // Vector
	arc( Cend_orient_ ); // Vector
	arc( Nend_pos_ ); // Vector
	arc( Cend_pos_ ); // Vector
	arc( mid_pos_ ); // Vector
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::fldsgn::topology::SS_Base );
CEREAL_REGISTER_TYPE( protocols::fldsgn::topology::SS_Base )


/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::fldsgn::topology::Loop::save( Archive & arc ) const {
	arc( cereal::base_class< class protocols::fldsgn::topology::SS_Base >( this ) );
	arc( CEREAL_NVP( type_ ) ); // String
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::fldsgn::topology::Loop::load( Archive & arc ) {
	arc( cereal::base_class< class protocols::fldsgn::topology::SS_Base >( this ) );
	arc( type_ ); // String
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::fldsgn::topology::Loop );
CEREAL_REGISTER_TYPE( protocols::fldsgn::topology::Loop )


/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::fldsgn::topology::Helix::save( Archive & arc ) const {
	arc( cereal::base_class< class protocols::fldsgn::topology::SS_Base >( this ) );
	arc( CEREAL_NVP( bend_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::fldsgn::topology::Helix::load( Archive & arc ) {
	arc( cereal::base_class< class protocols::fldsgn::topology::SS_Base >( this ) );
	arc( bend_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::fldsgn::topology::Helix );
CEREAL_REGISTER_TYPE( protocols::fldsgn::topology::Helix )


/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::fldsgn::topology::Strand::save( Archive & arc ) const {
	arc( cereal::base_class< class protocols::fldsgn::topology::SS_Base >( this ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::fldsgn::topology::Strand::load( Archive & arc ) {
	arc( cereal::base_class< class protocols::fldsgn::topology::SS_Base >( this ) );
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::fldsgn::topology::Strand );
CEREAL_REGISTER_TYPE( protocols::fldsgn::topology::Strand )


/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::fldsgn::topology::SS_Info2::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( bbpos_is_set_ ) ); // _Bool
	arc( CEREAL_NVP( secstruct_ ) ); // String
	arc( CEREAL_NVP( bb_pos_ ) ); // BB_Pos
	arc( CEREAL_NVP( strands_ ) ); // Strands
	arc( CEREAL_NVP( strand_id_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( helices_ ) ); // Helices
	arc( CEREAL_NVP( helix_id_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( loops_ ) ); // Loops
	arc( CEREAL_NVP( loop_id_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( ss_element_id_ ) ); // utility::vector1<Size>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::fldsgn::topology::SS_Info2::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( bbpos_is_set_ ); // _Bool
	arc( secstruct_ ); // String
	arc( bb_pos_ ); // BB_Pos
	arc( strands_ ); // Strands
	arc( strand_id_ ); // utility::vector1<Size>
	arc( helices_ ); // Helices
	arc( helix_id_ ); // utility::vector1<Size>
	arc( loops_ ); // Loops
	arc( loop_id_ ); // utility::vector1<Size>
	arc( ss_element_id_ ); // utility::vector1<Size>
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::fldsgn::topology::SS_Info2 );
CEREAL_REGISTER_TYPE( protocols::fldsgn::topology::SS_Info2 )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_fldsgn_topology_SS_Info2 )
#endif // SERIALIZATION



