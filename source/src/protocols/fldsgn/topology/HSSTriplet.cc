// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ./src/protocols/src/fldsgn/topology/HSSTriplet.cc
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@u.washington.edu )

// unit headers
#include <protocols/fldsgn/topology/HSSTriplet.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>

#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <ostream>
#include <utility/string_util.hh>
#include <boost/lexical_cast.hpp>
#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>

#include <utility/vector1.hh>


static thread_local basic::Tracer TR( "protocols.topology.HSSTriplet" );

namespace protocols {
namespace fldsgn {
namespace topology {

// @brief Auto-generated virtual destructor
HSSTriplet::~HSSTriplet() {}


/// @value constructor
HSSTriplet::HSSTriplet( String const & hss )
{
	utility::vector1< String > parts( utility::string_split( hss, ',' ) );
	runtime_assert( parts.size() == 2 );
	helix_ = boost::lexical_cast<Size>( parts[1] );

	utility::vector1< String > st( utility::string_split( parts[2], '-' ) );
	runtime_assert( st.size() == 2 );
	strand1_ = boost::lexical_cast<Size>( st[1] );
	strand2_ = boost::lexical_cast<Size>( st[2] );
}


/// @copy constructor
HSSTriplet::HSSTriplet( HSSTriplet const & hss ):
	ReferenceCount(),
	helix_( hss.helix_ ),
	strand1_( hss.strand1_ ),
	strand2_( hss.strand2_ )
{}

/// @brief IO Operator
std::ostream & operator<<(std::ostream & out, const HSSTriplet & s )
{
	out << "Helix:" << s.helix() << " Strand1:" << s.strand1() << " Strand2:" << s.strand2();
	return out;
}

/// @brief
core::Real
HSSTriplet::hs1_dist() const
{
	return hs1_dist_;
}

/// @brief
core::Real
HSSTriplet::hs2_dist() const
{
	return hs2_dist_;
}

/// @brief
core::Real
HSSTriplet::ss_dist() const
{
	return ss_dist_;
}

/// @brief reutrn distance between sheet ( defined by the 2 strands ) and helix
core::Real
HSSTriplet::hsheet_dist() const
{
	return hsheet_dist_;
}

/// @brief return distance between sheet ( defined by the 2 strands ) and helix
core::Real
HSSTriplet::hs_angle() const
{
	return hs_angle_;
}

/// @brief orientation between strands
std::string
HSSTriplet::ss_orient() const
{
	return ss_orient_;
}

/// @brief orientation between helix and 1st strand
std::string
HSSTriplet::hs1_orient() const
{
	return hs1_orient_;
}

/// @brief orientation between helix and 2nd strand
std::string
HSSTriplet::hs2_orient() const
{
	return hs2_orient_;
}

/// @brief orientation between helix and 2nd strand
bool
HSSTriplet::left_handed() const
{
	return left_handed_;
}

/// @brief calc geometry
void
HSSTriplet::calc_geometry( SS_Info2_COP const ssinfo )
{

	using protocols::fldsgn::topology::Helix;
	using protocols::fldsgn::topology::Strand;

	Helix  const & hx ( *ssinfo->helices()[ helix() ] );
	Strand const & s1 ( *ssinfo->strands()[ strand1() ] );
	Strand const & s2 ( *ssinfo->strands()[ strand2() ] );

	runtime_assert( s1.length() >= 2 );
	runtime_assert( s2.length() >= 2 );

	Real dot_ss  = s1.orient().dot( s2.orient() );
	Real dot_hs1 = hx.orient().dot( s1.orient() );
	Real dot_hs2 = hx.orient().dot( s2.orient() );

	Real sign_ss/*, sign_hs1, sign_hs2*/;

	if ( dot_ss < 0 ) {
		sign_ss = -1.0;
		ss_orient_ = "A";
	} else {
		sign_ss = 1.0;
		ss_orient_ = "P";
	}

	if ( dot_hs1 < 0 ) {
		//sign_hs1 = -1.0;  // set but never used ~Labonte
		hs1_orient_ = "A";
	} else {
		//sign_hs1 = 1.0;  // set but never used ~Labonte
		hs1_orient_ = "P";
	}

	if ( dot_hs2 < 0 ) {
		//sign_hs2 = -1.0;  // set but never used ~Labonte
		hs2_orient_ = "A";
	} else {
		//sign_hs2 = 1.0;  // set but never used ~Labonte
		hs2_orient_ = "P";
	}

	// get mid point of helix and strands
	Vector const hmid = hx.mid_pos();
	Vector const s1_mid = s1.mid_pos();
	Vector const s2_mid = s2.mid_pos();

	// calc_handness
	Vector const s12_vec = ( s2_mid - s1_mid ).normalized();
	Vector const hs1_vec = ( hmid - s1_mid ).normalized();
	Real d = hs1_vec.cross( s12_vec ).dot( s1.orient() );
	if ( d < 0 ) left_handed_ = true;

	// define beta sheet plane
	Vector const s1_chalf_mid = ( s1.mid_pos() + s1.Cend_pos() )/2.0;
	Vector const s2_nhalf_mid = ( s2.mid_pos() + s2.Nend_pos() )/2.0;
	Vector const v1 = ( s2_nhalf_mid - s1_chalf_mid ).normalized();

	// Vector const v2 = ( s1.orient() + sign*s2.orient() ).normalized();
	Vector const v2 = ( s1.Cend_orient() + sign_ss*s2.Nend_orient() ).normalized();

	Vector sheet_plane = v1.cross( v2 );

	// distance between helix and sheet
	hsheet_dist_ = sheet_plane.normalized().dot( hmid - s1_chalf_mid );

	// distance between helix and strand1
	hs1_dist_ = ( s1.mid_pos() - hmid ).length();

	// distance between helix and strand2
	hs2_dist_ = ( s2.mid_pos() - hmid ).length();

	// distance between strand1 and strand2
	ss_dist_ = ( s2.mid_pos() - s1.mid_pos() ).length();

	// calc angles between strands and helix
	Vector const h1 = hx.Nend_pos() - sheet_plane.dot( hx.Nend_pos() - s2_nhalf_mid )*sheet_plane;
	Vector const h2 = hx.Cend_pos() - sheet_plane.dot( hx.Cend_pos() - s2_nhalf_mid )*sheet_plane;
	Vector const hx_on_sheet = ( h1 - h2 ).normalized();
	Real ori = cross( hx_on_sheet, v2 ).dot( sheet_plane );
	hs_angle_ = numeric::conversions::degrees( angle_of( v2, hx_on_sheet ) );
	if ( ori < 0 ) hs_angle_ = -1*hs_angle_;

	// std::cout << ss_orient_ << " " << hs1_orient_ << " " << hs2_orient_ << " "
	// << hsheet_dist_ << " " << hs1_dist_ << " " << hs2_dist_ << " " << ss_dist_ << " "
	// << hs_angle_ << std::endl;

	geometry_is_initialized_ = true;

	// using namespace ObjexxFCL::format;
	// Size number( 0 );
	// for( int ii=1; ii<=25; ii++ ) {
	//  for( int jj=1; jj<=25; jj++ ) {
	//   ++number;
	//   Vector pos = v1*( ii - 12 ) + v2*( jj - 12 ) + s1_chalf_mid;
	//   std::cout << "ATOM  " << I(5,number) << ' ' << " CA " << ' ' << "ALA" << ' '
	//    << "B" << I(4,number) << "    "
	//    << F(8,3,pos.x()) << F(8,3,pos.y()) << F(8,3,pos.z())
	//    << F(6,2,1.0) << F(6,2,10.0) << std::endl;
	//  }
	// }
	// std::cout << "ATOM  " << I(5,1) << ' ' << " CA " << ' ' << "ALA" << ' '
	// << "C" << I(4,1) << "    "
	// << F(8,3,hmid.x()) << F(8,3,hmid.y()) << F(8,3,hmid.z())
	// << F(6,2,1.0) << F(6,2,10.0) << std::endl;

}


/// @brief default constructor
HSSTripletSet::HSSTripletSet()
{
	clear();
}


/// @brief value constructor
HSSTripletSet::HSSTripletSet( String const & s )
{
	clear();

	if ( s == "" ) {
		return;
	}

	utility::vector1< String > hsss( utility::string_split( s, ';' ) );
	for ( utility::vector1< String >::const_iterator iter = hsss.begin(); iter != hsss.end() ; ++iter ) {
		push_back( HSSTripletOP( new HSSTriplet( *iter ) ) );
	}
}


/// @brief value constructor
HSSTripletSet::HSSTripletSet( HSSTriplets const & s )
{
	for ( HSSTriplets::const_iterator it=s.begin(), ite=s.end(); it!= ite; ++it ) {
		HSSTripletOP const hss( *it );
		push_back( hss );
	}
}


/// @brief copy constructor
HSSTripletSet::HSSTripletSet( HSSTripletSet const & s ) :
	ReferenceCount(),
	hss_triplets_( s.hss_triplets_ ),
	helix2hss_( s.helix2hss_ )
{}


/// @brief destructor
HSSTripletSet::~HSSTripletSet(){}


/// @brief IO Operator
std::ostream & operator<<(std::ostream & out, const HSSTripletSet & s )
{
	out << "#### HSSTriplet Info " << std::endl;
	for ( HSSTriplets::const_iterator iter = s.hss_triplets().begin(),
			iter_end = s.hss_triplets().end(); iter != iter_end; ++iter ) {
		HSSTriplet const & hss( **iter );
		out << hss << std::endl;
	}
	return out;
}


/// @brief value constructor
void
HSSTripletSet::add_hsstriplets( HSSTriplets const & s )
{
	for ( HSSTriplets::const_iterator it=s.begin(), ite=s.end(); it!= ite; ++it ) {
		HSSTripletOP const hss( *it );
		push_back( hss );
	}
}

/// @brief puch back data
void
HSSTripletSet::push_back( HSSTripletOP const hsop )
{
	for ( std::map< Size, HSSTripletOP >::const_iterator it=helix2hss_.begin(),
			ite=helix2hss_.end(); it!=ite ; ++it ) {
		if ( it->first == hsop->helix() ) {
			TR <<  "Helix "  <<  it->first << " is already defined in HSSTriplet. " << std::endl;
			assert( false );
		}
	}

	hss_triplets_.push_back( hsop );
	helix2hss_[ hsop->helix() ] = hsop;
}

/// @brief clear data
void
HSSTripletSet::clear()
{
	hss_triplets_.clear();
	helix2hss_.clear();
}


/// @brief
HSSTripletOP
HSSTripletSet::hss_triplet( Size const helix )
{
	return helix2hss_[ helix ];
}


/// @brief return all data
HSSTriplets const &
HSSTripletSet::hss_triplets() const
{
	return hss_triplets_;
}


} // namespace topology
} // namespace fldsgn
} // namespace protocols
