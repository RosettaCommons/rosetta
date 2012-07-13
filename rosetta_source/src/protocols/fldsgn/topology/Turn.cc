 // -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// Copyright in the Rosetta software belongs to the developers and their institutions.
// For more information, see www.rosettacommons.org.

/// @file ./src/protocols/fldsgn/topology/Turn.cc
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@u.washington.edu )

// unit headers
#include <protocols/fldsgn/topology/Turn.hh>

// Project headers
#include <basic/Tracer.hh>
#include <protocols/fldsgn/topology/StrandPairing.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>

// utility headers
#include <utility/exit.hh>

// C++ headers
#include <iostream>
#include <map>

#include <utility/vector1.hh>
#include <utility/options/IntegerVectorOption.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>

static basic::Tracer TR( "protocols.topology.Turn" );

using namespace core;

namespace protocols {
namespace fldsgn {
namespace topology {


/// @brief default constructor
Turn::Turn() :
	ReferenceCount(),
	type_( "" ),
	length_( 0 ),
	begin_( 0 ),
	end_( 0 ),
	abego_( "" ),
	pleat_( 0 )
{}

/// @brief value constructor
Turn::Turn( SS_Info2_COP const ssinfo, StrandPairing const & sp ) :
ReferenceCount()
{
	initialize( ssinfo, sp );
}
	
/// @brief copy constructor
Turn::Turn( Turn const & s ):
	ReferenceCount(),
	type_( s.type_ ),
	length_( s.length_ ),
	begin_( s.begin_ ),
	end_( s.end_ ),
	abego_( s.abego_ ),
	pleat_( s.pleat_ )
{}

/// @brief get pleating
Size
Turn::calc_pleat( SS_Info2_COP const ssinfo, Size const r1, Size const r2 ) const
{

	Vector vCbCa = ( ssinfo->bb_pos().CB( r1 ) - ssinfo->bb_pos().CA( r1 ) ).normalized();
	Vector vCaCa = ( ssinfo->bb_pos().CA( r2 ) - ssinfo->bb_pos().CA( r1 ) ).normalized();
	Vector vNC   = ( ssinfo->bb_pos().C( r1 )  - ssinfo->bb_pos().N( r1 )  ).normalized();

	Real dot = vNC.cross( vCaCa ).dot( vCbCa );
	if( r1 < r2 ) { 
		if ( dot < 0 ) {
			return 2;
		} else {
			return 1;			
		}
	} else {
		if ( dot > 0 ) {
			return 1;
		} else {
			return 2;			
		}		
	}
}

/// @brief initialize this class
void
Turn::initialize( SS_Info2_COP const ssinfo, StrandPairing const & sp )
{
	type_   = "TEST";
	length_ = sp.end2() - sp.end1() - 1;
	begin_  = sp.end1() + 1;
	end_	= sp.end2() - 1;

	for( Size ii=begin_; ii<=end_; ++ii ) {
		abego_ += ssinfo->abego( ii );
	}

	core::chemical::ResidueType aa1 = *ssinfo()->bb_pos().residue_type( sp.end1() );
	core::chemical::ResidueType aa2 = *ssinfo()->bb_pos().residue_type( sp.end2() );
	if( aa1.aa() != core::chemical::aa_gly ) {
		pleat_ = calc_pleat( ssinfo, sp.end1(), sp.end2() );		
	} else {
		if( aa2.aa() != core::chemical::aa_gly ) {
			pleat_ = calc_pleat( ssinfo, sp.end2(), sp.end1() );		
		} else {
			pleat_ = 3;
		}
	}		
}

/// @brief default destructor
Turn::~Turn(){}

/// @brief print this
std::string	
Turn::get_info() const
{
	std::ostringstream info;
	info << type_ << " " << abego_ << " " << begin_ << "-" << end_ << " " << pleat_;
	return info.str();
}
	
/// @brief return strand pairing
std::ostream & operator<<( std::ostream & out, const Turn &s )
{
	out << s.type_ << ", " << s.begin_ << "-" << s.end_ << ", " << s.pleat_ << ", " << s.abego_;
	return out;
}

/////////////////////////////////////////////////////////////////////////////////////////////
/// @brief default constructor
TurnSet::TurnSet():
	max_length_turn_( 10 )
{}

/// @brief value constructor
TurnSet::TurnSet( Turns const & turns ):
	turns_( turns ),
	max_length_turn_( 10 )
{}
	
/// @brief value constructor
TurnSet::TurnSet( SS_Info2_COP const ssinfo, StrandPairingSetCOP const spairset ):
max_length_turn_( 10 )
{
	initialize( ssinfo, spairset );
}

/// @brief copy constructor
TurnSet::TurnSet( TurnSet const & s ) :
	ReferenceCount(),
	max_length_turn_( s.max_length_turn_ ),
	map_sspair_turn_( s.map_sspair_turn_ ),
	turns_( s.turns_ )
{}

/// @brief output turn information
std::ostream & operator<<( std::ostream & out, const TurnSet &s )
{
	out << "#### Turn Info " << std::endl;
	for( Size i=1; i<=s.size() ; i++ ){
		TurnCOP top ( s.turn( i ) );
		out << "## Turn " << i << ":";
		out << (*top) << std::endl;
	}	
	return out;
}

/// @brief add turn_ into turns_
void
TurnSet::push_back( TurnOP const sop )
{
	turns_.push_back( sop );
}

/// @brief clear data
void TurnSet::clear()
{
	turns_.clear();
}

/// @brief return turn_
TurnOP
TurnSet::turn( Size const s ) const
{
	runtime_assert( s <= turns_.size() );
	return turns_[ s ];
}
	
/// @brief return TurnOP	
TurnOP
TurnSet::sspair_turn( Size const s1, Size const s2 ) const
{
	return map_sspair_turn_[ s1 ][ s2 ];		
}


/// @brief initialize TurnSet
void
TurnSet::initialize( SS_Info2_COP const ssinfo, StrandPairingSetCOP const spairset )
{
	// do nothing if strand pairs < 2
	if( spairset->num_strands() < 2 ) {
		return;
	}

	// check spairset was finalized or not
	runtime_assert( spairset->finalized() );
	
	Size num_strands = ssinfo->strands().size();

	// intialize map strand pairings
	map_sspair_turn_.resize( num_strands );
	for( Size i=1; i<=num_strands; i++ ) {
		map_sspair_turn_[i].resize( num_strands );
		for( Size j=1; j<=num_strands; j++ ) {
			map_sspair_turn_[i][j] = 0;
		}
	}
		
	for( Size ii=1; ii<=spairset->size(); ++ii ) {
		StrandPairing spair( *spairset->strand_pairing( ii ) );
		if( spair.orient() == 'A' && (spair.s2() - spair.s1()) == 1 ) {
			Size len = spair.end2() - spair.end1() - 1;
			if( len <= max_length_turn_ ) {
				TurnOP top( new Turn( ssinfo, spair ) );
				turns_.push_back( top );
				map_sspair_turn_[ spair.s1() ][ spair.s2() ] = top;
				map_sspair_turn_[ spair.s2() ][ spair.s1() ] = top;
			}
		}		
	}	
} // initialize


} // namespace topology
} // namespace fldsgn
} // namespace protocols

