// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/src/fldsgn/topology/StrandPairing.cc
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@u.washington.edu )

// unit header
#include <protocols/fldsgn/topology/StrandPairing.hh>


// Project Headers
#include <protocols/fldsgn/topology/DimerPairing.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>

// utility headers
#include <utility>
#include <utility/string_util.hh>
#include <utility/stream_util.hh>

// C++ headers
#include <iostream>
#include <boost/lexical_cast.hpp>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>

static basic::Tracer TR( "protocols.topology.StrandPairing" );

using String = std::string;
using VecSize = utility::vector1<core::Size>;

namespace protocols {
namespace fldsgn {
namespace topology {

/// @brief default constructor
StrandPairing::StrandPairing():
	s1_( 0 ),
	s2_( 0 ),
	begin1_( 0 ),
	end1_( 0 ),
	begin2_( 0 ),
	end2_( 0 ),
	rgstr_shift_( 99 ),
	orient_( 'N' ),
	name_( "" ),
	bulges1_(),
	bulges2_()
{}


/// @brief value constructor
StrandPairing::StrandPairing(
	core::Size const s1,
	core::Size const s2,
	core::Size const b1,
	core::Size const b2,
	core::Size const p,
	Real const rs,
	char const o
):
	s1_( s1 ),
	s2_( s2 ),
	begin1_( b1 ),
	end1_( b1 ),
	begin2_( b2 ),
	end2_( b2 ),
	rgstr_shift_( rs ),
	orient_( o ),
	name_( "" ),
	bulges1_(),
	bulges2_()
{
	runtime_assert( s1 < s2 );
	runtime_assert( b1 < b2 );
	pleats1_.push_back( p );
	pleats2_.push_back( p );
	residue_pair_.insert( std::map< core::Size, core::Size >::value_type( b1, b2 ) );
	residue_pair_.insert( std::map< core::Size, core::Size >::value_type( b2, b1 ) );
	initialize();
}


/// @brief value constructor
StrandPairing::StrandPairing(
	core::Size const s1,
	core::Size const s2,
	Real const rs,
	char const o
):
	s1_( s1 ),
	s2_( s2 ),
	begin1_( 0 ),
	end1_( 0 ),
	begin2_( 0 ),
	end2_( 0 ),
	rgstr_shift_( rs ),
	orient_( o ),
	name_( "" ),
	bulges1_(),
	bulges2_()
{
	runtime_assert( s1 < s2 );
	initialize();
}


/// @brief value constructor
StrandPairing::StrandPairing( String const & spair ):
	begin1_( 0 ),
	end1_( 0 ),
	begin2_( 0 ),
	end2_( 0 ),
	name_( spair ),
	bulges1_(),
	bulges2_()
{
	utility::vector1< String > parts( utility::string_split( spair, '.' ) );
	runtime_assert( parts.size() == 3 );

	utility::vector1< String > st( utility::string_split( parts[1], '-' ) );
	s1_ = boost::lexical_cast<core::Size>( st[1] );
	s2_ = boost::lexical_cast<core::Size>( st[2] );
	runtime_assert( s1_ < s2_ );

	orient_ = parts[2][0];
	runtime_assert( orient_ == 'P' || orient_ == 'A' || orient_ == 'N'  );
	rgstr_shift_ = boost::lexical_cast<Real>( parts[3] );
}


/// @brief initialize StrandPairing
void
StrandPairing::initialize()
{
	using namespace boost;
	name_ = std::to_string(s1_) + '-' + std::to_string(s2_) + '.' + orient_ + '.'+ std::to_string(rgstr_shift_);
}


/// @brief default destructor
StrandPairing::~StrandPairing()= default;


/// @brief clone this object
StrandPairingOP
StrandPairing::clone()
{
	return utility::pointer::make_shared< StrandPairing >( *this );
}


/// @brief return name
std::ostream & operator<<( std::ostream & out, const StrandPairing &sp )
{
	out << sp.name();
	out << ": " << sp.begin1() << "-" << sp.end1() << "." << sp.begin2() << "-" << sp.end2();
	return out;
}


/// @brief
bool
StrandPairing::elongate( core::Size const r1, core::Size const r2, core::Size const p1, core::Size const p2 )
{
	runtime_assert( r2 > r1 );

	// set residue pair
	residue_pair_.insert( std::map< core::Size, core::Size >::value_type( r1, r2 ) );
	residue_pair_.insert( std::map< core::Size, core::Size >::value_type( r2, r1 ) );

	if ( begin1_ == 0 && begin2_ == 0 ) {
		begin1_ = r1;
		begin2_ = r2;
		end1_ = r1;
		end2_ = r2;
		pleats1_.push_back( p1 );
		pleats2_.push_back( p2 );
		return true;
	}

	if ( r1 < end1_ ) {
		return false;
	}
	if ( (orient_  ==  'P' && r2 < end2_) || (orient_ == 'A' && r2 > end2_) ) {
		return false;
	}

	// runtime_assert( r1 > end1_ );

	//  if( orient_  ==  'P' ){
	// runtime_assert( r2 > end2_ );
	// }else{
	// runtime_assert( r2 < end2_ );
	// }

	end1_ = r1;
	end2_ = r2;

	/*
	if ( size1() != size2() ) {
	has_bulge_ = true;
	}
	*/

	return true;
}


/// @brief
bool
StrandPairing::add_pair( core::Size const r1, core::Size const r2, char const orient, Real const rgstr )
{
	runtime_assert( r2 > r1 );

	if ( orient_ != orient ) return false;

	residue_pair_.insert( std::map< core::Size, core::Size >::value_type( r1, r2 ) );
	residue_pair_.insert( std::map< core::Size, core::Size >::value_type( r2, r1 ) );

	if ( begin1_ == 0 && begin2_ == 0 ) {
		begin1_ = r1;
		begin2_ = r2;
		end1_ = r1;
		end2_ = r2;
		return true;
	}

	if ( r1 <= begin1_ ) {
		begin1_ = r1;
		rgstr_shift_ = rgstr;
	} else if ( end1_ <= r1 ) {
		end1_ = r1;
	}

	if ( r2 <= begin2_ ) {
		begin2_ = r2;
		rgstr_shift_ = rgstr;
	} else if ( end2_ <= r2 ) {
		end2_ = r2;
	}

	/*
	if ( size1() != size2() ) {
	has_bulge_ = true;
	}
	*/
	return true;
}


/// @brief return length of 1st strand
StrandPairing::Size
StrandPairing::size1() const
{
	if ( end1_ >= begin1_ ) {
		return end1_ - begin1_ + 1;
	} else {
		return begin1_ - end1_ + 1;
	}
}


/// @brief return length of 2nd strand
StrandPairing::Size
StrandPairing::size2() const
{
	if ( end2_ >= begin2_ ) {
		return end2_ - begin2_ + 1;
	} else {
		return begin2_ - end2_ + 1;
	}
}


/// @brief return length of 2nd strand
bool
StrandPairing::is_parallel() const {
	if ( orient_ == 'P' ) {
		return true;
	} else {
		return false;
	}
}

/// @brief whether input residue is included in this StrandPairinge or not
bool
StrandPairing::is_member( core::Size const res ) {

	if ( begin1_ <= end1_ ) {
		if ( begin1_ <= res && res <= end1_ ) return true;
	} else {
		if ( end1_ <= res && res <= begin1_ ) return true;
	}

	if ( begin2_ <= end2_ ) {
		if ( begin2_ <= res && res <= end2_ ) return true;
	} else {
		if ( end2_ <= res && res <= begin2_ ) return true;
	}

	return false;

}

/// @brief return residue pairing
bool
StrandPairing::has_paired_residue( core::Size const res ) const
{
	return ( residue_pair_.find( res ) != residue_pair_.end() );
}

/// @brief return residue pairing
StrandPairing::Size
StrandPairing::residue_pair( core::Size const res ) const
{
	//runtime_assert( (begin1_ <= res && res <= end1_) || (begin2_ <= res && res <= end2_) ||
	//  (end1_ <= res && res <= begin1_) || (end2_ <= res && res <= begin2_) );
	auto it = residue_pair_.find( res );
	runtime_assert( it != residue_pair_.end() );
	return it->second;
}

std::set< core::Size >
compute_bulges( core::Size const strand_begin, core::Size const strand_end, utility::vector1< String > const & abego )
{
	std::set< core::Size > bulges;
	if ( abego.empty() ) return bulges;

	TR << "Computing bulges between " << strand_begin << " and " << strand_end << std::endl;
	// here it is assumed SS=E
	for ( core::Size resid=strand_begin; resid<=strand_end; ++resid ) {
		if ( abego[resid] == "A" ) bulges.insert( resid );
	}
	return bulges;
}


/// @brief reset begin1_, end1_, begin2_, end2_ based on ssinfo
/// @detailed abego is used for determining proper pairings in bulges
void
StrandPairing::redefine_begin_end( SS_Info2_COP const ss_info, utility::vector1< String > const & abego )
{
	if ( rgstr_shift_ == 99 ) {
		TR << "Skipping determining residue pairs, as all register shifts (99) are given as valid." << std::endl;
		return;
	}
	TR << "strand1=" << s1_ << "(" << ss_info->strand( s1_ )->begin() << "," << ss_info->strand( s1_ )->end() << ") " << std::endl;
	TR << "strand2=" << s2_ << "(" << ss_info->strand( s2_ )->begin() << "," << ss_info->strand( s2_ )->end() << ") " << std::endl;
	TR << "abego=" <<  abego<< std::endl;
	bulges1_ = compute_bulges( ss_info->strand( s1_ )->begin(), ss_info->strand( s1_ )->end(), abego );
	bulges2_ = compute_bulges( ss_info->strand( s2_ )->begin(), ss_info->strand( s2_ )->end(), abego );
	TR << "Bulges1: " << bulges1_ << " Bulges2: " << bulges2_ << std::endl;

	core::Size const s1_begin = ss_info->strand( s1_ )->begin();
	core::Size const s1_len = ss_info->strand( s1_ )->length() - bulges1_.size();
	core::Size const s2_begin = ss_info->strand( s2_ )->begin();
	core::Size const s2_end = ss_info->strand( s2_ )->end();
	core::Size const s2_len = ss_info->strand( s2_ )->length() - bulges2_.size();

	core::Size s1_pair_start, s2_pair_start, len;
	int inc = 1;

	if ( is_parallel() ) { // parallel
		if ( rgstr_shift_ >= 0 ) {
			s1_pair_start = s1_begin + rgstr_shift_;
			s2_pair_start = s2_begin;

			if ( s1_len >= (s2_len+rgstr_shift_) ) {
				//  i =========>
				//  j   =====>
				len = s2_len;
			} else {
				//  i =========>
				//  j   ==========>
				len = s1_len - rgstr_shift_;
			}
		} else {
			s1_pair_start = s1_begin;
			s2_pair_start = s2_begin - rgstr_shift_;
			if ( s1_len >= (s2_len+rgstr_shift_) ) {
				//  i      ==========>
				//  j   ==========>
				len = s2_len + rgstr_shift_;
			} else {
				//  i      =====>
				//  j   ==========>
				len = s1_len;
			}
		}

	} else { // anti parallel

		inc = -1;

		if ( rgstr_shift_ >= 0 ) {
			s1_pair_start = s1_begin + rgstr_shift_;
			s2_pair_start = s2_end;
			if ( s1_len >= (s2_len+rgstr_shift_) ) {
				//  i   ==========>
				//  j     <=====
				len = s2_len;
			} else {
				//  i   =========>
				//  j      <=========
				len = s1_len - rgstr_shift_;
			}

		} else {

			s1_pair_start = s1_begin;
			s2_pair_start = s2_end + rgstr_shift_;
			if ( s1_len >= (s2_len+rgstr_shift_) ) {
				//  i     =========>
				//  j  <=========
				len = s2_len + rgstr_shift_;
			} else {
				//  i     =========>
				//  j  <==============
				len = s1_len;
			}

		} // if( rgstr_shift_ >= 0 )
	} // if is_parallel ?


	for ( core::Size i=1; i<=len; i++ ) {
		TR.Debug << "elongating to include " << s1_pair_start << " " << s2_pair_start << std::endl;
		// if this residue is a bulge on strand 1, skip it
		if ( bulges1_.find( s1_pair_start ) != bulges1_.end() ) {
			++s1_pair_start;
			--i;
			continue;
		}
		// if this residue is a bulge on strand 2, skip it
		if ( bulges2_.find( s2_pair_start ) != bulges2_.end() ) {
			s2_pair_start += inc;
			--i;
			continue;
		}
		runtime_assert( s1_pair_start > 0 && s2_pair_start > 0 );
		if ( ! elongate( s1_pair_start, s2_pair_start, 0, 0 ) ) {
			TR.Fatal << "elongation failed ! " << std::endl;
			utility_exit_with_message("Elongation failed!");
		}

		++s1_pair_start;
		s2_pair_start += inc;
	}
	TR.Debug << "Done elongating. Size1=" << size1() << " Size2=" << size2() << std::endl;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief default constructor
StrandPairingSet::StrandPairingSet():
	spairset_name_( "" ),
	num_strands_( 0 ),
	finalized_( false ),
	empty_( utility::pointer::make_shared< StrandPairing >() )
{}


/// @brief value constructor
StrandPairingSet::StrandPairingSet( StrandPairings const & strand_pairings ):
	strand_pairings_( strand_pairings ),
	spairset_name_( "" ),
	num_strands_( 0 ),
	finalized_( false ),
	empty_( utility::pointer::make_shared< StrandPairing >() )
{
	finalize();
}


/// @brief value constructor
StrandPairingSet::StrandPairingSet( String const & spairstring, SS_Info2_COP const ssinfo ):
	spairset_name_( "" ),
	num_strands_( 0 ),
	finalized_( false ),
	empty_( utility::pointer::make_shared< StrandPairing >() )
{
	if ( spairstring == "" ) {
		return;
	}

	initialize_by_sspair_string( spairstring, ssinfo );
	finalize();
}

/// @brief value constructor
StrandPairingSet::StrandPairingSet( String const & spairstring, SS_Info2_COP const ssinfo, utility::vector1< String > const & abego ):
	spairset_name_( "" ),
	num_strands_( 0 ),
	finalized_( false ),
	empty_( utility::pointer::make_shared< StrandPairing >() )
{
	if ( spairstring == "" ) return;
	initialize_by_sspair_string_and_abego( spairstring, ssinfo, abego );
	finalize();
}

/// @brief value constructor
StrandPairingSet::StrandPairingSet( SS_Info2 const & ssinfo, DimerPairings const & dimer_pairs ):
	spairset_name_( "" ),
	num_strands_( 0 ),
	finalized_( false ),
	empty_( utility::pointer::make_shared< StrandPairing >() )
{
	initialize_by_dimer_pairs( ssinfo, dimer_pairs );
}

/// @brief destructor
StrandPairingSet::~StrandPairingSet()= default;

/// @brief clone this object
StrandPairingSetOP
StrandPairingSet::clone() const
{
	return utility::pointer::make_shared< StrandPairingSet >( *this );
}


/// @brief output detail of strand pairing set
std::ostream & operator<<( std::ostream & out, const StrandPairingSet &s )
{
	out << "#### StrandPairingSet Info " << std::endl;
	out << "# " << s.name() << std::endl;

	StrandPairings const & spairs( s.strand_pairings() );
	for ( auto const & spair : spairs ) {
		out << "# " << (*spair) << std::endl;
	}

	return out;
}


/// @brief add StrandPairingOP to StandPairingSet
void
StrandPairingSet::push_back( StrandPairingOP const sop )
{
	finalized_ = false;
	strand_pairings_.push_back( sop );
}


/// @brief add StrandPairingOP to StandPairingSet
void
StrandPairingSet::push_back_and_finalize( StrandPairingOP const sop )
{
	strand_pairings_.push_back( sop );
	finalize();
}


/// @brief clear data of this StrandParingSet
void
StrandPairingSet::clear()
{
	strand_pairings_.clear();
	spairset_name_  = "";
	finalized_ = false;
	map_strand_pairings_.clear();
	neighbor_strands_.clear();
}

/// @brief return num of strand_pairing
core::Size
StrandPairingSet::size() const
{
	return strand_pairings_.size();
}


/// @brief return spairset_name_
String
StrandPairingSet::name() const
{
	runtime_assert( finalized_ );
	return spairset_name_;
}


/// @brief return all strand pairings
StrandPairings const &
StrandPairingSet::strand_pairings() const
{
	return strand_pairings_;
}


/// @brief return one of the stand_pairings given a number
StrandPairingOP
StrandPairingSet::strand_pairing( core::Size const s ) const
{
	runtime_assert( s <= strand_pairings_.size() );
	return strand_pairings_[ s ];
}


/// @brief
StrandPairingOP
StrandPairingSet::strand_pairing( core::Size const s1, core::Size const s2 ) const
{
	runtime_assert( finalized_ );
	if ( s1 <= num_strands_ && s2 <= num_strands_ ) {
		return map_strand_pairings_[ s1 ][ s2 ];
	} else {
		return empty_;
	}
}

/// @brief
StrandPairingSet::VecSize const &
StrandPairingSet::neighbor_strands( core::Size const s ) const
{
	runtime_assert( finalized_ );
	runtime_assert( s <= num_strands_ );
	return neighbor_strands_[ s ];
}


/// @brief the name of StrandPairingSet without register shift
/// For example, 2kl8 of ferredoxin-like fold is described as 1-3.A;2-3.A;1-4.A
StrandPairingSet::String
StrandPairingSet::name_wo_rgstr() const
{
	String spairs = "";
	//std::ostringstream name;
	for ( auto const & strand_pairing : strand_pairings_ ) {
		StrandPairing const & spair( *strand_pairing );
		utility::vector1< String > sp( utility::string_split( spair.name(), '.' ) );
		runtime_assert( sp.size() == 3 );
		if ( spairs == "" ) {
			spairs = sp[ 1 ] + '.' + sp[ 2 ];
		} else {
			spairs = spairs + ';' + sp[ 1 ] + '.' + sp[ 2 ];
		}
	}
	return spairs;
}


/// @brief
bool pointer_sorter( StrandPairingCOP const a, StrandPairingCOP const b )
{
	return ( a->name() < b->name() );
}


/// @brief
void
StrandPairingSet::finalize()
{
	using VecSize = utility::vector1<core::Size>;

	finalized_ = true;

	// sort strand_parings_ by name_ of StrandPairingOP
	std::sort( strand_pairings_.begin(), strand_pairings_.end(), pointer_sorter );


	for ( StrandPairings::const_iterator it=strand_pairings_.begin(), ite=strand_pairings_.end(); it != ite; ++it ) {

		StrandPairing spair(**it);

		// set spairset_name_
		if ( spairset_name_ == "" ) {
			spairset_name_ = spair.name();
		} else {
			spairset_name_ += ';' + spair.name();
		}
		// find max number of strands
		if ( spair.s2() > num_strands_ ) {
			num_strands_ = spair.s2();
		}
	}

	map_strand_pairings_.resize( num_strands_ );
	for ( core::Size i=1; i<=num_strands_; i++ ) {

		// initialize neighbor_strands
		VecSize vec;
		neighbor_strands_.insert( std::map< core::Size, VecSize >::value_type( i, vec ) );

		// initialize map_strand_pairings_
		map_strand_pairings_[i].resize( num_strands_ );
		for ( core::Size j=1; j<=num_strands_; j++ ) {
			map_strand_pairings_[i][j] = empty_;
		}
	}

	for ( StrandPairings::const_iterator it=strand_pairings_.begin(), ite=strand_pairings_.end(); it != ite; ++it ) {
		StrandPairingOP const spair( *it );
		// define map_strand_pairings_
		map_strand_pairings_[ spair->s1() ][ spair->s2() ] = spair;
		map_strand_pairings_[ spair->s2() ][ spair->s1() ] = spair;
		// define neighbor_strands_
		VecSize & neighbor1 ( neighbor_strands_[ spair->s1() ] );
		VecSize & neighbor2 ( neighbor_strands_[ spair->s2() ] );
		neighbor1.push_back( spair->s2() );
		neighbor2.push_back( spair->s1() );
	}
	//TR << "#strands = " << num_strands_ << std::endl;
}


/// @brief
bool sort_by_length(
	protocols::fldsgn::topology::StrandPairingCOP const a,
	protocols::fldsgn::topology::StrandPairingCOP const b )
{
	return ( a->size1() >  b->size1() );
}

/// @brief
void
StrandPairingSet::drop_strand_pairs( StrandPairings const & drop_spairs )
{
	runtime_assert( drop_spairs.size() <= strand_pairings_.size() );

	StrandPairings new_spairs;
	for ( core::Size jj=1; jj<=strand_pairings_.size(); jj++ ) {
		bool drop( false );
		for ( core::Size ii=1; ii<=drop_spairs.size(); ii++ ) {
			if ( strand_pairings_[ jj ]->name() == drop_spairs[ ii ]->name() ) {
				drop = true;
				break;
			}
		}

		if ( ! drop ) {
			new_spairs.push_back( strand_pairings_[ jj ] );
		}
	}

	clear();
	strand_pairings_ = new_spairs;

	finalize();
}


/// @brief
void
StrandPairingSet::make_strand_neighbor_two()
{
	if ( ! finalized_ ) {
		finalize();
	}

	utility::vector1< utility::vector1< bool > > pairmap( num_strands_, utility::vector1< bool >( num_strands_, false ) );

	bool modified( false );
	StrandPairings drop_spairs;
	for ( core::Size ist=1; ist<=num_strands_; ist++ ) {

		StrandPairings spairs;
		if ( neighbor_strands( ist ).size() > 2 ) {

			modified = true;
			for ( auto
					it=neighbor_strands( ist ).begin(), ite=neighbor_strands( ist ).end(); it != ite; ++it ) {
				core::Size jst( *it );
				spairs.push_back( strand_pairing( ist, jst ) );
			}
			std::sort( spairs.begin(), spairs.end(), sort_by_length );

			for ( core::Size i=3; i<=spairs.size(); i++ ) {
				if ( pairmap[ spairs[ i ]->s1() ][ spairs[ i ]->s2() ] ) continue;
				pairmap[ spairs[ i ]->s1() ][ spairs[ i ]->s2() ] = true;
				pairmap[ spairs[ i ]->s2() ][ spairs[ i ]->s1() ] = true;
				drop_spairs.push_back( spairs[ i ] );
			}
		}

	}

	if ( modified ) {
		drop_strand_pairs( drop_spairs );
	}

}

void
StrandPairingSet::initialize_by_sspair_string( String const & spairstring, SS_Info2_COP const ssinfo )
{
	initialize_by_sspair_string_and_abego( spairstring, ssinfo, utility::vector1< String >() );
}

void
StrandPairingSet::initialize_by_sspair_string_and_abego( String const & spairstring, SS_Info2_COP const ssinfo, utility::vector1< String > const & abego )
{
	utility::vector1< String > spairs( utility::string_split( spairstring, ';' ) );
	for ( utility::vector1< String >::const_iterator iter = spairs.begin(); iter != spairs.end() ; ++iter ) {
		String spair( *iter );
		StrandPairingOP sp( new StrandPairing( spair ) );
		if ( ssinfo ) {
			sp->redefine_begin_end( ssinfo, abego );
		}
		push_back( sp );
	}
}


/// @brief initialize StrandPairingSet based on dimer_pairs ( under developed )
void
StrandPairingSet::initialize_by_dimer_pairs( SS_Info2 const & ssinfo, DimerPairings const & dimer_pairs )
{
	// set number of strands
	num_strands_ = ssinfo.strands().size();

	// intialize map strand pairings
	map_strand_pairings_.resize( num_strands_ );
	for ( core::Size i=1; i<=num_strands_; i++ ) {
		map_strand_pairings_[i].resize( num_strands_ );
		for ( core::Size j=1; j<=num_strands_; j++ ) {
			map_strand_pairings_[i][j] = empty_;
		}
	}

	for ( auto const & dimer_pair : dimer_pairs ) {

		DimerPairing const & dp ( *dimer_pair );

		if ( (dp.sign1() == 1 && dp.sign2() == 1) || (dp.sign1() == 2 && dp.sign2() == 2) ) continue;

		core::Size iaa = dp.res1();
		core::Size jaa = dp.res2();
		core::Size istrand = ssinfo.strand_id( iaa );
		core::Size jstrand = ssinfo.strand_id( jaa );
		core::Size ist_begin = ssinfo.strand( istrand )->begin();
		core::Size jst_begin = ssinfo.strand( jstrand )->begin();
		core::Size jst_length = ssinfo.strand( jstrand )->length();

		StrandPairingOP & spop = map_strand_pairings_[ istrand ][ jstrand ];
		if ( spop == nullptr ) {
			spop = utility::pointer::make_shared< StrandPairing >( istrand, jstrand, 0, dp.orient() );
			strand_pairings_.push_back( spop );
		}

		if ( dp.orient() == 'A' ) {

			Real rgstr = Real(iaa) - Real(ist_begin) - (Real(jst_length) - (Real(jaa) - Real(jst_begin)));
			spop->add_pair( iaa, jaa+1, dp.orient(), rgstr );
			spop->add_pair( iaa+1, jaa, dp.orient(), rgstr );

		} else if ( dp.orient() == 'P' ) {

			Real rgstr = Real(iaa) - Real(ist_begin) - ( Real(jaa) - Real(jst_begin) );
			spop->add_pair( iaa, jaa, dp.orient(), rgstr );
			spop->add_pair( iaa+1, jaa+1, dp.orient(), rgstr );

		} else {
			TR.Fatal << "Cannot understand DimerPairing orient value of '" << dp.orient() << "'" << std::endl;
			utility_exit_with_message("Invalide DimerPairings orient value.");
		}

	}  // for ( DimerPairings )

	finalize();

}


} // namespace topology
} // namespace fldsgn
} // namespace protocols
