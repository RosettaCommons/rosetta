// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rotamers/SingleResiduePeptoidLibrary.hh
/// @brief  SingleResiduePeptoidLibrary class
/// @brief  Similar to SingleResidueDunbrackLibrary class
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// Unit headers
#include <core/pack/rotamers/SingleResiduePeptoidLibrary.hh>

// Package headers
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/rotamers/RotamericSingleResiduePeptoidLibrary.hh>
#include <core/pack/rotamers/RotamericSingleResiduePeptoidLibrary.tmpl.hh>

// Project headers
#include <core/conformation/Residue.hh>

#include <core/graph/Graph.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <core/pack/task/PackerTask_.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>

#include <basic/basic.hh>
#include <basic/Tracer.hh>

#include <core/pack/task/PackerTask.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

//Auto Headers
#include <numeric/random/random.hh>
#include <boost/cstdint.hpp>

namespace core {
namespace pack {
namespace rotamers {

static THREAD_LOCAL basic::Tracer SRPL_TR( "core.pack.rotamers.SingleResiduePeptoidLibrary" );

// static const Real MIN_ROT_PROB = 1.e-8; - KAB - unused variable

Real const SingleResiduePeptoidLibrary::NEUTRAL_OMG = 180; // DOUG DOUG DOUG: Need to determine good values for peptoids
Real const SingleResiduePeptoidLibrary::NEUTRAL_PHI = -90; // DOUG DOUG DOUG: Need to determine good values for peptoids
Real const SingleResiduePeptoidLibrary::NEUTRAL_PSI = 180; // DOUG DOUG DOUG: Need to determine good values for peptoids

SingleResiduePeptoidLibrary::SingleResiduePeptoidLibrary(
	Size const n_rotameric_chi
) :
	dun02_( true ), /// DOUG DOUG DOUG Setting this always to true for now, if this makes it into trunk someone tell me
	aa_( chemical::aa_unk ), /// DOUG DOUG DOUG Setting this always to true for now, if this makes it into trunk someone tell me
	n_rotameric_chi_( n_rotameric_chi ),
	n_chi_bins_( n_rotameric_chi, 0 ),
	n_chi_products_( n_rotameric_chi, 1),
	n_packed_rots_( 0 ),
	n_possible_rots_( 0 ),
	prob_to_accumulate_buried_(2.00), // DOUG DOUG DOUG reset to 0.98
	prob_to_accumulate_nonburied_(2.00), // DOUG DOUG DOUG reset to 0.95
	packed_rotno_conversion_data_current_( false )
{
	// this builds on the hard coded hack bellow
	// since NCAAs are aa_unk we cannot hard code the information
	// alternativly it is added to the residue type paramater files
	read_options();
	if ( aa_ != chemical::aa_unk ) {
		n_rotameric_bins_for_aa( aa_, n_chi_bins_, dun02_ );
		for ( Size ii = n_rotameric_chi_; ii > 1; --ii ) {
			n_chi_products_[ ii - 1 ] = n_chi_products_[ ii ] * n_chi_bins_[ ii ];
		}
		n_possible_rots_ = n_rotameric_chi_ == 0 ? 0 : n_chi_products_[ 1 ] * n_chi_bins_[ 1 ];
		rotwell_exists_.resize( n_possible_rots_ );
		rotno_2_packed_rotno_.resize( n_possible_rots_ );
		std::fill( rotwell_exists_.begin(), rotwell_exists_.end(), false );
		std::fill( rotno_2_packed_rotno_.begin(), rotno_2_packed_rotno_.end(), 0 );
	}
}

/// @details Sets the number of bins for a particular chi angle, used for the NCAAs, info for CAAs is hardcoded bellow
void
SingleResiduePeptoidLibrary::set_n_chi_bins( utility::vector1< Size > const & n_chi_bins )
{
	// load rot size vector
	n_chi_bins_.resize( n_chi_bins.size() );
	for ( Size i = 1; i <= n_chi_bins.size(); ++i ) {
		n_chi_bins_[i] = n_chi_bins[i];
	}

	// basically the same as the ctor
	for ( Size ii = n_rotameric_chi_; ii > 1; --ii ) {
		n_chi_products_[ ii - 1 ] = n_chi_products_[ ii ] * n_chi_bins_[ ii ];
	}
	n_possible_rots_ = n_rotameric_chi_ == 0 ? 0 : n_chi_products_[ 1 ] * n_chi_bins_[ 1 ];
	rotwell_exists_.resize( n_possible_rots_ );
	rotno_2_packed_rotno_.resize( n_possible_rots_ );
	std::fill( rotwell_exists_.begin(), rotwell_exists_.end(), false );
	std::fill( rotno_2_packed_rotno_.begin(), rotno_2_packed_rotno_.end(), 0 );
}

SingleResiduePeptoidLibrary::~SingleResiduePeptoidLibrary() {}

void
SingleResiduePeptoidLibrary::read_options()
{
	using namespace basic::options;
	if ( option[ OptionKeys::packing::dunbrack_prob_buried ].user() ) {
		prob_to_accumulate_buried( option[ OptionKeys::packing::dunbrack_prob_buried ]() );
	}
	if ( option[ OptionKeys::packing::dunbrack_prob_nonburied ].user() ) {
		prob_to_accumulate_nonburied( option[ OptionKeys::packing::dunbrack_prob_nonburied ]() );
	}
}

/// @details the number of wells defined for the 08 library; includes
/// the number of wells for the semi-rotameric chi (arbitrarily chosen).
void
SingleResiduePeptoidLibrary::n_rotamer_bins_for_aa(
	chemical::AA const aa,
	RotVector & rot
)
{
	/// HARD CODED HACK -- SHOULD MOVE THIS TO BECOME A VIRTUAL FUNCTION CALL
	using namespace chemical;
	switch ( aa ) {
	case aa_ala : rot.resize( 0 ); break;
	case aa_cys : rot.resize( 1 ); rot[ 1 ] = 3; break;
	case aa_asp : rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 6; break;
	case aa_glu : rot.resize( 3 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = 6; break;
	case aa_phe : rot.resize( 2 ); rot[ 1 ] = 3; rot[ 2 ] = 6; break;
	case aa_gly : rot.resize( 0 ); break;
	case aa_his : rot.resize( 2 ); rot[ 1 ] = 3; rot[ 2 ] = 12; break;
	case aa_ile : rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
	case aa_lys : rot.resize( 4 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = rot[ 4 ] = 3; break;
	case aa_leu : rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
	case aa_met : rot.resize( 3 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = 3; break;
	case aa_asn : rot.resize( 2 ); rot[ 1 ] = 3; rot[ 2 ] = 12; break;
	case aa_pro : rot.resize( 3 ); rot[ 1 ] = 2; rot[ 2 ] = rot[ 3 ] = 1; break;
	case aa_gln : rot.resize( 3 ); rot[ 1 ] = rot[ 2 ] = 3; rot[ 3 ] = 12; break;
	case aa_arg : rot.resize( 4 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = rot[ 4 ] = 3; break;
	case aa_ser : rot.resize( 1 ); rot[ 1 ] = 3; break;
	case aa_thr : rot.resize( 1 ); rot[ 1 ] = 3; break;
	case aa_val : rot.resize( 1 ); rot[ 1 ] = 3; break;
	case aa_trp : rot.resize( 2 ); rot[ 1 ] = 3; rot[ 2 ] = 12; break;
	case aa_tyr : rot.resize( 2 ); rot[ 1 ] = 3; rot[ 2 ] = 6; break;
	default :
		rot.resize( 0 );
	}
}

/// @details For the rotameric chi -- not all chi are rotameric
void
SingleResiduePeptoidLibrary::n_rotameric_bins_for_aa(
	chemical::AA const aa,
	RotVector & rot,
	bool dun02
)
{
	if ( dun02 ) {
		n_rotamer_bins_for_aa_02( aa, rot );
		return;
	}


	using namespace chemical;
	switch ( aa ) {
	case aa_ala : rot.resize( 0 ); break;
	case aa_cys : rot.resize( 1 ); rot[ 1 ] = 3; break;
	case aa_asp : rot.resize( 1 ); rot[ 1 ] = 3; break;
	case aa_glu : rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
	case aa_phe : rot.resize( 1 ); rot[ 1 ] = 3; break;
	case aa_gly : rot.resize( 0 ); break;
	case aa_his : rot.resize( 1 ); rot[ 1 ] = 3; break;
	case aa_ile : rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
	case aa_lys : rot.resize( 4 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = rot[ 4 ] = 3; break;
	case aa_leu : rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
	case aa_met : rot.resize( 3 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = 3; break;
	case aa_asn : rot.resize( 1 ); rot[ 1 ] = 3; break;
	case aa_pro : rot.resize( 3 ); rot[ 1 ] = 2; rot[ 2 ] = rot[ 3 ] = 1; break;
	case aa_gln : rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
	case aa_arg : rot.resize( 4 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = rot[ 4 ] = 3; break;
	case aa_ser : rot.resize( 1 ); rot[ 1 ] = 3; break;
	case aa_thr : rot.resize( 1 ); rot[ 1 ] = 3; break;
	case aa_val : rot.resize( 1 ); rot[ 1 ] = 3; break;
	case aa_trp : rot.resize( 1 ); rot[ 1 ] = 3; break;
	case aa_tyr : rot.resize( 1 ); rot[ 1 ] = 3; break;
	default :
		rot.resize( 0 );
	}
}

/// @details To continue supporting the 2002 Dunbrack library, we
/// need to preserve the rotamer well definitions that it made.
void
SingleResiduePeptoidLibrary::n_rotamer_bins_for_aa_02(
	chemical::AA const aa,
	RotVector & rot
)
{
	using namespace chemical;
	switch ( aa ) {
	case aa_ala : rot.resize( 0 ); break;
	case aa_cys : rot.resize( 1 ); rot[ 1 ] = 3; break;
	case aa_asp : rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
	case aa_glu : rot.resize( 3 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = 3; break;
	case aa_phe : rot.resize( 2 ); rot[ 1 ] = 3; rot[ 2 ] = 2; break;
	case aa_gly : rot.resize( 0 ); break;
	case aa_his : rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
	case aa_ile : rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
	case aa_lys : rot.resize( 4 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = rot[ 4 ] = 3; break;
	case aa_leu : rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
	case aa_met : rot.resize( 3 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = 3; break;
	case aa_asn : rot.resize( 2 ); rot[ 1 ] = 3; rot[ 2 ] = 6; break;
	case aa_pro : rot.resize( 3 ); rot[ 1 ] = 2; rot[ 2 ] = rot[ 3 ] = 1; break;
	case aa_gln : rot.resize( 3 ); rot[ 1 ] = rot[ 2 ] = 3; rot[ 3 ] = 4; break;
	case aa_arg : rot.resize( 4 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = rot[ 4 ] = 3; break;
	case aa_ser : rot.resize( 1 ); rot[ 1 ] = 3; break;
	case aa_thr : rot.resize( 1 ); rot[ 1 ] = 3; break;
	case aa_val : rot.resize( 1 ); rot[ 1 ] = 3; break;
	case aa_trp : rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
	case aa_tyr : rot.resize( 2 ); rot[ 1 ] = 3; rot[ 2 ] = 2; break;
	default :
		rot.resize( 0 );
	}
}

namespace peptoid {
/// DOUG DOUG DOUG Should this function be part of the class, or is it a utility function just for this file, might need to have something like this for the NCAAs
/// 2002 Library hard code symmetry information.
Real
subtract_chi_angles(
	Real chi1,
	Real chi2,
	chemical::AA const & aa,
	int chino
)
{
	using namespace chemical;
	using basic::periodic_range;


	if ( chino == 2 ) {   // handle symmetry cases for chi2
		//  -- not for HIS in new dunbrack version
		if ( aa == aa_phe || aa == aa_tyr ) {
			return periodic_range((chi1-chi2),180.);
		} else if ( aa == aa_asp ) {
			return periodic_range((chi1-chi2),180.);
		}
	} else if ( chino == 3 ) {   // handle symmetry cases for chi3
		if ( aa == aa_glu ) {
			return periodic_range((chi1-chi2),180.);
		}
	}
	return periodic_range((chi1-chi2),360.);
}

}


/// @brief The base class needs to be informed about which rotamer wells
/// exist in order to create the rotwell to packed rot conversion data.
/// set_chi_nbins must be called first.
void
SingleResiduePeptoidLibrary::mark_rotwell_exists(
	utility::vector1< Size > const & rotwell
)
{
	debug_assert( ! packed_rotno_conversion_data_current_ );
	Size const rotno( rotwell_2_rotno( rotwell ) );
	rotwell_exists_[ rotno ] = true;
}

/// @brief After the derived class has marked all the rotwells that do exist,
/// the base class will create the rotwell to packerot conversion data.
void
SingleResiduePeptoidLibrary::declare_all_existing_rotwells_encountered()
{
	Size count = 0;
	for ( Size ii = 1; ii <= n_possible_rots_; ++ii ) {
		if ( rotwell_exists_[ ii ] ) ++count;
	}
	n_packed_rots_ = count;
	packed_rotno_2_rotno_.resize( n_packed_rots_ );
	packed_rotno_2_rotwell_.resize( n_packed_rots_ );
	count = 0;
	utility::vector1< Size > rotwell;
	for ( Size ii = 1; ii <= n_possible_rots_; ++ii ) {
		if ( rotwell_exists_[ ii ] ) {
			++count;
			packed_rotno_2_rotno_[ count ] = ii;
			rotno_2_packed_rotno_[ ii    ] = count;
			rotno_2_rotwell( ii, rotwell );
			packed_rotno_2_rotwell_[ count ] = rotwell;
		} else {
			rotno_2_packed_rotno_[ ii ] = 0;
		}
	}
	packed_rotno_conversion_data_current_ = true;

	// force deallocatation of rotamer_exists_ data
	// neither resize() nor clear() guarantee emptying; swap does.
	utility::vector1< bool > empty_vector;
	rotwell_exists_.swap( empty_vector );

}

/// Conversion functions

/// @brief Convert from the rotamer bin indices for each chi to the
/// (non-compact) "rotamer number"
Size
SingleResiduePeptoidLibrary::rotwell_2_rotno(
	utility::vector1< Size > const & rotwell
) const
{
	debug_assert( n_chi_products_.size() <= rotwell.size() );
	Size runsum = 1;
	for ( Size ii = 1; ii <= n_chi_products_.size(); ++ii ) {
		runsum += n_chi_products_[ ii ]*( rotwell[ ii ] - 1 );
	}
	return runsum;
}

/// @brief Convert from the rotamer bin indices for each chi to the
/// (non-compact) "rotamer number"
Size
SingleResiduePeptoidLibrary::rotwell_2_rotno( Size4 const & rotwell ) const
{
	debug_assert( n_chi_products_.size() <= rotwell.size() );
	Size runsum = 1;
	for ( Size ii = 1; ii <= n_chi_products_.size(); ++ii ) {
		runsum += n_chi_products_[ ii ]*( rotwell[ ii ] - 1 );
	}
	return runsum;

}


/// @brief Convert from the rotamer number to the compacted
/// "packed rotamer number".  Returns 0 if rotno has no corresponding packed rotno.
Size
SingleResiduePeptoidLibrary::rotno_2_packed_rotno( Size const rotno ) const
{
	debug_assert( packed_rotno_conversion_data_current_ );
	return rotno_2_packed_rotno_[ rotno ];
}


Size
SingleResiduePeptoidLibrary::rotwell_2_packed_rotno(
	utility::vector1< Size > const & rotwell
) const
{
	debug_assert( packed_rotno_conversion_data_current_ );
	return rotno_2_packed_rotno( rotwell_2_rotno( rotwell ) );
}

/// @brief Convert from the rotamer bin indices for each chi to the
/// compacted "packed rotamer number." Returns 0 if rotwell has no corresponding packed rotno
Size
SingleResiduePeptoidLibrary::rotwell_2_packed_rotno( Size4 const & rotwell ) const
{
	debug_assert( packed_rotno_conversion_data_current_ );
	return rotno_2_packed_rotno( rotwell_2_rotno( rotwell ) );
}


void
SingleResiduePeptoidLibrary::packed_rotno_2_rotwell(
	Size const packed_rotno,
	utility::vector1< Size > & rotwell
) const
{
	rotwell = packed_rotno_2_rotwell_[ packed_rotno ];
}

void
SingleResiduePeptoidLibrary::packed_rotno_2_rotwell(
	Size const packed_rotno,
	Size4 & rotwell
) const
{
	debug_assert( packed_rotno_2_rotwell_[ packed_rotno ].size() <= rotwell.size() );
	std::copy( packed_rotno_2_rotwell_[ packed_rotno ].begin(), packed_rotno_2_rotwell_[ packed_rotno ].end(), rotwell.begin() );
}

utility::vector1< Size > const &
SingleResiduePeptoidLibrary::packed_rotno_2_rotwell( Size const packed_rotno ) const
{
	return packed_rotno_2_rotwell_[ packed_rotno ];
}


void
SingleResiduePeptoidLibrary::write_to_binary( utility::io::ozstream & out ) const
{
	using namespace boost;
	/// 1. n_packed_rots_
	{
		boost::int32_t n_packed_rots( n_packed_rots_ );
		out.write( (char*) & n_packed_rots, sizeof( boost::int32_t ));
	}


	/// 2. rotno_2_packed_rotno_
	{
		boost::int32_t * rotno_2_packed_rotno = new boost::int32_t[ n_possible_rots_ ];
		for ( Size ii = 1; ii <= n_possible_rots_; ++ii ) rotno_2_packed_rotno[ ii - 1 ] = rotno_2_packed_rotno_[ ii ];
		out.write( (char*) rotno_2_packed_rotno, n_possible_rots_ * sizeof( boost::int32_t ) );
		delete [] rotno_2_packed_rotno; rotno_2_packed_rotno = 0;
	}

	/// 3. packed_rotno_2_rotno_
	{
		boost::int32_t * packed_rotno_2_rotno = new boost::int32_t[ n_packed_rots_ ];
		for ( Size ii = 1; ii <= n_packed_rots_; ++ii ) packed_rotno_2_rotno[ ii - 1 ] = packed_rotno_2_rotno_[ ii ];
		out.write( (char*) packed_rotno_2_rotno, n_packed_rots_ * sizeof( boost::int32_t ) );
		delete [] packed_rotno_2_rotno; packed_rotno_2_rotno = 0;
	}

	/// 4. packed_rotno_2_rotwell_
	{
		boost::int32_t * packed_rotno_2_rotwell = new boost::int32_t[ n_packed_rots_ * n_rotameric_chi_ ];
		Size count( 0 );
		for ( Size ii = 1; ii <= n_packed_rots_; ++ii ) {
			for ( Size jj = 1; jj <= n_rotameric_chi_; ++jj ) {
				packed_rotno_2_rotwell[ count ] = packed_rotno_2_rotwell_[ ii ][ jj ];
				++count;
			}
		}
		out.write( (char*) packed_rotno_2_rotwell, n_packed_rots_ * n_rotameric_chi_ * sizeof( boost::int32_t ) );
		delete [] packed_rotno_2_rotwell; packed_rotno_2_rotwell = 0;
	}

}

void
SingleResiduePeptoidLibrary::read_from_binary( utility::io::izstream & in )
{

	/// 1. n_packed_rots_
	{
		boost::int32_t n_packed_rots( 0 );
		in.read( (char*) & n_packed_rots, sizeof( boost::int32_t ));
		n_packed_rots_ = n_packed_rots;
	}

	/// 2. rotno_2_packed_rotno_
	{
		boost::int32_t * rotno_2_packed_rotno = new boost::int32_t[ n_possible_rots_ ];
		in.read( (char*) rotno_2_packed_rotno, n_possible_rots_ * sizeof( boost::int32_t ) );
		for ( Size ii = 1; ii <= n_possible_rots_; ++ii ) rotno_2_packed_rotno_[ ii ] = rotno_2_packed_rotno[ ii - 1 ];
		delete [] rotno_2_packed_rotno; rotno_2_packed_rotno = 0;
	}

	/// 3. packed_rotno_2_rotno_
	{
		boost::int32_t * packed_rotno_2_rotno = new boost::int32_t[ n_packed_rots_ ];
		in.read( (char*) packed_rotno_2_rotno, n_packed_rots_ * sizeof( boost::int32_t ) );
		packed_rotno_2_rotno_.resize( n_packed_rots_ );
		for ( Size ii = 1; ii <= n_packed_rots_; ++ii ) packed_rotno_2_rotno_[ ii ] = packed_rotno_2_rotno[ ii - 1 ];
		delete [] packed_rotno_2_rotno; packed_rotno_2_rotno = 0;
	}

	/// 4. packed_rotno_2_rotwell_
	{
		boost::int32_t * packed_rotno_2_rotwell = new boost::int32_t[ n_packed_rots_ * n_rotameric_chi_ ];
		in.read( (char*) packed_rotno_2_rotwell, n_packed_rots_ * n_rotameric_chi_ * sizeof( boost::int32_t ) );
		packed_rotno_2_rotwell_.resize( n_packed_rots_ );
		Size count( 0 );
		for ( Size ii = 1; ii <= n_packed_rots_; ++ii ) {
			packed_rotno_2_rotwell_[ ii ].resize( n_rotameric_chi_ );
			for ( Size jj = 1; jj <= n_rotameric_chi_; ++jj ) {
				packed_rotno_2_rotwell_[ ii ][ jj ] = packed_rotno_2_rotwell[ count ];
				++count;
			}
		}
		delete [] packed_rotno_2_rotwell;
	}
	packed_rotno_conversion_data_current_ = true;
}

/// @details not as fast as going through the packed_rotno_2_rotwell_
/// lookup table, but does the modulo converstion from a 1-based index
/// to a lexicographical index ordering.
///
/// if there are 3 chi, and 3 rotamer bins per chi, then
/// 21 would represent (3-1) * 3**2 + (1-1) * 3**1 + (3-1) * 3**0  + 1 = [ 3, 1, 3 ];
void
SingleResiduePeptoidLibrary::rotno_2_rotwell(
	Size const rotno,
	utility::vector1< Size > & rotwell
) const
{
	Size remainder = rotno - 1;
	rotwell.resize( n_rotameric_chi_ );
	std::fill( rotwell.begin(), rotwell.end(), 1 );
	for ( Size ii = 1; ii <= n_rotameric_chi_; /* no increment */ ) {
		if ( remainder > n_chi_products_[ ii ] ) {
			remainder -= n_chi_products_[ ii ];
			++rotwell[ ii ];
		} else {
			++ii;
		}
	}
}

Real
SingleResiduePeptoidLibrary::probability_to_accumulate_while_building_rotamers(
	bool buried
) const
{
	return ( buried ? prob_to_accumulate_buried_ : prob_to_accumulate_nonburied_ );
}

/// @brief setters for accumulation probability cutoff (to support externally-controlled option dependence)
void
SingleResiduePeptoidLibrary::prob_to_accumulate( Real buried, Real nonburied )
{
	prob_to_accumulate_buried( buried );
	prob_to_accumulate_nonburied( nonburied );
}
void
SingleResiduePeptoidLibrary::prob_to_accumulate_buried( Real buried )
{
	if ( buried <= 0. || buried > 1.0 ) utility_exit_with_message("illegal probability");
	prob_to_accumulate_buried_ = buried;
}
void
SingleResiduePeptoidLibrary::prob_to_accumulate_nonburied( Real nonburied )
{
	if ( nonburied <= 0. || nonburied > 1.0 ) utility_exit_with_message("illegal probability");
	prob_to_accumulate_nonburied_ = nonburied;
}

Size SingleResiduePeptoidLibrary::memory_usage_in_bytes() const
{
	return memory_usage_static() + memory_usage_dynamic();
}

Size SingleResiduePeptoidLibrary::memory_usage_dynamic() const
{
	Size total = 0;
	total += n_chi_bins_.size() * sizeof( Size );
	total += n_chi_products_.size() * sizeof( Size );
	total += rotwell_exists_.size() * sizeof( Size );
	total += rotno_2_packed_rotno_.size() * sizeof( Size );
	total += packed_rotno_2_rotno_.size() * sizeof( Size );
	for ( Size ii = 1; ii <= packed_rotno_2_rotwell_.size(); ++ii ) {
		total += packed_rotno_2_rotwell_[ ii ].size() * sizeof( Size );
	}
	total += packed_rotno_2_rotwell_.size() * sizeof( utility::vector1< Size > );
	return total;
}


/// @details forces instantiation of virtual functions for templated
/// derived classes... part of the uglyness of mixing templates and
/// polymorphism.  Never invoke this function.
void
SingleResiduePeptoidLibrary::hokey_template_workaround()
{
	debug_assert( false );
	utility_exit_with_message(
		"ERROR: SingleResiduePeptoidLibrary::hokey_template_workaround should never be called!");

	RotamericSingleResiduePeptoidLibrary< ONE, THREE >   rsrdl_1;
	RotamericSingleResiduePeptoidLibrary< TWO, THREE >   rsrdl_2;
	RotamericSingleResiduePeptoidLibrary< THREE, THREE > rsrdl_3;
	RotamericSingleResiduePeptoidLibrary< FOUR, THREE >  rsrdl_4;

	chemical::ResidueType rt( NULL, NULL, NULL, NULL );
	conformation::Residue rsd( rt, true );
	RotamerLibraryScratchSpace scratch;

	rsrdl_1.nchi();
	rsrdl_2.nchi();
	rsrdl_3.nchi();
	rsrdl_4.nchi();

	rsrdl_1.n_rotamer_bins();
	rsrdl_2.n_rotamer_bins();
	rsrdl_3.n_rotamer_bins();
	rsrdl_4.n_rotamer_bins();

	rsrdl_1.rotamer_energy( rsd, scratch );
	rsrdl_2.rotamer_energy( rsd, scratch );
	rsrdl_3.rotamer_energy( rsd, scratch );
	rsrdl_4.rotamer_energy( rsd, scratch );

	rsrdl_1.rotamer_energy_deriv( rsd, scratch );
	rsrdl_2.rotamer_energy_deriv( rsd, scratch );
	rsrdl_3.rotamer_energy_deriv( rsd, scratch );
	rsrdl_4.rotamer_energy_deriv( rsd, scratch );

	pose::Pose pose; scoring::ScoreFunction sfxn; pack::task::PackerTask_ task( pose );
	core::graph::GraphOP png;
	chemical::ResidueTypeCOP cr;
	utility::vector1< utility::vector1< Real > > ecs;
	RotamerVector rv;

	rsrdl_1.best_rotamer_energy( rsd, true, scratch );
	rsrdl_2.best_rotamer_energy( rsd, true, scratch );
	rsrdl_3.best_rotamer_energy( rsd, true, scratch );
	rsrdl_4.best_rotamer_energy( rsd, true, scratch );

	rsrdl_1.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	rsrdl_2.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	rsrdl_3.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	rsrdl_4.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );

	utility::io::ozstream ozs;
	rsrdl_1.write_to_file( ozs );
	rsrdl_2.write_to_file( ozs );
	rsrdl_3.write_to_file( ozs );
	rsrdl_4.write_to_file( ozs );

	utility::io::ozstream os;
	rsrdl_1.write_to_binary( os );
	rsrdl_2.write_to_binary( os );
	rsrdl_3.write_to_binary( os );
	rsrdl_4.write_to_binary( os );

	utility::io::izstream is;
	rsrdl_1.read_from_binary( is );
	rsrdl_2.read_from_binary( is );
	rsrdl_3.read_from_binary( is );
	rsrdl_4.read_from_binary( is );

	rsrdl_1.memory_usage_in_bytes();
	rsrdl_2.memory_usage_in_bytes();
	rsrdl_3.memory_usage_in_bytes();
	rsrdl_4.memory_usage_in_bytes();

	utility::vector1< Real > chi; utility::vector1< Size > rot;
	rsrdl_1.get_rotamer_from_chi( chi, rot );
	rsrdl_2.get_rotamer_from_chi( chi, rot );
	rsrdl_3.get_rotamer_from_chi( chi, rot );
	rsrdl_4.get_rotamer_from_chi( chi, rot );

	rsrdl_1.memory_usage_dynamic();
	rsrdl_2.memory_usage_dynamic();
	rsrdl_3.memory_usage_dynamic();
	rsrdl_4.memory_usage_dynamic();

	rsrdl_1.memory_usage_static();
	rsrdl_2.memory_usage_static();
	rsrdl_3.memory_usage_static();
	rsrdl_4.memory_usage_static();

	rsrdl_1.get_all_rotamer_samples( 0.0, 0.0, 0.0 );
	rsrdl_2.get_all_rotamer_samples( 0.0, 0.0, 0.0 );
	rsrdl_3.get_all_rotamer_samples( 0.0, 0.0, 0.0 );
	rsrdl_4.get_all_rotamer_samples( 0.0, 0.0, 0.0 );

	ChiVector chiv;

	rsrdl_1.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	rsrdl_2.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	rsrdl_3.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	rsrdl_4.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );

	rsrdl_1.get_probability_for_rotamer( 0.0, 0.0, 0.0, 1 );
	rsrdl_2.get_probability_for_rotamer( 0.0, 0.0, 0.0, 1 );
	rsrdl_3.get_probability_for_rotamer( 0.0, 0.0, 0.0, 1 );
	rsrdl_4.get_probability_for_rotamer( 0.0, 0.0, 0.0, 1 );

	rsrdl_1.get_rotamer( 0.0, 0.0, 0.0, 1 );
	rsrdl_2.get_rotamer( 0.0, 0.0, 0.0, 1 );
	rsrdl_3.get_rotamer( 0.0, 0.0, 0.0, 1 );
	rsrdl_4.get_rotamer( 0.0, 0.0, 0.0, 1 );
}


} // namespace rotamers
} // namespace scoring
} // namespace core

