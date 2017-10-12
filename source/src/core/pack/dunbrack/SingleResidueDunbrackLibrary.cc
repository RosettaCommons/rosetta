// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Unit headers
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>

// Package headers
///#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.tmpl.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>

#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.tmpl.hh>
#include <core/pack/dunbrack/SemiRotamericSingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/SemiRotamericSingleResidueDunbrackLibrary.tmpl.hh>

// Project headers
#include <core/conformation/Residue.hh>

#include <utility/graph/Graph.hh>

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
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/random/random.hh>


namespace core {
namespace pack {
namespace dunbrack {

static THREAD_LOCAL basic::Tracer TR( "core.pack.dunbrack" );

Real const SingleResidueDunbrackLibrary::NEUTRAL_PHI = -90; // R++ value.  Roland Dunbrack suggests -60.
Real const SingleResidueDunbrackLibrary::NEUTRAL_PSI = 130; // R++ value.  Roland Dunbrack suggests  60.

core::Real const SingleResidueDunbrackLibrary::ANGLE_DELTA( 1e-6 );
core::Real const SingleResidueDunbrackLibrary::PROB_DELTA( 1e-6 );
core::Real const SingleResidueDunbrackLibrary::ENERGY_DELTA( 1e-6 );
core::Real const SingleResidueDunbrackLibrary::COEF_DELTA( 1e-6 );

SingleResidueDunbrackLibrary::SingleResidueDunbrackLibrary(
	AA const aa,
	Size const n_rotameric_chi,
	bool dun02,
	bool use_bicubic,
	bool dun_entropy_correction,
	core::Real prob_buried,
	core::Real prob_nonburied
) :
	dun02_( dun02 ),
	use_bicubic_( use_bicubic ),
	dun_entropy_correction_( dun_entropy_correction ),
	aa_( aa ),
	n_rotameric_chi_( n_rotameric_chi ),
	n_chi_bins_( n_rotameric_chi, 0 ),
	n_chi_products_( n_rotameric_chi, 1),
	n_packed_rots_( 0 ),
	n_possible_rots_( 0 ),
	prob_to_accumulate_buried_( prob_buried ),
	prob_to_accumulate_nonburied_( prob_nonburied ),
	packed_rotno_conversion_data_current_( false )
{
	// this builds on the hard coded hack bellow
	// since NCAAs are aa_unk we cannot hard code the information
	// alternativly it is added to the residue type paramater files
	if ( aa_ != chemical::aa_unk ) {
		n_rotameric_bins_for_aa( aa_, n_chi_bins_, dun02 );
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
SingleResidueDunbrackLibrary::set_n_chi_bins( utility::vector1< Size > const & n_chi_bins )
{
	// load rot size vector
	n_chi_bins_.resize( n_chi_bins.size() );
	for ( Size i = 1; i <= n_chi_bins.size(); ++i ) n_chi_bins_[i] = n_chi_bins[i];

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

/// @details the number of wells defined for the 08 library; includes
/// the number of wells for the semi-rotameric chi (arbitrarily chosen).
void
SingleResidueDunbrackLibrary::n_rotamer_bins_for_aa(
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
SingleResidueDunbrackLibrary::n_rotameric_bins_for_aa(
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
SingleResidueDunbrackLibrary::n_rotamer_bins_for_aa_02(
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


/// @brief The base class needs to be informed about which rotamer wells
/// exist in order to create the rotwell to packed rot conversion data.
/// set_chi_nbins must be called first.
void
SingleResidueDunbrackLibrary::mark_rotwell_exists(
	utility::vector1< Size > const & rotwell
)
{
	debug_assert( ! packed_rotno_conversion_data_current_ );
	Size const rotno( rotwell_2_rotno( rotwell ) );
	rotwell_exists_[ rotno ] = true;
}

/// @brief After the derived class has marked all the rotwells that do exist,
/// the base class will create the rotwell to packedrot conversion data.
void
SingleResidueDunbrackLibrary::declare_all_existing_rotwells_encountered()
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
SingleResidueDunbrackLibrary::rotwell_2_rotno(
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
SingleResidueDunbrackLibrary::rotwell_2_rotno( Size4 const & rotwell ) const
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
SingleResidueDunbrackLibrary::rotno_2_packed_rotno( Size const rotno ) const
{
	debug_assert( packed_rotno_conversion_data_current_ );
	if ( rotno < 1 || rotno > rotno_2_packed_rotno_.size() ) return 0;
	return rotno_2_packed_rotno_[ rotno ];
}


Size
SingleResidueDunbrackLibrary::rotwell_2_packed_rotno(
	utility::vector1< Size > const & rotwell
) const
{
	debug_assert( packed_rotno_conversion_data_current_ );
	return rotno_2_packed_rotno( rotwell_2_rotno( rotwell ) );
}

/// @brief Convert from the rotamer bin indices for each chi to the
/// compacted "packed rotamer number." Returns 0 if rotwell has no corresponding packed rotno
Size
SingleResidueDunbrackLibrary::rotwell_2_packed_rotno( Size4 const & rotwell ) const
{
	debug_assert( packed_rotno_conversion_data_current_ );
	return rotno_2_packed_rotno( rotwell_2_rotno( rotwell ) );
}


void
SingleResidueDunbrackLibrary::packed_rotno_2_rotwell(
	Size const packed_rotno,
	utility::vector1< Size > & rotwell
) const
{
	rotwell = packed_rotno_2_rotwell_[ packed_rotno ];
}

void
SingleResidueDunbrackLibrary::packed_rotno_2_rotwell(
	Size const packed_rotno,
	Size4 & rotwell
) const
{
	debug_assert( packed_rotno_2_rotwell_[ packed_rotno ].size() <= rotwell.size() );
	std::copy( packed_rotno_2_rotwell_[ packed_rotno ].begin(), packed_rotno_2_rotwell_[ packed_rotno ].end(), rotwell.begin() );
}

utility::vector1< Size > const &
SingleResidueDunbrackLibrary::packed_rotno_2_rotwell( Size const packed_rotno ) const
{
	return packed_rotno_2_rotwell_[ packed_rotno ];
}


void
SingleResidueDunbrackLibrary::write_to_binary( utility::io::ozstream & out ) const
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
SingleResidueDunbrackLibrary::read_from_binary( utility::io::izstream & in )
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

/// @brief Comparison operator, mainly intended to use in ASCII/binary comparsion tests
/// Values tested should parallel those used in the read_from_binary() function.
bool
SingleResidueDunbrackLibrary::operator ==( SingleResidueRotamerLibrary const & rhs) const {
	// Raw pointer okay, we're just using it to check for conversion
	SingleResidueDunbrackLibrary const * ptr( dynamic_cast< SingleResidueDunbrackLibrary const * > ( &rhs ) );
	if ( ptr == 0 ) {
		TR << "In comparison operator: right-hand side is not a SingleResidueDunbrackLibrary." << std::endl;
		return false;
	}

	SingleResidueDunbrackLibrary const & other( dynamic_cast< SingleResidueDunbrackLibrary const & > ( rhs ) );

	// We don't call the parent == operator, as there's no meaningful data there to compare.
	bool equal( true );

	if ( aa_ != other.aa_ ) {
		TR.Debug << "Comparison failure: " << core::chemical::name_from_aa(aa_)
			<< " does not match " << core::chemical::name_from_aa(other.aa_) << std::endl;
		return false; // Major data mismatch - don't bother reporting others.
	}
	if ( dun02_ != other.dun02_ ) {
		TR.Debug << "Comparison failure: Dun02 setting mismatch: " << dun02_ << " vs. " << other.dun02_ << std::endl;
		return false; // Major data mismatch - don't bother reporting others.
	}

	/// Non-binary loaded data - check for equality for sanity purposes
	if ( n_rotameric_chi_ != other.n_rotameric_chi_ ) {
		TR.Debug << "Comparison failure: n_rotameric_chi: " << n_rotameric_chi_ << " vs. " << other.n_rotameric_chi_ << std::endl;
		equal = false;
	}
	if ( n_possible_rots_ != other.n_possible_rots_ ) {
		TR.Debug << "Comparison failure: n_possible_rots: " << n_possible_rots_ << " vs. " << other.n_possible_rots_ << std::endl;
		equal = false;
	}
	if ( packed_rotno_conversion_data_current_ != other.packed_rotno_conversion_data_current_ ) {
		TR.Debug << "Comparison failure: packed_rotno_conversion_data_current: " << packed_rotno_conversion_data_current_ << " vs. " << other.packed_rotno_conversion_data_current_ << std::endl;
		equal = false;
	}
	if ( ! numeric::equal_by_epsilon( prob_to_accumulate_buried_, other.prob_to_accumulate_buried_, PROB_DELTA) ) {
		TR.Debug << "Comparison failure: prob_to_accumulate_buried: " << prob_to_accumulate_buried_ << " vs. " << other.prob_to_accumulate_buried_ << std::endl;
		equal = false;
	}
	if ( ! numeric::equal_by_epsilon( prob_to_accumulate_nonburied_, other.prob_to_accumulate_nonburied_, PROB_DELTA) ) {
		TR.Debug << "Comparison failure: prob_to_accumulate_nonburied: " << prob_to_accumulate_nonburied_ << " vs. " << other.prob_to_accumulate_nonburied_ << std::endl;
		equal = false;
	}

	if ( n_chi_bins_.size() != other.n_chi_bins_.size() ) {
		TR.Debug << "Comparison failure: n_chi_bins vector length: " << n_chi_bins_.size() << " vs. " << other.n_chi_bins_.size() << std::endl;
		equal = false;
	} else {
		for ( core::Size ii(1); ii <= n_chi_bins_.size(); ++ii ) {
			if ( n_chi_bins_[ii] != other.n_chi_bins_[ii] ) {
				TR.Debug << "Comparison failure: n_chi_bins: " << n_chi_bins_[ii] << " vs. " << other.n_chi_bins_[ii] << std::endl;
				equal = false;
			}
		}
	}
	if ( n_chi_products_.size() != other.n_chi_products_.size() ) {
		TR.Debug << "Comparison failure: n_chi_products vector length: " << n_chi_products_.size() << " vs. " << other.n_chi_products_.size() << std::endl;
		equal = false;
	} else {
		for ( core::Size ii(1); ii <= n_chi_products_.size(); ++ii ) {
			if ( n_chi_products_[ii] != other.n_chi_products_[ii] ) {
				TR.Debug << "Comparison failure: n_chi_products: " << n_chi_products_[ii] << " vs. " << other.n_chi_products_[ii] << std::endl;
				equal = false;
			}
		}
	}
	// // For some reason, rotwell_exists_ is not matching on direct re-load comparison (empty vector on ASCII reload ?)
	//if( rotwell_exists_.size() != other.rotwell_exists_.size() ) {
	//  TR.Debug << "Comparison failure: rotwell_exists vector length: " << rotwell_exists_.size() << " vs. " << other.rotwell_exists_.size() << std::endl;
	//  equal = false;
	//} else {
	// for( core::Size ii(1); ii <= rotwell_exists_.size(); ++ii ) {
	//  if( rotwell_exists_[ii] != other.rotwell_exists_[ii] ) {
	//   TR.Debug << "Comparison failure: rotwell_exists: " << rotwell_exists_[ii] << " vs. " << other.rotwell_exists_[ii] << std::endl;
	//   equal = false;
	//  }
	// }
	//}

	/// Binary-loaded data.

	/// 1. n_packed_rots_
	if ( n_packed_rots_ != other.n_packed_rots_ ) {
		TR.Debug << "Comparsion failure in " << core::chemical::name_from_aa(aa_ )
			<< ": n_packed_rots is unequal. " << n_packed_rots_ << " vs. " << other.n_packed_rots_ << std::endl;
		equal = false;
	}

	/// 2. rotno_2_packed_rotno_
	debug_assert( n_possible_rots_ == other.n_possible_rots_ ); // This doesn't vary?
	for ( Size ii = 1; ii <= n_possible_rots_; ++ii ) {
		if ( rotno_2_packed_rotno_[ ii ] != other.rotno_2_packed_rotno_[ ii ] ) {
			TR.Debug << "Comparsion failure in " << core::chemical::name_from_aa(aa_ )
				<< ": rotno_2_packed_rotno " << ii << " - "
				<< rotno_2_packed_rotno_[ ii ] << " vs. " << other.rotno_2_packed_rotno_[ ii ] << std::endl;
			equal = false;
		}
	}

	/// 3. packed_rotno_2_rotno_
	if ( n_packed_rots_ == other.n_packed_rots_ ) {
		for ( Size ii = 1; ii <= n_packed_rots_; ++ii ) {
			if ( packed_rotno_2_rotno_[ ii ] != other.packed_rotno_2_rotno_[ ii ] ) {
				TR.Debug << "Comparsion failure in " << core::chemical::name_from_aa(aa_ )
					<< ": packed_rotno_2_rotno " << ii << " - "
					<< packed_rotno_2_rotno_[ ii ] << " vs. " << other.packed_rotno_2_rotno_[ ii ] << std::endl;
				equal = false;
			}
		}
	}

	/// 4. packed_rotno_2_rotwell_
	debug_assert( n_rotameric_chi_ == other.n_rotameric_chi_ );
	if ( n_packed_rots_ == other.n_packed_rots_ ) {
		for ( Size ii = 1; ii <= n_packed_rots_; ++ii ) {
			for ( Size jj = 1; jj <= n_rotameric_chi_; ++jj ) {
				if ( packed_rotno_2_rotwell_[ ii ][ jj ] != other.packed_rotno_2_rotwell_[ ii ][ jj ] ) {
					TR.Debug << "Comparsion failure in " << core::chemical::name_from_aa(aa_)
						<< ": packed_rotno_2_rotwell " << ii << " " << jj << " - "
						<< packed_rotno_2_rotwell_[ ii ][ jj ] << " vs. " << other.packed_rotno_2_rotwell_[ ii ][ jj ] << std::endl;
					equal = false;
				}
			}
		}
	}

	return equal;
}

/// @details not as fast as going through the packed_rotno_2_rotwell_
/// lookup table, but does the modulo converstion from a 1-based index
/// to a lexicographical index ordering.
///
/// if there are 3 chi, and 3 rotamer bins per chi, then
/// 21 would represent (3-1) * 3**2 + (1-1) * 3**1 + (3-1) * 3**0  + 1 = [ 3, 1, 3 ];
void
SingleResidueDunbrackLibrary::rotno_2_rotwell(
	Size const rotno,
	utility::vector1< Size > & rotwell
) const
{
	Size remainder = rotno - 1;
	rotwell.resize( n_rotameric_chi_ );
	std::fill( rotwell.begin(), rotwell.end(), 1 );
	for ( Size ii = 1; ii <= n_rotameric_chi_; /* no increment */ ) {
		if ( remainder >= n_chi_products_[ ii ] ) {
			remainder -= n_chi_products_[ ii ];
			++rotwell[ ii ];
		} else {
			++ii;
		}
	}
}

Real
SingleResidueDunbrackLibrary::probability_to_accumulate_while_building_rotamers(
	bool buried
) const
{
	return ( buried ? prob_to_accumulate_buried_ : prob_to_accumulate_nonburied_ );
}

/// @brief setters for accumulation probability cutoff (to support externally-controlled option dependence)
void
SingleResidueDunbrackLibrary::prob_to_accumulate( Real buried, Real nonburied )
{
	prob_to_accumulate_buried( buried );
	prob_to_accumulate_nonburied( nonburied );
}
void
SingleResidueDunbrackLibrary::prob_to_accumulate_buried( Real buried )
{
	if ( buried <= 0. || buried > 1.0 ) utility_exit_with_message("illegal probability");
	prob_to_accumulate_buried_ = buried;
}
void
SingleResidueDunbrackLibrary::prob_to_accumulate_nonburied( Real nonburied )
{
	if ( nonburied <= 0. || nonburied > 1.0 ) utility_exit_with_message("illegal probability");
	prob_to_accumulate_nonburied_ = nonburied;
}

Size SingleResidueDunbrackLibrary::memory_usage_in_bytes() const
{
	return memory_usage_static() + memory_usage_dynamic();
}

Size SingleResidueDunbrackLibrary::memory_usage_dynamic() const
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
SingleResidueDunbrackLibrary::hokey_template_workaround()
{
	utility_exit_with_message(
		"ERROR: SingleResidueDunbrackLibrary::hokey_template_workaround should never be called!");

	// Don't look at this #define--read the comment below it.

	#define INIT( CHI, BB ) \
RotamericSingleResidueDunbrackLibrary< CHI, BB > rsrdl_ ## CHI ## _ ## BB( chemical::aa_ala, false, true, true, 1.0, 1.0 ); \
SemiRotamericSingleResidueDunbrackLibrary< CHI, BB > srsrdl_ ## CHI ## _ ## BB( chemical::aa_ala, true, true, false, true, true, 1.0, 1.0 ); \
PackedDunbrackRotamer< CHI, BB, Real > prot_ ## CHI ## _ ## BB; \
rsrdl_ ## CHI ## _ ## BB.nchi(); \
srsrdl_ ## CHI ## _ ## BB.nchi(); \
rsrdl_ ## CHI ## _ ## BB.n_rotamer_bins(); \
srsrdl_ ## CHI ## _ ## BB.n_rotamer_bins(); \
rsrdl_ ## CHI ## _ ## BB.rotamer_energy( rsd, scratch ); \
srsrdl_ ## CHI ## _ ## BB.rotamer_energy( rsd, scratch ); \
rsrdl_ ## CHI ## _ ## BB.rotamer_energy_deriv( rsd, scratch ); \
srsrdl_ ## CHI ## _ ## BB.rotamer_energy_deriv( rsd, scratch ); \
rsrdl_ ## CHI ## _ ## BB.best_rotamer_energy( rsd, true, scratch ); \
srsrdl_ ## CHI ## _ ## BB.best_rotamer_energy( rsd, true, scratch ); \
rsrdl_ ## CHI ## _ ## BB.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv ); \
srsrdl_ ## CHI ## _ ## BB.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv ); \
rsrdl_ ## CHI ## _ ## BB.write_to_file( ozs ); \
srsrdl_ ## CHI ## _ ## BB.write_to_file( ozs ); \
rsrdl_ ## CHI ## _ ## BB.write_to_binary( os ); \
srsrdl_ ## CHI ## _ ## BB.write_to_binary( os ); \
rsrdl_ ## CHI ## _ ## BB.read_from_binary( is ); \
srsrdl_ ## CHI ## _ ## BB.read_from_binary( is ); \
rsrdl_ ## CHI ## _ ## BB.memory_usage_in_bytes(); \
srsrdl_ ## CHI ## _ ## BB.memory_usage_in_bytes(); \
rsrdl_ ## CHI ## _ ## BB.get_rotamer_from_chi( chi, rot ); \
srsrdl_ ## CHI ## _ ## BB.get_rotamer_from_chi( chi, rot ); \
rsrdl_ ## CHI ## _ ## BB.find_another_representative_for_unlikely_rotamer( rsd, rotwell ); \
srsrdl_ ## CHI ## _ ## BB.find_another_representative_for_unlikely_rotamer( rsd, rotwell ); \
rsrdl_ ## CHI ## _ ## BB.interpolate_rotamers( rsd, scratch, i, prot_ ## CHI ## _ ## BB ); \
srsrdl_ ## CHI ## _ ## BB.interpolate_rotamers( rsd, scratch, i, prot_ ## CHI ## _ ## BB ); \
rsrdl_ ## CHI ## _ ## BB.memory_usage_dynamic(); \
srsrdl_ ## CHI ## _ ## BB.memory_usage_dynamic(); \
rsrdl_ ## CHI ## _ ## BB.memory_usage_static(); \
srsrdl_ ## CHI ## _ ## BB.memory_usage_static(); \
rsrdl_ ## CHI ## _ ## BB.get_all_rotamer_samples( bb_FIVE ); \
srsrdl_ ## CHI ## _ ## BB.get_all_rotamer_samples( bb_FIVE ); \
rsrdl_ ## CHI ## _ ## BB.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );  \
srsrdl_ ## CHI ## _ ## BB.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true ); \
rsrdl_ ## CHI ## _ ## BB.get_probability_for_rotamer( bb_ ## BB, 1 ); \
srsrdl_ ## CHI ## _ ## BB.get_probability_for_rotamer( bb_ ## BB, 1 ); \
rsrdl_ ## CHI ## _ ## BB.get_rotamer( bb_ ## BB, 1 ); \
srsrdl_ ## CHI ## _ ## BB.get_rotamer( bb_ ## BB, 1 ); \
srsrdl_ ## CHI ## _ ## BB.interpolate_nrchi_values( sizevec_ ## BB, sizevec_ ## BB, realvec_ ## BB, 1, realvec );

	// If you are adding support for rotamers with more than four chis or five
	// backbone dihedrals, good for you. Have you created const Sizes for the
	// additional numbers in core/pack/dunbrack/DunbrackRotamer.fwd.hh yet?
	// If you have, add the new combinations to this clever FOREACH define.

	#define FOREACH_Chi_BB(INIT) \
INIT(   ONE,   ONE ) \
INIT(   ONE,   TWO ) \
INIT(   ONE, THREE ) \
INIT(   ONE,  FOUR ) \
INIT(   ONE,  FIVE ) \
INIT(   TWO,   ONE ) \
INIT(   TWO,   TWO ) \
INIT(   TWO, THREE ) \
INIT(   TWO,  FOUR ) \
INIT(   TWO,  FIVE ) \
INIT( THREE,   ONE ) \
INIT( THREE,   TWO ) \
INIT( THREE, THREE ) \
INIT( THREE,  FOUR ) \
INIT( THREE,  FIVE ) \
INIT(  FOUR,   ONE ) \
INIT(  FOUR,   TWO ) \
INIT(  FOUR, THREE ) \
INIT(  FOUR,  FOUR ) \
INIT(  FOUR,  FIVE ) \
INIT(  FIVE,   ONE ) \
INIT(  FIVE,   TWO ) \
INIT(  FIVE, THREE ) \
INIT(  FIVE,  FOUR ) \
INIT(  FIVE,  FIVE )


	// If it's a new number of backbone angles you have added, add lines here.
	// And increase the size of "realvec" too.
	utility::fixedsizearray1< Real,   ONE > bb_ONE;
	utility::fixedsizearray1< Real,   TWO > bb_TWO;
	utility::fixedsizearray1< Real, THREE > bb_THREE;
	utility::fixedsizearray1< Real,  FOUR > bb_FOUR;
	utility::fixedsizearray1< Real,  FIVE > bb_FIVE;

	utility::fixedsizearray1< Size,   ONE > sizevec_ONE;
	utility::fixedsizearray1< Real,   ONE > realvec_ONE;
	utility::fixedsizearray1< Size,   TWO > sizevec_TWO;
	utility::fixedsizearray1< Real,   TWO > realvec_TWO;
	utility::fixedsizearray1< Size, THREE > sizevec_THREE;
	utility::fixedsizearray1< Real, THREE > realvec_THREE;
	utility::fixedsizearray1< Size,  FOUR > sizevec_FOUR;
	utility::fixedsizearray1< Real,  FOUR > realvec_FOUR;
	utility::fixedsizearray1< Size,  FIVE > sizevec_FIVE;
	utility::fixedsizearray1< Real,  FIVE > realvec_FIVE;
	utility::vector1< Real > realvec( FIVE );

	chemical::ResidueType rt( NULL, NULL, NULL, NULL );
	conformation::Residue rsd( rt, true );
	RotamerLibraryScratchSpace scratch;
	Size4 rotwell;
	Size i(0);

	pose::Pose pose;
	scoring::ScoreFunction sfxn;
	pack::task::PackerTask_ task( pose );

	utility::graph::GraphOP png;
	chemical::ResidueTypeCOP cr;
	utility::vector1< utility::vector1< Real > > ecs;
	rotamers::RotamerVector rv;
	ChiVector chiv;

	utility::io::ozstream os;
	utility::io::ozstream ozs;
	utility::io::izstream is;

	utility::vector1< Real > chi; utility::vector1< Size > rot;

	FOREACH_Chi_BB( INIT );
}


} // namespace dunbrack
} // namespace scoring
} // namespace core
