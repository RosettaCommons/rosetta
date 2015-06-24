// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/random/random.hh>


namespace core {
namespace pack {
namespace dunbrack {

static thread_local basic::Tracer TR( "core.pack.dunbrack" );

Real const SingleResidueDunbrackLibrary::NEUTRAL_PHI = -90; // R++ value.  Roland Dunbrack suggests -60.
Real const SingleResidueDunbrackLibrary::NEUTRAL_PSI = 130; // R++ value.  Roland Dunbrack suggests  60.

core::Real const SingleResidueDunbrackLibrary::ANGLE_DELTA( 1e-6 );
core::Real const SingleResidueDunbrackLibrary::PROB_DELTA( 1e-6 );
core::Real const SingleResidueDunbrackLibrary::ENERGY_DELTA( 1e-6 );
core::Real const SingleResidueDunbrackLibrary::COEF_DELTA( 1e-6 );

SingleResidueDunbrackLibrary::SingleResidueDunbrackLibrary(
	AA const aa,
	Size const n_rotameric_chi,
	bool dun02
) :
	dun02_( dun02 ),
	aa_( aa ),
	n_rotameric_chi_( n_rotameric_chi ),
	n_chi_bins_( n_rotameric_chi, 0 ),
	n_chi_products_( n_rotameric_chi, 1),
	n_packed_rots_( 0 ),
	n_possible_rots_( 0 ),
	prob_to_accumulate_buried_(0.98),
	prob_to_accumulate_nonburied_(0.95),
	packed_rotno_conversion_data_current_( false )
{
	// this builds on the hard coded hack bellow
	// since NCAAs are aa_unk we cannot hard code the information
	// alternativly it is added to the residue type paramater files
	read_options();
	if (aa_ != chemical::aa_unk) {
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

SingleResidueDunbrackLibrary::~SingleResidueDunbrackLibrary() {}

void
SingleResidueDunbrackLibrary::read_options()
{
	using namespace basic::options;
	if ( option[ OptionKeys::packing::dunbrack_prob_buried ].user() )
		prob_to_accumulate_buried( option[ OptionKeys::packing::dunbrack_prob_buried ]() );
	if ( option[ OptionKeys::packing::dunbrack_prob_nonburied ].user() )
		prob_to_accumulate_nonburied( option[ OptionKeys::packing::dunbrack_prob_nonburied ]() );
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
		case aa_ala: rot.resize( 0 ); break;
		case aa_cys: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_asp: rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 6; break;
		case aa_glu: rot.resize( 3 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = 6; break;
		case aa_phe: rot.resize( 2 ); rot[ 1 ] = 3; rot[ 2 ] = 6; break;
		case aa_gly: rot.resize( 0 ); break;
		case aa_his: rot.resize( 2 ); rot[ 1 ] = 3; rot[ 2 ] = 12; break;
		case aa_ile: rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
		case aa_lys: rot.resize( 4 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = rot[ 4 ] = 3; break;
		case aa_leu: rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
		case aa_met: rot.resize( 3 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = 3; break;
		case aa_asn: rot.resize( 2 ); rot[ 1 ] = 3; rot[ 2 ] = 12; break;
		case aa_pro: rot.resize( 3 ); rot[ 1 ] = 2; rot[ 2 ] = rot[ 3 ] = 1; break;
		case aa_gln: rot.resize( 3 ); rot[ 1 ] = rot[ 2 ] = 3; rot[ 3 ] = 12; break;
		case aa_arg: rot.resize( 4 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = rot[ 4 ] = 3; break;
		case aa_ser: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_thr: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_val: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_trp: rot.resize( 2 ); rot[ 1 ] = 3; rot[ 2 ] = 12; break;
		case aa_tyr: rot.resize( 2 ); rot[ 1 ] = 3; rot[ 2 ] = 6; break;
		default:
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
		case aa_ala: rot.resize( 0 ); break;
		case aa_cys: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_asp: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_glu: rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
		case aa_phe: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_gly: rot.resize( 0 ); break;
		case aa_his: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_ile: rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
		case aa_lys: rot.resize( 4 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = rot[ 4 ] = 3; break;
		case aa_leu: rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
		case aa_met: rot.resize( 3 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = 3; break;
		case aa_asn: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_pro: rot.resize( 3 ); rot[ 1 ] = 2; rot[ 2 ] = rot[ 3 ] = 1; break;
		case aa_gln: rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
		case aa_arg: rot.resize( 4 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = rot[ 4 ] = 3; break;
		case aa_ser: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_thr: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_val: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_trp: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_tyr: rot.resize( 1 ); rot[ 1 ] = 3; break;
		default:
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
		case aa_ala: rot.resize( 0 ); break;
		case aa_cys: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_asp: rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
		case aa_glu: rot.resize( 3 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = 3; break;
		case aa_phe: rot.resize( 2 ); rot[ 1 ] = 3; rot[ 2 ] = 2; break;
		case aa_gly: rot.resize( 0 ); break;
		case aa_his: rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
		case aa_ile: rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
		case aa_lys: rot.resize( 4 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = rot[ 4 ] = 3; break;
		case aa_leu: rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
		case aa_met: rot.resize( 3 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = 3; break;
		case aa_asn: rot.resize( 2 ); rot[ 1 ] = 3; rot[ 2 ] = 6; break;
		case aa_pro: rot.resize( 3 ); rot[ 1 ] = 2; rot[ 2 ] = rot[ 3 ] = 1; break;
		case aa_gln: rot.resize( 3 ); rot[ 1 ] = rot[ 2 ] = 3; rot[ 3 ] = 4; break;
		case aa_arg: rot.resize( 4 ); rot[ 1 ] = rot[ 2 ] = rot[ 3 ] = rot[ 4 ] = 3; break;
		case aa_ser: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_thr: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_val: rot.resize( 1 ); rot[ 1 ] = 3; break;
		case aa_trp: rot.resize( 2 ); rot[ 1 ] = rot[ 2 ] = 3; break;
		case aa_tyr: rot.resize( 2 ); rot[ 1 ] = 3; rot[ 2 ] = 2; break;
		default:
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
	if (rotno < 1 || rotno > rotno_2_packed_rotno_.size() ) return 0;
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
	if( ptr == 0 ) {
		TR << "In comparison operator: right-hand side is not a SingleResidueDunbrackLibrary." << std::endl;
		return false;
	}

	SingleResidueDunbrackLibrary const & other( dynamic_cast< SingleResidueDunbrackLibrary const & > ( rhs ) );

	// We don't call the parent == operator, as there's no meaningful data there to compare.
	bool equal( true );

	if( aa_ != other.aa_ ) {
		TR.Debug << "Comparison failure: " << core::chemical::name_from_aa(aa_)
				<< " does not match " << core::chemical::name_from_aa(other.aa_) << std::endl;
		return false; // Major data mismatch - don't bother reporting others.
	}
	if( dun02_ != other.dun02_ ) {
		TR.Debug << "Comparison failure: Dun02 setting mismatch: " << dun02_ << " vs. " << other.dun02_ << std::endl;
		return false; // Major data mismatch - don't bother reporting others.
	}

	/// Non-binary loaded data - check for equality for sanity purposes
	if( n_rotameric_chi_ != other.n_rotameric_chi_ ) {
		TR.Debug << "Comparison failure: n_rotameric_chi: " << n_rotameric_chi_ << " vs. " << other.n_rotameric_chi_ << std::endl;
		equal = false;
	}
	if( n_possible_rots_ != other.n_possible_rots_ ) {
		TR.Debug << "Comparison failure: n_possible_rots: " << n_possible_rots_ << " vs. " << other.n_possible_rots_ << std::endl;
		equal = false;
	}
	if( packed_rotno_conversion_data_current_ != other.packed_rotno_conversion_data_current_ ) {
		TR.Debug << "Comparison failure: packed_rotno_conversion_data_current: " << packed_rotno_conversion_data_current_ << " vs. " << other.packed_rotno_conversion_data_current_ << std::endl;
		equal = false;
	}
	if( ! numeric::equal_by_epsilon( prob_to_accumulate_buried_, other.prob_to_accumulate_buried_, PROB_DELTA) ) {
		TR.Debug << "Comparison failure: prob_to_accumulate_buried: " << prob_to_accumulate_buried_ << " vs. " << other.prob_to_accumulate_buried_ << std::endl;
		equal = false;
	}
	if( ! numeric::equal_by_epsilon( prob_to_accumulate_nonburied_, other.prob_to_accumulate_nonburied_, PROB_DELTA) ) {
		TR.Debug << "Comparison failure: prob_to_accumulate_nonburied: " << prob_to_accumulate_nonburied_ << " vs. " << other.prob_to_accumulate_nonburied_ << std::endl;
		equal = false;
	}

	if( n_chi_bins_.size() != other.n_chi_bins_.size() ) {
			TR.Debug << "Comparison failure: n_chi_bins vector length: " << n_chi_bins_.size() << " vs. " << other.n_chi_bins_.size() << std::endl;
			equal = false;
	} else {
		for( core::Size ii(1); ii <= n_chi_bins_.size(); ++ii ) {
			if( n_chi_bins_[ii] != other.n_chi_bins_[ii] ) {
				TR.Debug << "Comparison failure: n_chi_bins: " << n_chi_bins_[ii] << " vs. " << other.n_chi_bins_[ii] << std::endl;
				equal = false;
			}
		}
	}
	if( n_chi_products_.size() != other.n_chi_products_.size() ) {
			TR.Debug << "Comparison failure: n_chi_products vector length: " << n_chi_products_.size() << " vs. " << other.n_chi_products_.size() << std::endl;
			equal = false;
	} else {
		for( core::Size ii(1); ii <= n_chi_products_.size(); ++ii ) {
			if( n_chi_products_[ii] != other.n_chi_products_[ii] ) {
				TR.Debug << "Comparison failure: n_chi_products: " << n_chi_products_[ii] << " vs. " << other.n_chi_products_[ii] << std::endl;
				equal = false;
			}
		}
	}
	// // For some reason, rotwell_exists_ is not matching on direct re-load comparison (empty vector on ASCII reload ?)
	//if( rotwell_exists_.size() != other.rotwell_exists_.size() ) {
	//		TR.Debug << "Comparison failure: rotwell_exists vector length: " << rotwell_exists_.size() << " vs. " << other.rotwell_exists_.size() << std::endl;
	//		equal = false;
	//} else {
	//	for( core::Size ii(1); ii <= rotwell_exists_.size(); ++ii ) {
	//		if( rotwell_exists_[ii] != other.rotwell_exists_[ii] ) {
	//			TR.Debug << "Comparison failure: rotwell_exists: " << rotwell_exists_[ii] << " vs. " << other.rotwell_exists_[ii] << std::endl;
	//			equal = false;
	//		}
	//	}
	//}

	/// Binary-loaded data.

	/// 1. n_packed_rots_
	if( n_packed_rots_ != other.n_packed_rots_ ) {
		TR.Debug << "Comparsion failure in " << core::chemical::name_from_aa(aa_ )
				<< ": n_packed_rots is unequal. " << n_packed_rots_ << " vs. " << other.n_packed_rots_ << std::endl;
		equal = false;
	}

	/// 2. rotno_2_packed_rotno_
	assert( n_possible_rots_ == other.n_possible_rots_ ); // This doesn't vary?
	for ( Size ii = 1; ii <= n_possible_rots_; ++ii ) {
		if( rotno_2_packed_rotno_[ ii ] != other.rotno_2_packed_rotno_[ ii ] ) {
			TR.Debug << "Comparsion failure in " << core::chemical::name_from_aa(aa_ )
					<< ": rotno_2_packed_rotno " << ii << " - "
					<< rotno_2_packed_rotno_[ ii ] << " vs. " << other.rotno_2_packed_rotno_[ ii ] << std::endl;
			equal = false;
		}
	}

	/// 3. packed_rotno_2_rotno_
	if( n_packed_rots_ == other.n_packed_rots_ ) {
		for ( Size ii = 1; ii <= n_packed_rots_; ++ii ) {
			if( packed_rotno_2_rotno_[ ii ] != other.packed_rotno_2_rotno_[ ii ] ) {
				TR.Debug << "Comparsion failure in " << core::chemical::name_from_aa(aa_ )
						<< ": packed_rotno_2_rotno " << ii << " - "
						<< packed_rotno_2_rotno_[ ii ] << " vs. " << other.packed_rotno_2_rotno_[ ii ] << std::endl;
				equal = false;
			}
		}
	}

	/// 4. packed_rotno_2_rotwell_
	assert( n_rotameric_chi_ == other.n_rotameric_chi_ );
	if( n_packed_rots_ == other.n_packed_rots_ ) {
		for ( Size ii = 1; ii <= n_packed_rots_; ++ii ) {
			for ( Size jj = 1; jj <= n_rotameric_chi_; ++jj ) {
				if( packed_rotno_2_rotwell_[ ii ][ jj ] != other.packed_rotno_2_rotwell_[ ii ][ jj ] ) {
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
debug_assert( false );
	utility_exit_with_message(
	"ERROR: SingleResidueDunbrackLibrary::hokey_template_workaround should never be called!");

	RotamericSingleResidueDunbrackLibrary< ONE, ONE >   rsrdl_11( chemical::aa_ala, false );
	RotamericSingleResidueDunbrackLibrary< ONE, TWO >   rsrdl_12( chemical::aa_ala, false );
	RotamericSingleResidueDunbrackLibrary< ONE, THREE > rsrdl_13( chemical::aa_ala, false );
	RotamericSingleResidueDunbrackLibrary< ONE, FOUR >  rsrdl_14( chemical::aa_ala, false );
	RotamericSingleResidueDunbrackLibrary< ONE, FIVE >  rsrdl_15( chemical::aa_ala, false );
	RotamericSingleResidueDunbrackLibrary< TWO, ONE >   rsrdl_21( chemical::aa_ala, false );
	RotamericSingleResidueDunbrackLibrary< TWO, TWO >   rsrdl_22( chemical::aa_ala, false );
	RotamericSingleResidueDunbrackLibrary< TWO, THREE > rsrdl_23( chemical::aa_ala, false );
	RotamericSingleResidueDunbrackLibrary< TWO, FOUR >  rsrdl_24( chemical::aa_ala, false );
	RotamericSingleResidueDunbrackLibrary< TWO, FIVE >  rsrdl_25( chemical::aa_ala, false );
	RotamericSingleResidueDunbrackLibrary< THREE, ONE >   rsrdl_31( chemical::aa_ala, false );
	RotamericSingleResidueDunbrackLibrary< THREE, TWO >   rsrdl_32( chemical::aa_ala, false );
	RotamericSingleResidueDunbrackLibrary< THREE, THREE > rsrdl_33( chemical::aa_ala, false );
	RotamericSingleResidueDunbrackLibrary< THREE, FOUR >  rsrdl_34( chemical::aa_ala, false );
	RotamericSingleResidueDunbrackLibrary< THREE, FIVE >  rsrdl_35( chemical::aa_ala, false );
	RotamericSingleResidueDunbrackLibrary< FOUR, ONE >   rsrdl_41( chemical::aa_ala, false );
	RotamericSingleResidueDunbrackLibrary< FOUR, TWO >   rsrdl_42( chemical::aa_ala, false );
	RotamericSingleResidueDunbrackLibrary< FOUR, THREE > rsrdl_43( chemical::aa_ala, false );
	RotamericSingleResidueDunbrackLibrary< FOUR, FOUR >  rsrdl_44( chemical::aa_ala, false );
	RotamericSingleResidueDunbrackLibrary< FOUR, FIVE >  rsrdl_45( chemical::aa_ala, false );

	SemiRotamericSingleResidueDunbrackLibrary< ONE, ONE >   srsrdl_11( chemical::aa_ala, true, true ); //e.g. asn,phe
	SemiRotamericSingleResidueDunbrackLibrary< TWO, ONE >   srsrdl_12( chemical::aa_ala, true, true ); // e.g. glu
	SemiRotamericSingleResidueDunbrackLibrary< ONE, TWO >   srsrdl_21( chemical::aa_ala, true, true ); //e.g. asn,phe
	SemiRotamericSingleResidueDunbrackLibrary< TWO, TWO >   srsrdl_22( chemical::aa_ala, true, true ); // e.g. glu
	SemiRotamericSingleResidueDunbrackLibrary< ONE, THREE >   srsrdl_31( chemical::aa_ala, true, true ); //e.g. asn,phe
	SemiRotamericSingleResidueDunbrackLibrary< TWO, THREE >   srsrdl_32( chemical::aa_ala, true, true ); // e.g. glu
	SemiRotamericSingleResidueDunbrackLibrary< ONE, FOUR >   srsrdl_41( chemical::aa_ala, true, true ); //e.g. asn,phe
	SemiRotamericSingleResidueDunbrackLibrary< TWO, FOUR >   srsrdl_42( chemical::aa_ala, true, true ); // e.g. glu
	SemiRotamericSingleResidueDunbrackLibrary< ONE, FIVE >   srsrdl_51( chemical::aa_ala, true, true ); //e.g. asn,phe
	SemiRotamericSingleResidueDunbrackLibrary< TWO, FIVE >   srsrdl_52( chemical::aa_ala, true, true ); // e.g. glu
	// three and four do not exist... they could in the future if needed.

	chemical::ResidueType rt( NULL, NULL, NULL, NULL );
	conformation::Residue rsd( rt, true );
	RotamerLibraryScratchSpace scratch;
	Size4 rotwell;
	Size i(0);
	utility::fixedsizearray1< Real, ONE > bb1;
	utility::fixedsizearray1< Real, TWO > bb2;
	utility::fixedsizearray1< Real, THREE > bb3;
	utility::fixedsizearray1< Real, FOUR > bb4;
	utility::fixedsizearray1< Real, FIVE > bb5;

	PackedDunbrackRotamer< ONE, ONE, Real > prot11;
	PackedDunbrackRotamer< ONE, TWO, Real > prot12;
	PackedDunbrackRotamer< ONE, THREE, Real > prot13;
	PackedDunbrackRotamer< ONE, FOUR, Real > prot14;
	PackedDunbrackRotamer< ONE, FIVE, Real > prot15;
	PackedDunbrackRotamer< TWO, ONE, Real > prot21;
	PackedDunbrackRotamer< TWO, TWO, Real > prot22;
	PackedDunbrackRotamer< TWO, THREE, Real > prot23;
	PackedDunbrackRotamer< TWO, FOUR, Real > prot24;
	PackedDunbrackRotamer< TWO, FIVE, Real > prot25;
	PackedDunbrackRotamer< THREE, ONE, Real > prot31;
	PackedDunbrackRotamer< THREE, TWO, Real > prot32;
	PackedDunbrackRotamer< THREE, THREE, Real > prot33;
	PackedDunbrackRotamer< THREE, FOUR, Real > prot34;
	PackedDunbrackRotamer< THREE, FIVE, Real > prot35;
	PackedDunbrackRotamer< FOUR, ONE, Real > prot41;
	PackedDunbrackRotamer< FOUR, TWO, Real > prot42;
	PackedDunbrackRotamer< FOUR, THREE, Real > prot43;
	PackedDunbrackRotamer< FOUR, FOUR, Real > prot44;
	PackedDunbrackRotamer< FOUR, FIVE, Real > prot45;

	rsrdl_11.nchi();
	rsrdl_12.nchi();
	rsrdl_13.nchi();
	rsrdl_14.nchi();
	rsrdl_15.nchi();
	srsrdl_11.nchi();
	srsrdl_12.nchi();
	rsrdl_21.nchi();
	rsrdl_22.nchi();
	rsrdl_23.nchi();
	rsrdl_24.nchi();
	rsrdl_25.nchi();
	srsrdl_21.nchi();
	srsrdl_22.nchi();
	rsrdl_31.nchi();
	rsrdl_32.nchi();
	rsrdl_33.nchi();
	rsrdl_34.nchi();
	rsrdl_35.nchi();
	srsrdl_31.nchi();
	srsrdl_32.nchi();
	rsrdl_41.nchi();
	rsrdl_42.nchi();
	rsrdl_43.nchi();
	rsrdl_44.nchi();
	rsrdl_45.nchi();
	srsrdl_41.nchi();
	srsrdl_42.nchi();
	srsrdl_51.nchi();
	srsrdl_52.nchi();

	rsrdl_11.n_rotamer_bins();
	rsrdl_12.n_rotamer_bins();
	rsrdl_13.n_rotamer_bins();
	rsrdl_14.n_rotamer_bins();
	rsrdl_15.n_rotamer_bins();
	srsrdl_11.n_rotamer_bins();
	srsrdl_12.n_rotamer_bins();
	rsrdl_21.n_rotamer_bins();
	rsrdl_22.n_rotamer_bins();
	rsrdl_23.n_rotamer_bins();
	rsrdl_24.n_rotamer_bins();
	rsrdl_25.n_rotamer_bins();
	srsrdl_21.n_rotamer_bins();
	srsrdl_22.n_rotamer_bins();
	rsrdl_31.n_rotamer_bins();
	rsrdl_32.n_rotamer_bins();
	rsrdl_33.n_rotamer_bins();
	rsrdl_34.n_rotamer_bins();
	rsrdl_35.n_rotamer_bins();
	srsrdl_31.n_rotamer_bins();
	srsrdl_32.n_rotamer_bins();
	rsrdl_41.n_rotamer_bins();
	rsrdl_42.n_rotamer_bins();
	rsrdl_43.n_rotamer_bins();
	rsrdl_44.n_rotamer_bins();
	rsrdl_45.n_rotamer_bins();
	srsrdl_41.n_rotamer_bins();
	srsrdl_42.n_rotamer_bins();
	srsrdl_51.n_rotamer_bins();
	srsrdl_52.n_rotamer_bins();

	rsrdl_11.rotamer_energy( rsd, scratch );
	rsrdl_12.rotamer_energy( rsd, scratch );
	rsrdl_13.rotamer_energy( rsd, scratch );
	rsrdl_14.rotamer_energy( rsd, scratch );
	rsrdl_15.rotamer_energy( rsd, scratch );
	srsrdl_11.rotamer_energy( rsd, scratch );
	srsrdl_12.rotamer_energy( rsd, scratch );
	rsrdl_21.rotamer_energy( rsd, scratch );
	rsrdl_22.rotamer_energy( rsd, scratch );
	rsrdl_23.rotamer_energy( rsd, scratch );
	rsrdl_24.rotamer_energy( rsd, scratch );
	rsrdl_25.rotamer_energy( rsd, scratch );
	srsrdl_21.rotamer_energy( rsd, scratch );
	srsrdl_22.rotamer_energy( rsd, scratch );
	rsrdl_31.rotamer_energy( rsd, scratch );
	rsrdl_32.rotamer_energy( rsd, scratch );
	rsrdl_33.rotamer_energy( rsd, scratch );
	rsrdl_34.rotamer_energy( rsd, scratch );
	rsrdl_35.rotamer_energy( rsd, scratch );
	srsrdl_31.rotamer_energy( rsd, scratch );
	srsrdl_32.rotamer_energy( rsd, scratch );
	rsrdl_41.rotamer_energy( rsd, scratch );
	rsrdl_42.rotamer_energy( rsd, scratch );
	rsrdl_43.rotamer_energy( rsd, scratch );
	rsrdl_44.rotamer_energy( rsd, scratch );
	rsrdl_45.rotamer_energy( rsd, scratch );
	srsrdl_41.rotamer_energy( rsd, scratch );
	srsrdl_42.rotamer_energy( rsd, scratch );
	srsrdl_51.rotamer_energy( rsd, scratch );
	srsrdl_52.rotamer_energy( rsd, scratch );

	rsrdl_11.rotamer_energy_deriv( rsd, scratch );
	rsrdl_12.rotamer_energy_deriv( rsd, scratch );
	rsrdl_13.rotamer_energy_deriv( rsd, scratch );
	rsrdl_14.rotamer_energy_deriv( rsd, scratch );
	rsrdl_15.rotamer_energy_deriv( rsd, scratch );
	srsrdl_11.rotamer_energy_deriv( rsd, scratch );
	srsrdl_12.rotamer_energy_deriv( rsd, scratch );
	rsrdl_21.rotamer_energy_deriv( rsd, scratch );
	rsrdl_22.rotamer_energy_deriv( rsd, scratch );
	rsrdl_23.rotamer_energy_deriv( rsd, scratch );
	rsrdl_24.rotamer_energy_deriv( rsd, scratch );
	rsrdl_25.rotamer_energy_deriv( rsd, scratch );
	srsrdl_21.rotamer_energy_deriv( rsd, scratch );
	srsrdl_22.rotamer_energy_deriv( rsd, scratch );
	rsrdl_31.rotamer_energy_deriv( rsd, scratch );
	rsrdl_32.rotamer_energy_deriv( rsd, scratch );
	rsrdl_33.rotamer_energy_deriv( rsd, scratch );
	rsrdl_34.rotamer_energy_deriv( rsd, scratch );
	rsrdl_35.rotamer_energy_deriv( rsd, scratch );
	srsrdl_31.rotamer_energy_deriv( rsd, scratch );
	srsrdl_32.rotamer_energy_deriv( rsd, scratch );
	rsrdl_41.rotamer_energy_deriv( rsd, scratch );
	rsrdl_42.rotamer_energy_deriv( rsd, scratch );
	rsrdl_43.rotamer_energy_deriv( rsd, scratch );
	rsrdl_44.rotamer_energy_deriv( rsd, scratch );
	rsrdl_45.rotamer_energy_deriv( rsd, scratch );
	srsrdl_41.rotamer_energy_deriv( rsd, scratch );
	srsrdl_42.rotamer_energy_deriv( rsd, scratch );
	srsrdl_51.rotamer_energy_deriv( rsd, scratch );
	srsrdl_52.rotamer_energy_deriv( rsd, scratch );

	rsrdl_11.best_rotamer_energy( rsd, true, scratch );
	rsrdl_12.best_rotamer_energy( rsd, true, scratch );
	rsrdl_13.best_rotamer_energy( rsd, true, scratch );
	rsrdl_14.best_rotamer_energy( rsd, true, scratch );
	rsrdl_15.best_rotamer_energy( rsd, true, scratch );
	srsrdl_11.best_rotamer_energy( rsd, true, scratch );
	srsrdl_12.best_rotamer_energy( rsd, true, scratch );
	rsrdl_21.best_rotamer_energy( rsd, true, scratch );
	rsrdl_22.best_rotamer_energy( rsd, true, scratch );
	rsrdl_23.best_rotamer_energy( rsd, true, scratch );
	rsrdl_24.best_rotamer_energy( rsd, true, scratch );
	rsrdl_25.best_rotamer_energy( rsd, true, scratch );
	srsrdl_21.best_rotamer_energy( rsd, true, scratch );
	srsrdl_22.best_rotamer_energy( rsd, true, scratch );
	rsrdl_31.best_rotamer_energy( rsd, true, scratch );
	rsrdl_32.best_rotamer_energy( rsd, true, scratch );
	rsrdl_33.best_rotamer_energy( rsd, true, scratch );
	rsrdl_34.best_rotamer_energy( rsd, true, scratch );
	rsrdl_35.best_rotamer_energy( rsd, true, scratch );
	srsrdl_31.best_rotamer_energy( rsd, true, scratch );
	srsrdl_32.best_rotamer_energy( rsd, true, scratch );
	rsrdl_41.best_rotamer_energy( rsd, true, scratch );
	rsrdl_42.best_rotamer_energy( rsd, true, scratch );
	rsrdl_43.best_rotamer_energy( rsd, true, scratch );
	rsrdl_44.best_rotamer_energy( rsd, true, scratch );
	rsrdl_45.best_rotamer_energy( rsd, true, scratch );
	srsrdl_41.best_rotamer_energy( rsd, true, scratch );
	srsrdl_42.best_rotamer_energy( rsd, true, scratch );
	srsrdl_51.best_rotamer_energy( rsd, true, scratch );
	srsrdl_52.best_rotamer_energy( rsd, true, scratch );

	pose::Pose pose;
	scoring::ScoreFunction sfxn;
	pack::task::PackerTask_ task( pose );

	core::graph::GraphOP png;
	chemical::ResidueTypeCOP cr;
	utility::vector1< utility::vector1< Real > > ecs;
	rotamers::RotamerVector rv;

	rsrdl_11.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	rsrdl_12.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	rsrdl_13.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	rsrdl_14.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	rsrdl_15.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	srsrdl_11.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	srsrdl_12.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	rsrdl_21.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	rsrdl_22.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	rsrdl_23.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	rsrdl_24.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	rsrdl_25.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	srsrdl_21.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	srsrdl_22.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	rsrdl_31.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	rsrdl_32.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	rsrdl_33.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	rsrdl_34.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	rsrdl_35.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	srsrdl_31.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	srsrdl_32.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	rsrdl_41.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	rsrdl_42.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	rsrdl_43.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	rsrdl_44.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	rsrdl_45.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	srsrdl_41.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	srsrdl_42.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	srsrdl_51.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );
	srsrdl_52.fill_rotamer_vector( pose, sfxn, task, png, cr, rsd, ecs, true, rv );

	utility::io::ozstream ozs;
	rsrdl_11.write_to_file( ozs );
	rsrdl_12.write_to_file( ozs );
	rsrdl_13.write_to_file( ozs );
	rsrdl_14.write_to_file( ozs );
	rsrdl_15.write_to_file( ozs );
	srsrdl_11.write_to_file( ozs );
	srsrdl_12.write_to_file( ozs );
	rsrdl_21.write_to_file( ozs );
	rsrdl_22.write_to_file( ozs );
	rsrdl_23.write_to_file( ozs );
	rsrdl_24.write_to_file( ozs );
	rsrdl_25.write_to_file( ozs );
	srsrdl_21.write_to_file( ozs );
	srsrdl_22.write_to_file( ozs );
	rsrdl_31.write_to_file( ozs );
	rsrdl_32.write_to_file( ozs );
	rsrdl_33.write_to_file( ozs );
	rsrdl_34.write_to_file( ozs );
	rsrdl_35.write_to_file( ozs );
	srsrdl_31.write_to_file( ozs );
	srsrdl_32.write_to_file( ozs );
	rsrdl_41.write_to_file( ozs );
	rsrdl_42.write_to_file( ozs );
	rsrdl_43.write_to_file( ozs );
	rsrdl_44.write_to_file( ozs );
	rsrdl_45.write_to_file( ozs );
	srsrdl_41.write_to_file( ozs );
	srsrdl_42.write_to_file( ozs );
	srsrdl_51.write_to_file( ozs );
	srsrdl_52.write_to_file( ozs );

	utility::io::ozstream os;
	rsrdl_11.write_to_binary( os );
	rsrdl_12.write_to_binary( os );
	rsrdl_13.write_to_binary( os );
	rsrdl_14.write_to_binary( os );
	rsrdl_15.write_to_binary( os );
	srsrdl_11.write_to_binary( os );
	srsrdl_12.write_to_binary( os );
	rsrdl_21.write_to_binary( os );
	rsrdl_22.write_to_binary( os );
	rsrdl_23.write_to_binary( os );
	rsrdl_24.write_to_binary( os );
	rsrdl_25.write_to_binary( os );
	srsrdl_21.write_to_binary( os );
	srsrdl_22.write_to_binary( os );
	rsrdl_31.write_to_binary( os );
	rsrdl_32.write_to_binary( os );
	rsrdl_33.write_to_binary( os );
	rsrdl_34.write_to_binary( os );
	rsrdl_35.write_to_binary( os );
	srsrdl_31.write_to_binary( os );
	srsrdl_32.write_to_binary( os );
	rsrdl_41.write_to_binary( os );
	rsrdl_42.write_to_binary( os );
	rsrdl_43.write_to_binary( os );
	rsrdl_44.write_to_binary( os );
	rsrdl_45.write_to_binary( os );
	srsrdl_41.write_to_binary( os );
	srsrdl_42.write_to_binary( os );
	srsrdl_51.write_to_binary( os );
	srsrdl_52.write_to_binary( os );

	utility::io::izstream is;
	rsrdl_11.read_from_binary( is );
	rsrdl_12.read_from_binary( is );
	rsrdl_13.read_from_binary( is );
	rsrdl_14.read_from_binary( is );
	rsrdl_15.read_from_binary( is );
	srsrdl_11.read_from_binary( is );
	srsrdl_12.read_from_binary( is );
	rsrdl_21.read_from_binary( is );
	rsrdl_22.read_from_binary( is );
	rsrdl_23.read_from_binary( is );
	rsrdl_24.read_from_binary( is );
	rsrdl_25.read_from_binary( is );
	srsrdl_21.read_from_binary( is );
	srsrdl_22.read_from_binary( is );
	rsrdl_31.read_from_binary( is );
	rsrdl_32.read_from_binary( is );
	rsrdl_33.read_from_binary( is );
	rsrdl_34.read_from_binary( is );
	rsrdl_35.read_from_binary( is );
	srsrdl_31.read_from_binary( is );
	srsrdl_32.read_from_binary( is );
	rsrdl_41.read_from_binary( is );
	rsrdl_42.read_from_binary( is );
	rsrdl_43.read_from_binary( is );
	rsrdl_44.read_from_binary( is );
	rsrdl_45.read_from_binary( is );
	srsrdl_41.read_from_binary( is );
	srsrdl_42.read_from_binary( is );
	srsrdl_51.read_from_binary( is );
	srsrdl_52.read_from_binary( is );

	rsrdl_11.memory_usage_in_bytes();
	rsrdl_12.memory_usage_in_bytes();
	rsrdl_13.memory_usage_in_bytes();
	rsrdl_14.memory_usage_in_bytes();
	rsrdl_15.memory_usage_in_bytes();
	srsrdl_11.memory_usage_in_bytes();
	srsrdl_12.memory_usage_in_bytes();
	rsrdl_21.memory_usage_in_bytes();
	rsrdl_22.memory_usage_in_bytes();
	rsrdl_23.memory_usage_in_bytes();
	rsrdl_24.memory_usage_in_bytes();
	rsrdl_25.memory_usage_in_bytes();
	srsrdl_21.memory_usage_in_bytes();
	srsrdl_22.memory_usage_in_bytes();
	rsrdl_31.memory_usage_in_bytes();
	rsrdl_32.memory_usage_in_bytes();
	rsrdl_33.memory_usage_in_bytes();
	rsrdl_34.memory_usage_in_bytes();
	rsrdl_35.memory_usage_in_bytes();
	srsrdl_31.memory_usage_in_bytes();
	srsrdl_32.memory_usage_in_bytes();
	rsrdl_41.memory_usage_in_bytes();
	rsrdl_42.memory_usage_in_bytes();
	rsrdl_43.memory_usage_in_bytes();
	rsrdl_44.memory_usage_in_bytes();
	rsrdl_45.memory_usage_in_bytes();
	srsrdl_41.memory_usage_in_bytes();
	srsrdl_42.memory_usage_in_bytes();
	srsrdl_51.memory_usage_in_bytes();
	srsrdl_52.memory_usage_in_bytes();

	utility::vector1< Real > chi; utility::vector1< Size > rot;
	rsrdl_11.get_rotamer_from_chi( chi, rot );
	rsrdl_12.get_rotamer_from_chi( chi, rot );
	rsrdl_13.get_rotamer_from_chi( chi, rot );
	rsrdl_14.get_rotamer_from_chi( chi, rot );
	rsrdl_15.get_rotamer_from_chi( chi, rot );
	srsrdl_11.get_rotamer_from_chi( chi, rot );
	srsrdl_12.get_rotamer_from_chi( chi, rot );
	rsrdl_21.get_rotamer_from_chi( chi, rot );
	rsrdl_22.get_rotamer_from_chi( chi, rot );
	rsrdl_23.get_rotamer_from_chi( chi, rot );
	rsrdl_24.get_rotamer_from_chi( chi, rot );
	rsrdl_25.get_rotamer_from_chi( chi, rot );
	srsrdl_21.get_rotamer_from_chi( chi, rot );
	srsrdl_22.get_rotamer_from_chi( chi, rot );
	rsrdl_31.get_rotamer_from_chi( chi, rot );
	rsrdl_32.get_rotamer_from_chi( chi, rot );
	rsrdl_33.get_rotamer_from_chi( chi, rot );
	rsrdl_34.get_rotamer_from_chi( chi, rot );
	rsrdl_35.get_rotamer_from_chi( chi, rot );
	srsrdl_31.get_rotamer_from_chi( chi, rot );
	srsrdl_32.get_rotamer_from_chi( chi, rot );
	rsrdl_41.get_rotamer_from_chi( chi, rot );
	rsrdl_42.get_rotamer_from_chi( chi, rot );
	rsrdl_43.get_rotamer_from_chi( chi, rot );
	rsrdl_44.get_rotamer_from_chi( chi, rot );
	rsrdl_45.get_rotamer_from_chi( chi, rot );
	srsrdl_41.get_rotamer_from_chi( chi, rot );
	srsrdl_42.get_rotamer_from_chi( chi, rot );
	srsrdl_51.get_rotamer_from_chi( chi, rot );
	srsrdl_52.get_rotamer_from_chi( chi, rot );

	rsrdl_11.find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	rsrdl_12.find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	rsrdl_13.find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	rsrdl_14.find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	rsrdl_15.find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	srsrdl_11.find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	srsrdl_12.find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	rsrdl_21.find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	rsrdl_22.find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	rsrdl_23.find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	rsrdl_24.find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	rsrdl_25.find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	srsrdl_21.find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	srsrdl_22.find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	rsrdl_31.find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	rsrdl_32.find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	rsrdl_33.find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	rsrdl_34.find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	rsrdl_35.find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	srsrdl_31.find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	srsrdl_32.find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	rsrdl_41.find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	rsrdl_42.find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	rsrdl_43.find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	rsrdl_44.find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	rsrdl_45.find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	srsrdl_41.find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	srsrdl_42.find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	srsrdl_51.find_another_representative_for_unlikely_rotamer( rsd, rotwell );
	srsrdl_52.find_another_representative_for_unlikely_rotamer( rsd, rotwell );

	rsrdl_11.interpolate_rotamers( rsd, scratch, i, prot11 );
	rsrdl_12.interpolate_rotamers( rsd, scratch, i, prot12 );
	rsrdl_13.interpolate_rotamers( rsd, scratch, i, prot13 );
	rsrdl_14.interpolate_rotamers( rsd, scratch, i, prot14 );
	rsrdl_15.interpolate_rotamers( rsd, scratch, i, prot15 );
	srsrdl_11.interpolate_rotamers( rsd, scratch, i, prot11 );
	srsrdl_21.interpolate_rotamers( rsd, scratch, i, prot12 );
	srsrdl_31.interpolate_rotamers( rsd, scratch, i, prot13 );
	srsrdl_41.interpolate_rotamers( rsd, scratch, i, prot14 );
	srsrdl_51.interpolate_rotamers( rsd, scratch, i, prot15 );
	rsrdl_21.interpolate_rotamers( rsd, scratch, i, prot21 );
	rsrdl_22.interpolate_rotamers( rsd, scratch, i, prot22 );
	rsrdl_23.interpolate_rotamers( rsd, scratch, i, prot23 );
	rsrdl_24.interpolate_rotamers( rsd, scratch, i, prot24 );
	rsrdl_25.interpolate_rotamers( rsd, scratch, i, prot25 );
	srsrdl_12.interpolate_rotamers( rsd, scratch, i, prot21 );
	srsrdl_22.interpolate_rotamers( rsd, scratch, i, prot22 );
	srsrdl_32.interpolate_rotamers( rsd, scratch, i, prot23 );
	srsrdl_42.interpolate_rotamers( rsd, scratch, i, prot24 );
	srsrdl_52.interpolate_rotamers( rsd, scratch, i, prot25 );
	rsrdl_31.interpolate_rotamers( rsd, scratch, i, prot31 );
	rsrdl_32.interpolate_rotamers( rsd, scratch, i, prot32 );
	rsrdl_33.interpolate_rotamers( rsd, scratch, i, prot33 );
	rsrdl_34.interpolate_rotamers( rsd, scratch, i, prot34 );
	rsrdl_35.interpolate_rotamers( rsd, scratch, i, prot35 );
	rsrdl_41.interpolate_rotamers( rsd, scratch, i, prot41 );
	rsrdl_42.interpolate_rotamers( rsd, scratch, i, prot42 );
	rsrdl_43.interpolate_rotamers( rsd, scratch, i, prot43 );
	rsrdl_44.interpolate_rotamers( rsd, scratch, i, prot44 );
	rsrdl_45.interpolate_rotamers( rsd, scratch, i, prot45 );

	rsrdl_11.memory_usage_dynamic();
	rsrdl_12.memory_usage_dynamic();
	rsrdl_13.memory_usage_dynamic();
	rsrdl_14.memory_usage_dynamic();
	rsrdl_15.memory_usage_dynamic();
	srsrdl_11.memory_usage_dynamic();
	srsrdl_12.memory_usage_dynamic();
	rsrdl_21.memory_usage_dynamic();
	rsrdl_22.memory_usage_dynamic();
	rsrdl_23.memory_usage_dynamic();
	rsrdl_24.memory_usage_dynamic();
	rsrdl_25.memory_usage_dynamic();
	srsrdl_21.memory_usage_dynamic();
	srsrdl_22.memory_usage_dynamic();
	rsrdl_31.memory_usage_dynamic();
	rsrdl_32.memory_usage_dynamic();
	rsrdl_33.memory_usage_dynamic();
	rsrdl_34.memory_usage_dynamic();
	rsrdl_35.memory_usage_dynamic();
	srsrdl_31.memory_usage_dynamic();
	srsrdl_32.memory_usage_dynamic();
	rsrdl_41.memory_usage_dynamic();
	rsrdl_42.memory_usage_dynamic();
	rsrdl_43.memory_usage_dynamic();
	rsrdl_44.memory_usage_dynamic();
	rsrdl_45.memory_usage_dynamic();
	srsrdl_41.memory_usage_dynamic();
	srsrdl_42.memory_usage_dynamic();
	srsrdl_51.memory_usage_dynamic();
	srsrdl_52.memory_usage_dynamic();

	rsrdl_11.memory_usage_static();
	rsrdl_12.memory_usage_static();
	rsrdl_13.memory_usage_static();
	rsrdl_14.memory_usage_static();
	rsrdl_15.memory_usage_static();
	srsrdl_11.memory_usage_static();
	srsrdl_12.memory_usage_static();
	rsrdl_21.memory_usage_static();
	rsrdl_22.memory_usage_static();
	rsrdl_23.memory_usage_static();
	rsrdl_24.memory_usage_static();
	rsrdl_25.memory_usage_static();
	srsrdl_21.memory_usage_static();
	srsrdl_22.memory_usage_static();
	rsrdl_31.memory_usage_static();
	rsrdl_32.memory_usage_static();
	rsrdl_33.memory_usage_static();
	rsrdl_34.memory_usage_static();
	rsrdl_35.memory_usage_static();
	srsrdl_31.memory_usage_static();
	srsrdl_32.memory_usage_static();
	rsrdl_41.memory_usage_static();
	rsrdl_42.memory_usage_static();
	rsrdl_43.memory_usage_static();
	rsrdl_44.memory_usage_static();
	rsrdl_45.memory_usage_static();
	srsrdl_41.memory_usage_static();
	srsrdl_42.memory_usage_static();
	srsrdl_51.memory_usage_static();
	srsrdl_52.memory_usage_static();

	/*rsrdl_11.get_all_rotamer_samples( bb1 );
	rsrdl_12.get_all_rotamer_samples( bb2 );
	rsrdl_13.get_all_rotamer_samples( bb3 );
	rsrdl_14.get_all_rotamer_samples( bb4 );
	rsrdl_15.get_all_rotamer_samples( bb5 );
	srsrdl_11.get_all_rotamer_samples( bb1 );
	srsrdl_12.get_all_rotamer_samples( bb1 );
	rsrdl_21.get_all_rotamer_samples( bb1 );
	rsrdl_22.get_all_rotamer_samples( bb2 );
	rsrdl_23.get_all_rotamer_samples( bb3 );
	rsrdl_24.get_all_rotamer_samples( bb4 );
	rsrdl_25.get_all_rotamer_samples( bb5 );
	srsrdl_21.get_all_rotamer_samples( bb2 );
	srsrdl_22.get_all_rotamer_samples( bb2 );
	rsrdl_31.get_all_rotamer_samples( bb1 );
	rsrdl_32.get_all_rotamer_samples( bb2 );
	rsrdl_33.get_all_rotamer_samples( bb3 );
	rsrdl_34.get_all_rotamer_samples( bb4 );
	rsrdl_35.get_all_rotamer_samples( bb5 );
	srsrdl_31.get_all_rotamer_samples( bb3 );
	srsrdl_32.get_all_rotamer_samples( bb3 );
	rsrdl_41.get_all_rotamer_samples( bb1 );
	rsrdl_42.get_all_rotamer_samples( bb2 );
	rsrdl_43.get_all_rotamer_samples( bb3 );
	rsrdl_44.get_all_rotamer_samples( bb4 );
	rsrdl_45.get_all_rotamer_samples( bb5 );
	srsrdl_41.get_all_rotamer_samples( bb4 );
	srsrdl_42.get_all_rotamer_samples( bb4 );
	srsrdl_51.get_all_rotamer_samples( bb5 );
	srsrdl_52.get_all_rotamer_samples( bb5 );
*/
	rsrdl_11.get_all_rotamer_samples( bb5 );
	rsrdl_12.get_all_rotamer_samples( bb5 );
	rsrdl_13.get_all_rotamer_samples( bb5 );
	rsrdl_14.get_all_rotamer_samples( bb5 );
	rsrdl_15.get_all_rotamer_samples( bb5 );
	srsrdl_11.get_all_rotamer_samples( bb5 );
	srsrdl_12.get_all_rotamer_samples( bb5 );
	rsrdl_21.get_all_rotamer_samples( bb5 );
	rsrdl_22.get_all_rotamer_samples( bb5 );
	rsrdl_23.get_all_rotamer_samples( bb5 );
	rsrdl_24.get_all_rotamer_samples( bb5 );
	rsrdl_25.get_all_rotamer_samples( bb5 );
	srsrdl_21.get_all_rotamer_samples( bb5 );
	srsrdl_22.get_all_rotamer_samples( bb5 );
	rsrdl_31.get_all_rotamer_samples( bb5 );
	rsrdl_32.get_all_rotamer_samples( bb5 );
	rsrdl_33.get_all_rotamer_samples( bb5 );
	rsrdl_34.get_all_rotamer_samples( bb5 );
	rsrdl_35.get_all_rotamer_samples( bb5 );
	srsrdl_31.get_all_rotamer_samples( bb5 );
	srsrdl_32.get_all_rotamer_samples( bb5 );
	rsrdl_41.get_all_rotamer_samples( bb5 );
	rsrdl_42.get_all_rotamer_samples( bb5 );
	rsrdl_43.get_all_rotamer_samples( bb5 );
	rsrdl_44.get_all_rotamer_samples( bb5 );
	rsrdl_45.get_all_rotamer_samples( bb5 );
	srsrdl_41.get_all_rotamer_samples( bb5 );
	srsrdl_42.get_all_rotamer_samples( bb5 );
	srsrdl_51.get_all_rotamer_samples( bb5 );
	srsrdl_52.get_all_rotamer_samples( bb5 );

	ChiVector chiv;

	rsrdl_11.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	rsrdl_12.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	rsrdl_13.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	rsrdl_14.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	rsrdl_15.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	srsrdl_11.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	srsrdl_12.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	rsrdl_21.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	rsrdl_22.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	rsrdl_23.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	rsrdl_24.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	rsrdl_25.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	srsrdl_21.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	srsrdl_22.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	rsrdl_31.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	rsrdl_32.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	rsrdl_33.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	rsrdl_34.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	rsrdl_35.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	srsrdl_31.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	srsrdl_32.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	rsrdl_41.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	rsrdl_42.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	rsrdl_43.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	rsrdl_44.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	rsrdl_45.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	srsrdl_41.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	srsrdl_42.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	srsrdl_51.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );
	srsrdl_52.assign_random_rotamer_with_bias( rsd, pose, scratch, numeric::random::rg(), chiv, true );

	rsrdl_11.get_probability_for_rotamer( bb1, 1 );
	rsrdl_12.get_probability_for_rotamer( bb2, 1 );
	rsrdl_13.get_probability_for_rotamer( bb3, 1 );
	rsrdl_14.get_probability_for_rotamer( bb4, 1 );
	rsrdl_15.get_probability_for_rotamer( bb5, 1 );
	srsrdl_11.get_probability_for_rotamer( bb1, 1 );
	srsrdl_12.get_probability_for_rotamer( bb1, 1 );
	rsrdl_21.get_probability_for_rotamer( bb1, 1 );
	rsrdl_22.get_probability_for_rotamer( bb2, 1 );
	rsrdl_23.get_probability_for_rotamer( bb3, 1 );
	rsrdl_24.get_probability_for_rotamer( bb4, 1 );
	rsrdl_25.get_probability_for_rotamer( bb5, 1 );
	srsrdl_21.get_probability_for_rotamer( bb2, 1 );
	srsrdl_22.get_probability_for_rotamer( bb2, 1 );
	rsrdl_31.get_probability_for_rotamer( bb1, 1 );
	rsrdl_32.get_probability_for_rotamer( bb2, 1 );
	rsrdl_33.get_probability_for_rotamer( bb3, 1 );
	rsrdl_34.get_probability_for_rotamer( bb4, 1 );
	rsrdl_35.get_probability_for_rotamer( bb5, 1 );
	srsrdl_31.get_probability_for_rotamer( bb3, 1 );
	srsrdl_32.get_probability_for_rotamer( bb3, 1 );
	rsrdl_41.get_probability_for_rotamer( bb1, 1 );
	rsrdl_42.get_probability_for_rotamer( bb2, 1 );
	rsrdl_43.get_probability_for_rotamer( bb3, 1 );
	rsrdl_44.get_probability_for_rotamer( bb4, 1 );
	rsrdl_45.get_probability_for_rotamer( bb5, 1 );
	srsrdl_41.get_probability_for_rotamer( bb4, 1 );
	srsrdl_42.get_probability_for_rotamer( bb4, 1 );
	srsrdl_51.get_probability_for_rotamer( bb5, 1 );
	srsrdl_52.get_probability_for_rotamer( bb5, 1 );

	rsrdl_11.get_rotamer( bb1, 1 );
	rsrdl_12.get_rotamer( bb2, 1 );
	rsrdl_13.get_rotamer( bb3, 1 );
	rsrdl_14.get_rotamer( bb4, 1 );
	rsrdl_15.get_rotamer( bb5, 1 );
	srsrdl_11.get_rotamer( bb1, 1 );
	srsrdl_12.get_rotamer( bb1, 1 );
	rsrdl_21.get_rotamer( bb1, 1 );
	rsrdl_22.get_rotamer( bb2, 1 );
	rsrdl_23.get_rotamer( bb3, 1 );
	rsrdl_24.get_rotamer( bb4, 1 );
	rsrdl_25.get_rotamer( bb5, 1 );
	srsrdl_21.get_rotamer( bb2, 1 );
	srsrdl_22.get_rotamer( bb2, 1 );
	rsrdl_31.get_rotamer( bb1, 1 );
	rsrdl_32.get_rotamer( bb2, 1 );
	rsrdl_33.get_rotamer( bb3, 1 );
	rsrdl_34.get_rotamer( bb4, 1 );
	rsrdl_35.get_rotamer( bb5, 1 );
	srsrdl_31.get_rotamer( bb3, 1 );
	srsrdl_32.get_rotamer( bb3, 1 );
	rsrdl_41.get_rotamer( bb1, 1 );
	rsrdl_42.get_rotamer( bb2, 1 );
	rsrdl_43.get_rotamer( bb3, 1 );
	rsrdl_44.get_rotamer( bb4, 1 );
	rsrdl_45.get_rotamer( bb5, 1 );
	srsrdl_41.get_rotamer( bb4, 1 );
	srsrdl_42.get_rotamer( bb4, 1 );
	srsrdl_51.get_rotamer( bb5, 1 );
	srsrdl_52.get_rotamer( bb5, 1 );

	utility::fixedsizearray1< Size, 1 > sizevec1;
	utility::fixedsizearray1< Real, 1 > realvec1;
	utility::fixedsizearray1< Size, 2 > sizevec2;
	utility::fixedsizearray1< Real, 2 > realvec2;
	utility::fixedsizearray1< Size, 3 > sizevec3;
	utility::fixedsizearray1< Real, 3 > realvec3;
	utility::fixedsizearray1< Size, 4 > sizevec4;
	utility::fixedsizearray1< Real, 4 > realvec4;
	utility::fixedsizearray1< Size, 5 > sizevec5;
	utility::fixedsizearray1< Real, 5 > realvec5;
	utility::vector1< Real > realvec( 5 );
	
	srsrdl_11.interpolate_nrchi_values( sizevec1, sizevec1, realvec1, 1, realvec );
	srsrdl_12.interpolate_nrchi_values( sizevec1, sizevec1, realvec1, 1, realvec );
	srsrdl_21.interpolate_nrchi_values( sizevec2, sizevec2, realvec2, 1, realvec );
	srsrdl_22.interpolate_nrchi_values( sizevec2, sizevec2, realvec2, 1, realvec );
	srsrdl_31.interpolate_nrchi_values( sizevec3, sizevec3, realvec3, 1, realvec );
	srsrdl_32.interpolate_nrchi_values( sizevec3, sizevec3, realvec3, 1, realvec );
	srsrdl_41.interpolate_nrchi_values( sizevec4, sizevec4, realvec4, 1, realvec );
	srsrdl_42.interpolate_nrchi_values( sizevec4, sizevec4, realvec4, 1, realvec );
	srsrdl_51.interpolate_nrchi_values( sizevec5, sizevec5, realvec5, 1, realvec );
	srsrdl_52.interpolate_nrchi_values( sizevec5, sizevec5, realvec5, 1, realvec );
	
}


} // namespace dunbrack
} // namespace scoring
} // namespace core

