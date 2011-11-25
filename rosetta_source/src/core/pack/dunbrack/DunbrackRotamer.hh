// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/dunbrack/DunbrackRotamer.hh
/// @brief  Fixed-sized dunbrack rotamer object
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_core_pack_dunbrack_DunbrackRotamer_hh
#define INCLUDED_core_pack_dunbrack_DunbrackRotamer_hh

// Unit headers
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>

// Package headers
#include <core/pack/dunbrack/ChiSet.fwd.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.fwd.hh>
// AUTO-REMOVED #include <core/scoring/types.hh>
#include <core/types.hh>

// Project headers
#include <core/chemical/ResidueType.fwd.hh>
// AUTO-REMOVED #include <core/conformation/Residue.fwd.hh>
// AUTO-REMOVED #include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/RotamerSampleOptions.hh>

// Numeric headers
#include <numeric/random/random.fwd.hh>

// Utility headers
#include <utility/fixedsizearray1.hh>
#include <utility/pointer/ReferenceCount.hh>
// AUTO-REMOVED #include <utility/vector1.hh>
// AUTO-REMOVED #include <utility/exit.hh>

#include <utility/vector1_bool.hh>



namespace core {
namespace pack {
namespace dunbrack {

/// @brief A class who's size is known at compile time.  A vector of these objects
/// will occupy a perfectly contiguous region of memory, with no extra space
/// allocated for pointers.  These objects may be rapidly allocated and deallocated
/// on the stack.
///
/// @details S for Size, P for Precision

template < Size S, class P >
class DunbrackRotamerMeanSD {
public:
	DunbrackRotamerMeanSD() :
		rotamer_probability_( P( 0.0 ) )
	{
		for ( Size ii = 1; ii <= S; ++ii ) {
			chi_mean_[ ii ] = P( 0.0 );
			chi_sd_[ ii ]   = P( 0.0 );
		}
	}

	DunbrackRotamerMeanSD(
		typename utility::vector1< P > const & chimean_in,
		typename utility::vector1< P > const & chisd_in,
		P const prob_in
	) :
		rotamer_probability_( prob_in )
	{
		for ( Size ii = 1; ii <= S; ++ii ) {
			chi_mean_[ ii ] = chimean_in[ ii ];
			chi_sd_  [ ii ] = chisd_in  [ ii ];
		}
	}

	P chi_mean ( Size which_chi ) const {
		return chi_mean_[ which_chi ];
	}

	P chi_sd ( Size which_chi ) const {
		return chi_sd_[ which_chi ];
	}

	P
	rotamer_probability() const {
		return rotamer_probability_;
	}

	P &
	chi_mean( Size which_chi ) {
		return chi_mean_[ which_chi ];
	}

	P &
	chi_sd( Size which_chi ) {
		return chi_sd_[ which_chi ];
	}

	P &
	rotamer_probability() {
		return rotamer_probability_;
	}

	void chi_mean ( Size which_chi, P chi_mean_in ) {
		chi_mean_[ which_chi ] = chi_mean_in;
	}

	void chi_sd ( Size which_chi, P chi_sd_in ) {
		chi_sd_[ which_chi ] = chi_sd_in;
	}

	void rotamer_probability( P rotprob_in ) {
		rotamer_probability_ = rotprob_in;
	}

private:
	utility::fixedsizearray1< P, S > chi_mean_;
	utility::fixedsizearray1< P, S > chi_sd_;
	P rotamer_probability_;
};

template < Size S, class P >
class DunbrackRotamer : public DunbrackRotamerMeanSD< S, P > {
public:
	typedef DunbrackRotamerMeanSD< S, P > parent;

public:

	DunbrackRotamer() :
		parent()
	{
		for ( Size ii = 1; ii <= S; ++ii ) {
			rotwell_[ ii ] = 0;
		}
	}

	DunbrackRotamer(
		typename utility::vector1< P > const & chimean_in,
		typename utility::vector1< P > const & chisd_in,
		P const prob_in,
		typename utility::vector1< Size > const & rotwell_in
	) :
		parent( chimean_in, chisd_in, prob_in )
	{
		for ( Size ii = 1; ii <= S; ++ii ) {
			rotwell_[ ii ] = rotwell_in[ ii ];
		}
	}

	Size rotwell( Size which_chi ) const {
		return rotwell_[ which_chi ];
	}


	Size &
	rotwell( Size which_chi ) {
		return rotwell_[ which_chi ];
	}


	void rotwell( Size which_chi, Size rotwell_in ) {
		rotwell_[ which_chi ] = rotwell_in;
	}

private:
	utility::fixedsizearray1< Size, S > rotwell_;
};


template < Size S, class P >
class PackedDunbrackRotamer : public DunbrackRotamerMeanSD< S, P > {
public:
	typedef DunbrackRotamerMeanSD< S, P > parent;

public:

	PackedDunbrackRotamer() :
		parent(),
		packed_rotno_( 0 )
	{}

	PackedDunbrackRotamer(
		typename utility::vector1< P > const & chimean_in,
		typename utility::vector1< P > const & chisd_in,
		P const prob_in,
		Size const packed_rotno_in
	) :
		parent( chimean_in, chisd_in, prob_in ),
		packed_rotno_( packed_rotno_in )
	{}

	PackedDunbrackRotamer(

		DunbrackRotamer< S, P > const & sibling,
		Size const packed_rotno_in
	) :
		parent( sibling ),
		packed_rotno_( packed_rotno_in )
	{}


	Size &
	packed_rotno() {
		return packed_rotno_;
	}

	void
	packed_rotno( Size packed_rotno_in ) {
		packed_rotno_ = packed_rotno_in;
	}

	Size
	packed_rotno() const {
		return packed_rotno_;
	}


private:

	Size packed_rotno_;

};

class DunbrackRotamerSampleData
{
public:
	DunbrackRotamerSampleData();
	DunbrackRotamerSampleData( bool is_nrchi_sample );
	~DunbrackRotamerSampleData();

public:
	/// Setters

	void set_nrchi_sample( bool setting );
	void set_nchi( Size nchi );
	void set_rotwell(  Size chi_index, Size rotwell );
	void set_rotwell( utility::vector1< Size > const & );
	void set_chi_mean( Size chi_index, Real mean    );
	void set_chi_sd(   Size chi_index, Real sd      );
	void set_prob( Real probability );

	void set_nrchi_lower_boundary( Real low );
	void set_nrchi_upper_boundary( Real high );
	void set_nrchi_probability( Real nrchi_prob );

	bool  nrchi_sample() const { return nrchi_sample_; }
	bool  chi_is_nonrotameric( Size chi ) const { return nrchi_sample_ && chi == nchi_; }
	Size  nchi() const { return nchi_; }
	Size4 const & rot_well() const { return rot_well_; }
	Real4 const & chi_mean() const { return chi_mean_; }
	Real4 const & chi_sd() const { return chi_sd_; }
	Real  probability() const { return probability_; }

	Real  nrchi_lower_boundary() const { assert( nrchi_sample_ ); return nrchi_lower_boundary_; }
	Real  nrchi_upper_boundary() const { assert( nrchi_sample_ ); return nrchi_upper_boundary_; }
	Real  nrchi_probability() const { assert( nrchi_sample_ ); return nrchi_probability_; }

	void  assign_random_chi( utility::vector1< Real > & chi_angles, numeric::random::RandomGenerator & RG, core::Real factor=1.0) const;
	Real  chi_probability( utility::vector1< Real > const & chi_angles, core::Real factor=1.0 ) const;

private:
	bool  nrchi_sample_; /// Does this sample describe a semi-rotameric residue?
	Size  nchi_;         /// The number of chi values the library describes; Between 1 & 4.
	Size4 rot_well_;     /// The integer description of the rotamer conformation
	Real4 chi_mean_;     /// The chi angle at the center of the distribution
	Real4 chi_sd_;       /// Standard deviation from the center for each chi.
	Real  probability_;  /// The probability of finding the rotamer in this rot_well.

	Real  nrchi_lower_boundary_; /// Requires nrchi_sample_.  The lower boundary of the non-rotameric chi well, in degrees.
	Real  nrchi_upper_boundary_; /// Requires nrchi_sample_.  The upper boundary of the non-rotameric chi well, in degrees.
	Real  nrchi_probability_;    /// Requires nrchi_sample_.  The probability that the non-rotameric chi ends up in this well given its rotameric-chi assignment.
};


/// @brief a simple class for passing data around in virtual function
/// calls of the rotamer creating process.  Derived classes will be simple
/// containers for interpolated rotameric data that 1) has to be available
/// to the derived class when building rotamers and 2) cannot be stored as
/// member data in the derived class in a thread-safe manner.  Derived classes
/// of the RotamerBuildingData can be declared on the stack, passed into
/// the RotamericSingleResidueDunbrackLibrary::build_chi_sets function,
/// and then in the (virtual) chisamples_for_rotamer function, the derived classes
/// may be downcast.
class RotamerBuildingData : public utility::pointer::ReferenceCount
{
public:
	virtual ~RotamerBuildingData() = 0;
};



/// Should this be here?

void
expand_proton_chi(
	pack::task::ExtraRotSample ex_samp_level,
	chemical::ResidueTypeCAP concrete_residue,
	Size proton_chi,
	utility::vector1< ChiSetOP > & chi_set_vector
);

void interpolate_rotamers(
	DunbrackRotamer< FOUR > const & rot00,
	DunbrackRotamer< FOUR > const & rot10,
	DunbrackRotamer< FOUR > const & rot01,
	DunbrackRotamer< FOUR > const & rot11,
	Real phi_err, Real psi_err, Real binrange,
	Size nchi_aa,
	DunbrackRotamer< FOUR, Real > & interpolated_rotamer
);

template < Size S >
DunbrackRotamer< S, Real >
increase_rotamer_precision(
	DunbrackRotamer< S, DunbrackReal > const & original_rotamer
)
{
	DunbrackRotamer< S, Real > new_rotamer;
	for ( Size ii = 1; ii <= S; ++ii ) {
		new_rotamer.chi_mean( ii ) = static_cast< Real > ( original_rotamer.chi_mean( ii ) );
		new_rotamer.chi_sd( ii )   = static_cast< Real > ( original_rotamer.chi_sd( ii ) );
		new_rotamer.rotwell( ii )   = original_rotamer.rotwell( ii );
	}
	new_rotamer.rotamer_probability() = static_cast< Real > ( original_rotamer.rotamer_probability() );
	return new_rotamer;
}


} //dunbrack
} //scoring
} // core

#endif
