// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/dunbrack/SingleResidueDunbrackLibrary.hh
/// @brief  SingleResidueDunbrackLibrary class
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_pack_dunbrack_SingleResidueDunbrackLibrary_hh
#define INCLUDED_core_pack_dunbrack_SingleResidueDunbrackLibrary_hh

// Unit Headers
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.fwd.hh>

// Package Headers
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.fwd.hh>
#include <core/chemical/AA.hh>

// Utility Headers
#include <utility/assert.hh>
#include <utility/io/izstream.fwd.hh>
#include <utility/io/ozstream.fwd.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>

// Numeric Headers
#include <numeric/numeric.functions.hh>

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace dunbrack {

class SingleResidueDunbrackLibrary : public rotamers::SingleResidueRotamerLibrary {
public:
	typedef chemical::AA AA;

	SingleResidueDunbrackLibrary& operator=(SingleResidueDunbrackLibrary const & other) = delete;
	SingleResidueDunbrackLibrary ( SingleResidueDunbrackLibrary const & ) = delete;
public:

	//std::string who_am_i() override { return "SingleResidueDunbrackLibrary"; }


	/// constants

	/// A good phi and psi values to use for terminal peptide residues where they cannont be completely defined
	static Real const PEPTIDE_NEUTRAL_PHI;
	static Real const PEPTIDE_NEUTRAL_PSI;

	/// A good omega, phi, and psi values to use for terminal peptoid residues where they cannont be completely defined
	static Real const PEPTOID_NEUTRAL_OMG;
	static Real const PEPTOID_NEUTRAL_PHI;
	static Real const PEPTOID_NEUTRAL_PSI;

	/// @brief Precision measures for comparsions.
	static Real const ANGLE_DELTA;
	static Real const PROB_DELTA;
	static Real const ENERGY_DELTA;
	static Real const COEF_DELTA; // intepolation coefficients

public:

	/// c-tor
	SingleResidueDunbrackLibrary(
		chemical::ResidueType const & rt,
		Size const n_rotameric_chi,
		bool dun02,
		bool use_bicubic,
		bool dun_entropy_correction,
		core::Real prob_buried, // 0.98
		core::Real prob_nonburied // 0.95
	);

public:

	virtual void write_to_binary( utility::io::ozstream & out ) const;
	virtual void read_from_binary( utility::io::izstream & in );

	/// @brief Return all of the rotamer sample data given a particular phi/psi.
	/// For N-terminus residues, hand in the phi value SingleResidueDunbrackLibrary::PHI_NEUTRAL and
	/// for C-terminus residues, hand in the psi value SingleResidueDunbrackLibrary::PSI_NEUTRAL.
	/// The returned samples should be in semi-decrasing order by probability; semi, because the
	/// rotamers are constructed in sorted order by their probability in the lower phi-psi bin that
	/// the input phi/psi perscribes.
	/*virtual
	utility::vector1< DunbrackRotamerSampleData >
	get_all_rotamer_samples(
	Real phi,
	Real psi
	) const = 0;*/

	virtual
	utility::vector1< DunbrackRotamerSampleData >
	get_all_rotamer_samples(
		Real5 bbs
	) const = 0;

	/// @brief Return the probability for a particular rotamer where rotamers are
	/// indexed in order of decreasing probability (or something very close to
	/// decreasing probability).
	virtual
	Real
	get_probability_for_rotamer(
		Real phi,
		Real psi,
		Size rot_ind
	) const = 0;

	/*virtual
	Real
	get_probability_for_rotamer(
	Real5 bbs,
	Size rot_ind
	) const = 0;*/

	virtual
	DunbrackRotamerSampleData
	get_rotamer(
		Real phi,
		Real psi,
		Size rot_ind
	) const = 0;

	/*virtual
	DunbrackRotamerSampleData
	get_rotamer(
	utility::vector1< Real > bbs,
	Size rot_ind
	) const = 0;*/

	virtual
	Real
	get_phi_from_rsd(
		conformation::Residue const & rsd
	) const = 0;

	virtual
	Real
	get_psi_from_rsd(
		conformation::Residue const & rsd
	) const = 0;

	//virtual
	//Real
	//get_IV_from_rsd(
	// conformation::Residue const & rsd,
	// Size bbn
	//) const = 0;
	//virtual
	//utility::vector1< Real >
	//get_IVs_from_rsd(
	// conformation::Residue const & rsd
	//) const = 0;

public:
	/// Virtual functions the derived classes must implement

	/// @brief Derived classes should invoke base class function as well.
	virtual Size memory_usage_in_bytes() const;

	/// @brief The number of chi represented by the library.
	virtual
	Size nchi() const = 0;

	/// @brief the number of backbone dihedrals represented by the library
	virtual
	Size nbb() const = 0;

	virtual
	Size n_rotamer_bins() const = 0;

	/// @brief Tell the base class the number of chi bins for each rotameric
	/// chi dimension
	void
	set_n_chi_bins( utility::vector1< Size > const & );

protected:
	/// Read access for the derived class
	bool dun02() const { return dun02_; }

	/// Read access for the derived class
	bool use_bicubic() const { return use_bicubic_; }

	/// Read access for the derived class
	bool dun_entropy_correction() const { return dun_entropy_correction_; }

	/// Worker functions available to the derived classes

	virtual Size memory_usage_static() const = 0;
	virtual Size memory_usage_dynamic() const;

	/// @brief Read access to the n_chi_bins_ vector
	utility::vector1< Size > const &
	n_chi_bins() const {
		return n_chi_bins_;
	}

	/// @brief The base class needs to be informed about which rotamer wells
	/// exist in order to create the rotwell to packed rot conversion data.
	/// set_chi_nbins must be called first.
	void
	mark_rotwell_exists( utility::vector1< Size > const & rotwell );

	/// @brief After the derived class has marked all the rotwells that do exist,
	/// the base class will create the rotwell to packerot conversion data.
	void
	declare_all_existing_rotwells_encountered();

	/// @brief The number of existing rotamers
	Size
	n_packed_rots() const {
		return n_packed_rots_;
	}

	/// @brief The number of possible rotamers -- product of the chi_nbins_ array
	Size
	n_possible_rots() const {
		return n_possible_rots_;
	}

public:

	/// @brief Convert a vector of chi angles (degrees) into a integer vector of rotamer wells.
	/// Derived class should be consistent, but may be arbitrary in how wells divide angle space.
	virtual
	void
	get_rotamer_from_chi(
		ChiVector const & chi,
		RotVector & rot ) const = 0;

public:
	/// Conversion functions

	/// @brief Convert from the rotamer bin indices for each chi to the
	/// (non-compact) "rotamer number"
	Size
	rotwell_2_rotno( utility::vector1< Size > const & rotwell ) const;

	/// @brief Convert from the rotamer bin indices for each chi to the
	/// (non-compact) "rotamer number"
	Size
	rotwell_2_rotno( Size4 const & rotwell ) const;

	/// @brief Convert from the rotamer number to the compacted
	/// "packed rotamer number".  Returns 0 if rotno has no corresponding packed rotno.
	Size
	rotno_2_packed_rotno( Size const rotno ) const;

	/// @brief Convert from the rotamer bin indices for each chi to the
	/// compacted "packed rotamer number." Returns 0 if rotwell has no corresponding packed rotno
	Size
	rotwell_2_packed_rotno( utility::vector1< Size > const & rotwell ) const;

	/// @brief Convert from the rotamer bin indices for each chi to the
	/// compacted "packed rotamer number." Returns 0 if rotwell has no corresponding packed rotno
	Size
	rotwell_2_packed_rotno( Size4 const & rotwell ) const;

	/// @brief Convert from the packed rotamer number to the rotamer well
	void
	packed_rotno_2_rotwell( Size const packed_rotno, utility::vector1< Size > & rotwell ) const;

	void
	packed_rotno_2_rotwell(
		Size const packed_rotno,
		Size4 & rotwell
	) const;

	utility::vector1< Size > const &
	packed_rotno_2_rotwell( Size const packed_rotno ) const;

	/// @brief Convert from the rotamer number to the rotamer well
	void
	rotno_2_rotwell( Size const rotno, utility::vector1< Size > & rotwell ) const;

	/// @brief, Turns out, when non-rotameric chi are taken out of the picture,
	/// all remaining chi are binned the same way, except proline. Valid only for
	/// Dun10 libraries.  For D-amino acids, chi must be inverted before passing
	/// to this function.
	inline
	Size
	bin_rotameric_chi(
		Real chi,
		Size which_chi
	) const {

		// I'm removing this debug_assert because it is enforced earlier -- in the
		// function that calls it -- and furthermore because this assert NEEDS to have
		// a more complex condition: it has to accept NCAA rotlibs. So we need to be
		// able to compare aa against last_canonical or last_D or something, and here
		// we don't have aa with us.
		//debug_assert( ! dun02_ );
		//debug_assert( -180.0 <= chi && chi <= 180.0 );
		if ( !(-180.0 <= chi && chi <= 180.0 ) ) {
			throw CREATE_EXCEPTION(utility::excn::RangeError,
				"chi angle must be between -180 and 180: "+utility::to_string(chi));
		}

		if ( aa_ == chemical::aa_pro || aa_ == chemical::aa_dpr /*D-proline*/ || aa_ == chemical::aa_b3p || aa_ == chemical::ou3_pro  ) {
			if ( which_chi == 1 ) {
				if ( chi > 0 ) { return 1; }
				else { return 2; }
			} else {
				return 1;
			}
		}

		if ( ( chi >= 0.0 ) && ( chi <= 120.0 ) ) { return 1; }
		else if ( std::abs(chi) >= 120.0 ) { return 2; }
		else /*if ( ( chi_i >= -120.0 ) && ( chi_i <= 0.0 ) )*/ { return 3; }
	}

	inline
	void bin_angle(
		Real const angle_start,
		Real const angle_step,
		Real const ASSERT_ONLY( angle_range ),
		Size const nbins,
		Real const ang,
		Size & bin_lower,
		Size & bin_upper,
		Real & angle_alpha
	) const {
		/// very, very rarely, periodic_range( angle, 360 ) will return 180 instead of -180.
		/// though it is supposed to return values in the range [-180, 180).
		debug_assert( angle_start <= ang && ang <= angle_start + angle_range );
		// temp -- by merge-time we won't call this function at all.
		//debug_assert( std::abs( nbins * angle_step - angle_range ) < 1e-15 );

		Real real_bin_lower = ( ang - angle_start ) / angle_step;
		Size bin_prev = static_cast< Size > ( real_bin_lower );
		bin_lower = 1 + numeric::mod( bin_prev, nbins );
		bin_upper = numeric::mod( bin_lower, nbins ) + 1;
		angle_alpha = ( (ang - angle_start ) - ( bin_prev * angle_step ) ) / angle_step;
	}

	inline
	void bin_angle(
		Real const angle_start,
		utility::vector1< core::Size > const & bin_equivs,
		Real const angle_step,
		Real const ASSERT_ONLY( angle_range ),
		Size const nbins,
		Real const ang,
		Size & bin_lower,
		Size & bin_upper,
		Real & angle_alpha
	) const {
		/// very, very rarely, periodic_range( angle, 360 ) will return 180 instead of -180.
		/// though it is supposed to return values in the range [-180, 180).
		debug_assert( angle_start <= ang && ang <= angle_start + angle_range );
		// temp -- by merge-time we won't call this function at all.
		//debug_assert( std::abs( nbins * angle_step - angle_range ) < 1e-15 );

		Real real_bin_lower = ( ang - angle_start ) / angle_step;

		Size bin_prev_raw = static_cast< Size > ( real_bin_lower );
		Size bin_lower_raw = 1 + numeric::mod( bin_prev_raw, nbins );
		Size bin_upper_raw = numeric::mod( bin_lower_raw, nbins ) + 1;

		bin_lower = bin_equivs[ bin_lower_raw ];
		// In odd cases if bin_prev_raw is zero we need to wrap. Or vice-versa. (It has to be on [0, 35].)
		Size bin_prev = bin_prev_raw == 0 ? 0 /*numeric::mod( bin_equivs[ nbins ], nbins )*/ : bin_equivs[ bin_prev_raw ];
		bin_upper = bin_equivs[ bin_upper_raw ];


		if ( bin_prev == bin_lower || bin_lower == bin_upper ) { //== bin_equivs[ bin_prev_raw - 1 ] ) {
			angle_alpha = 0;
		} else {
			// raw? bc if we are here raw may be equivalent... SLASH BETTER.
			angle_alpha = ( ( ang - angle_start ) - ( bin_prev_raw * angle_step ) ) / angle_step;
		}

	}


	/*
	/// @brief This is not the right place for this code, but the numeric interpolation library
	/// uselessly indexes by 0 and the basic functions aren't inlined...
	inline
	void bin_angle(
	Real const angle_start,// same as bin_lowers[1] - small_angle_step
	utility::vector1< Real > const & bin_lowers, // could be nonuniform
	Real const small_angle_step, // when stuff isn't being skipped
	Real const angle_range,
	Size const nbins,
	Real const ang,
	Size & bin_lower,
	Size & bin_upper,
	Real & angle_alpha
	) const {
	/// very, very rarely, periodic_range( angle, 360 ) will return 180 instead of -180.
	/// though it is supposed to return values in the range [-180, 180).
	// this just isn't true anymore
	//debug_assert( angle_start <= ang && ang <= angle_start + angle_range );
	// angle_step is now nonuniform. You would have to add up all the differences between
	// consecutive bins. angle_range now cannot be provided; you would have to take the difference
	// between the first and last bin_lower, plus bin_step. this is equivalent to the above
	// so this assert would pass by definition.
	//debug_assert( std::abs( nbins * angle_step - angle_range ) < 1e-15 );
	//Real angle_start = bin_lowers[1] - small_angle_step;
	//Real angle_range = bin_lowers[ bin_lowers.size() ] - bin_lowers[ 1 ] + small_angle_step;
	Real real_ang = ang - angle_start;//numeric::mod(ang - angle_start, angle_range);
	// We want to return real_ang on (0, 360]
	if ( real_ang <= 0 ) real_ang += angle_range;

	std::cout << "angle_start angle_range real_ang " << angle_start << " " << angle_range << " " << real_ang << std::endl;
	for ( bin_upper = 1; bin_upper <= bin_lowers.size(); ++bin_upper ) {
	if ( bin_lowers[ bin_upper ] - angle_start > real_ang ) break;
	}
	//if ( bin_lowers[ bin_upper ] - angle_start == real_ang ) ++bin_upper;
	//std::cout << "bin_upper " << bin_upper << std::endl;
	//if ( bin_lowers[ bin_upper - 1 ] == real_ang ) bin_lower = bin_upper - 1;

	// These bin_uppers appear to be consistently 1 too small.
	// The issue is that we're basically solving for bin_prev+1 and bin_prev?
	++bin_upper;

	bin_upper = bin_upper % nbins;
	if ( bin_upper == 0 ) bin_upper = nbins;

	if ( bin_upper == 1 && real_ang != 0 ) {
	bin_lower = nbins;
	} else if ( bin_upper == 1 && real_ang == 0 ) {
	bin_lower = 1; bin_upper = 2;
	} else  {
	bin_lower = bin_upper - 1;
	}

	std::cout << "bin_upper " << bin_upper << " " << bin_lowers[ bin_upper ] - angle_start
	<< " bin_lower " << bin_lower << " " << bin_lowers[ bin_lower ] - angle_start << std::endl;
	//if ( angle_alpha > 1 ) { // sign you're in the middle of a 'skip'
	if ( bin_lowers[ bin_lower ] - angle_start - small_angle_step > real_ang || std::abs( bin_lowers[ bin_upper ] - bin_lowers[ bin_lower ] ) > small_angle_step ) {
	// 'snap' to lower or upper, respectively.
	// A skip that happens around the wrap (300 to 60, for example) would
	// fail this condition, I think.
	if ( ( real_ang - ( bin_lowers[ bin_lower ] - angle_start - small_angle_step ) ) / std::abs( bin_lowers[ bin_upper ] - bin_lowers[ bin_lower ] ) >= 0.5 ) {
	bin_lower = bin_upper;
	angle_alpha = 0;
	} else {
	bin_upper = bin_lower;
	angle_alpha = 1;
	}
	//angle_alpha = 0;
	} else {
	angle_alpha = ( real_ang - ( bin_lowers[ bin_lower ] - angle_start - small_angle_step ) ) / small_angle_step;
	}
	}
	*/

public:
	/// @brief The amino acid this library is representing
	AA
	aa() const {
		return aa_;
	}

	/// @brief When creating rotamer, what position in the CDF should one build until?
	/// Unlikely rotamers ( < 0.5 %) are numerous, but are very infrequently useful.
	Real
	probability_to_accumulate_while_building_rotamers( bool buried ) const;

	/// @brief setters for accumulation probability cutoff (to support externally-controlled option dependence)
	void prob_to_accumulate( Real, Real );
	void prob_to_accumulate_buried( Real );
	void prob_to_accumulate_nonburied( Real );

	/// @brief Extract the number of rotamer bins (vector of bin counts for each chi) from the ResidueType
	/// and store it in rot.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu) to make this general and to remove hard-coded
	/// information about canonical amino acids.
	static void n_rotamer_bins_for_aa(
		chemical::ResidueType const & rt,
		RotVector & rot,
		bool const dun02=false
	);

	/// @details The number of wells for canonical amino acids defined for the 08 library; includes
	/// the number of wells for the semi-rotameric chi (arbitrarily chosen).
	static void n_rotamer_bins_for_aa(
		chemical::AA const aa,
		RotVector & rot
	);

	/// @details For the rotameric chi -- not all chi are rotameric
	static void n_rotameric_bins_for_aa(
		chemical::AA const aa,
		RotVector & rot,
		bool dun02
	);

	/// @details To continue supporting the 2002 Dunbrack library, we
	/// need to preserve the rotamer well definitions for canonical amino acids that it made.
	static void n_rotamer_bins_for_aa_02(
		chemical::AA const aa,
		RotVector & rot
	);

	/// @brief Comparison operator, mainly intended to use in ASCII/binary comparsion tests
	/// Values tested should parallel those used in the read_from_binary() function.
	virtual
	bool
	operator ==( rotamers::SingleResidueRotamerLibrary const & ) const;

private:

	/// @brief This function forces the instantiation of virtual templated methods in the derived classes.
	/// Functions like this one are necessary when combining polymorphism and templates.  Though
	/// these functions must be compiled, they need never be called. Do not call this function.
	void hokey_template_workaround();


private:
	////////////////////////
	// Options settings
	bool const dun02_; // Are we using the 2002 definitions for rotamer wells?
	bool use_bicubic_; // Are we using a bicubic interpolation?
	bool dun_entropy_correction_; // Are we applying the entropy correction?

	/////////////
	/// data
	AA const aa_;
	Size const n_rotameric_chi_;
	utility::vector1< Size > n_chi_bins_;
	utility::vector1< Size > n_chi_products_; // n_chi_products_[ i ] = prod( j in i+1 to nchi, n_chi_bins_[ j ] );

	Size n_packed_rots_;
	Size n_possible_rots_; // prod( i in 1 to nchi, n_chi_bins_[ i ] );

	Real prob_to_accumulate_buried_, prob_to_accumulate_nonburied_;

	utility::vector1< bool > rotwell_exists_;
	bool packed_rotno_conversion_data_current_;

	utility::vector1< Size > rotno_2_packed_rotno_;
	utility::vector1< Size > packed_rotno_2_rotno_;
	utility::vector1< utility::vector1< Size > > packed_rotno_2_rotwell_;

};


} // dunbrack
} // pack
} // core

#endif // INCLUDED_core_pack_dunbrack_SingleResidueDunbrackLibrary_HH
