// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/Ramachandran.hh
/// @brief  Ramachandran potential class delcaration
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_scoring_Ramachandran_hh
#define INCLUDED_core_scoring_Ramachandran_hh

// Unit Headers
#include <core/scoring/Ramachandran.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/ppo_torsion_bin.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/interpolation/spline/Bicubic_spline.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray4D.hh>

// C++ Headers
#include <map>


namespace core {
namespace scoring {

enum Rama_Table_Type {
	//When adding an effect to this enum:
	// 1. Add its name to the get_ramatable_name_by_type() function.
	// 2. Add an option to basic/options/options_rosetta.py for the file that corresponds to this rama table.
	// 3. Add the option name and this enum to get_custom_rama_table_filename().
	// 4. If you're adding a D-amino acid type, be sure to update all the functions that check for flat_d_aa_ramatable and substitute the appropriate L-amino acid tables.
	flat_l_aa_ramatable=1,
	flat_d_aa_ramatable,
	flat_symm_dl_aa_ramatable,
	flat_symm_gly_ramatable,
	flat_symm_pro_ramatable,
	flat_l_aa_ramatable_stringent,
	flat_d_aa_ramatable_stringent,
	flat_symm_dl_aa_ramatable_stringent,
	flat_symm_gly_ramatable_stringent,
	flat_symm_pro_ramatable_stringent,
	unknown_ramatable_type, //Keep this second-to-last
	end_of_ramatable_type_list=unknown_ramatable_type //Keep this last
};

class Ramachandran : public utility::pointer::ReferenceCount
{
public:
	typedef pose::Pose Pose;
	typedef chemical::AA AA;

public:
	Ramachandran();

	Ramachandran(
		std::string const & rama_map_filename,
		bool use_bicubic_interpolation
	);

	virtual ~Ramachandran() ; // auto-removing definition from header{}

	/// @brief Given a custom Rama table type, get the filename of the database file containing
	/// the data.
	/// @details Accesses the options system for this information.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	std::string
	get_custom_rama_table_filename (
		Rama_Table_Type const type
	) const;

	/// @brief Given a Rama_Table_Type name, return the type.
	/// @details Calls get_rama_table_name_by_type().
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	Rama_Table_Type
	get_ramatable_type_by_name (
		std::string const &name
	) const;

	/// @brief Given a Rama_Table_Type, return the name.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	std::string
	get_ramatable_name_by_type (
		Rama_Table_Type const type
	) const;

	bool
	is_canonical_d_aminoacid(
		AA const res_aa
	) const;

	AA
	get_l_equivalent(
		AA const d_aa
	) const;

	Real
	eval_rama_score_residue(
		conformation::Residue const & res
	) const;

	Real
	eval_rama_score_residue(
		AA const res_aa,
		Real const phi,
		Real const psi
	) const;

	void
	eval_rama_score_residue(
		conformation::Residue const & res,
		Real & rama,
		Real & drama_dphi,
		Real & drama_dpsi
	) const;

	void
	eval_rama_score_residue(
		AA const res_aa,
		Real const phi,
		Real const psi,
		Real & rama,
		Real & drama_dphi,
		Real & drama_dpsi
	) const;

	void
	eval_rama_score_residue(
		bool use_bicubic_interpolation,
		bool rama_not_squared,
		AA const res_aa,
		Real const phi,
		Real const psi,
		Real & rama,
		Real & drama_dphi,
		Real & drama_dpsi
	) const;

	/// @brief Function to evaluate the rama score for residues whose connection partners are not necessarily adjacent
	/// in linear sequence (e.g. for backbone-cyclized peptides).  Note that this assumes that rama is being used only
	/// for scoring alpha-amino acids (L- or D-).
	/// @author Vikram K. Mulligan
	void
	eval_rama_score_residue_nonstandard_connection(
		core::pose::Pose const & mypose,
		conformation::Residue const & res,
		Real & rama,
		Real & drama_dphi,
		Real & drama_dpsi
	) const;

	/// @brief Generate random values for phi and psi, biased by the Ramachandran plot of a particular amino acid.
	/// @param[in] res_aa The amino acid in question.
	/// @param[out] phi Output phi value.
	/// @param[out] psi Output psi value.
	void
	random_phipsi_from_rama(
		AA const res_aa,
		Real & phi,
		Real & psi
	) const;

	/// @brief Return a phi/psi pair picked uniformly from the regions of rama
	/// space with nonzero weight.  Sampling with this method will not give a
	/// rama distribution; it will give a flat distribution in only the allowed
	/// regions of rama space.
	void
	uniform_phipsi_from_allowed_rama(
		AA const res_aa,
		Real & phi,
		Real & psi
	) const;

	/// @brief Return true if the given phi/psi pair is in the allowed space
	/// sampled by uniform_phipsi_from_allowed_rama.
	bool
	phipsi_in_allowed_rama(
		AA const res_aa,
		Real phi,
		Real psi
	) const;

	/// @brief Return false if the given phi/psi pair is in the allowed space
	/// sampled by uniform_phipsi_from_allowed_rama.
	bool
	phipsi_in_forbidden_rama(
		AA const res_aa,
		Real phi,
		Real psi
	) const;


	/// @brief functions for torsion-bin specific but otherwise random phi/psi angles
	/// @author Amelie Stein
	void
	random_phipsi_from_rama_by_torsion_bin(
		AA res_aa,
		Real & phi,
		Real & psi,
		conformation::ppo_torsion_bin torsion_bin
	) const;

	void
	get_entries_per_torsion_bin( AA const res_aa, std::map< conformation::ppo_torsion_bin, core::Size > & tb_frequencies ) const;


	///////////////////////////////
	// unused??
	void
	eval_rama_score_all(
		Pose & pose,
		ScoreFunction const & scorefxn
	) const;

	void
	write_rama_score_all(
		Pose const & pose
	) const;

	//Real get_rama_score_residue_deriv( int res, Pose const & a_pose, ProteinTorsion torsion ) const;
	//void eval_procheck_rama( Pose const & a_pose,
	// Real & favorable, Real & allowed, Real & generous ) const;

	/// @brief Function to do a quick check that the upper connection is seqpos+1 and the lower connection is seqpos-1.  Returns true if this is so, false otherwise.
	/// @author Vikram K. Mulligan
	bool is_normally_connected ( conformation::Residue const & res ) const;

	Size n_phi_bins() const;
	Size n_psi_bins() const;

	Real rama_probability( core::chemical::AA aa, Real phi, Real psi ) const;
	Real minimum_sampling_probability() const;

	utility::vector1< Real > const & cdf_for_aa( core::chemical::AA aa ) const;

	utility::vector1< Real > const & cdf_for_aa_for_torsion_bin(
		chemical::AA aa,
		conformation::ppo_torsion_bin
	) const;

	/// @brief Pick a random phi, psi value from a custom Rama table.
	/// @details The custom Rama table is lazily loaded, so this function
	/// is necessarily non-const.  By default, only the 20 canonical Rama
	/// tables are loaded.
	/// @param[in] type The type of custom rama table (an enum value).
	/// @param[out] phi Randomly-drawn phi value, biased by the custom rama
	/// table.
	/// @param[out] psi Randomly-drawn psi value, biased by the custom rama
	/// table.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void
	draw_random_phi_psi_from_extra_cdf(
		Rama_Table_Type const type,
		Real & phi,
		Real & psi
	);

private:

	void read_rama(
		std::string const & rama_map_filename,
		bool use_bicubic_interpolation);

	/// @brief Load a custom Ramachandran table, in addition to the 20x3 standard
	/// ones that are always loaded.
	/// @detailed Intended for sampling with alternative Ramachandran distributions.  Custom
	/// tables are lazily loaded so as not to add to total Rosetta memory footprint.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void load_custom_rama_table( Rama_Table_Type const type );

	/// @brief If the -symmetric_gly_tables option is used, symmetrize the aa_gly table.
	/// @details By default, the gly table is asymmetric because it is based on statistics from the PDB (which disproportionately put glycine
	/// in the D-amino acid region of Ramachandran space).  However, the intrinsic propensities of glycine make it equally inclined to favour
	/// right- or left-handed conformation.  (Glycine is achrial, and can't have a preference.)  Must be called AFTER gly table load, but prior
	/// to bicubic interpolation setup.  Note that symmetrization is based on the probability table, and carries over to the energy table; the
	/// counts table is left as-is (asymmetric).
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void symmetrize_gly_table();

	void read_rama_map_file ( utility::io::izstream * iunit );
	//MaximCode
	void read_rama_map_file_shapovalov ( utility::io::izstream * iunit );
	void initialize_rama_sampling_tables();

	void init_rama_sampling_table( conformation::ppo_torsion_bin torsion_bin );
	void init_uniform_sampling_table();

	/// @brief Pick a random phi and psi value given a cumulative distribution function.
	/// @param[in] cdf The cumulative distribution function.
	/// @param[out] phi The output phi value.
	/// @param[out] psi The output psi value.
	/// @details A random bin from the cumulative distribution function is chosen first.  Then uniform
	/// randomness is added within the chosen bin to produce the final phi and psi values.
	void
	draw_random_phi_psi_from_cdf(
		utility::vector1< Real > const & cdf,
		Real & phi,
		Real & psi
	) const;

	/// @brief Has a particular custom Rama probability table been loaded?
	/// @details Custom Rama tables are lazily loaded.
	inline bool has_custom_rama_probability_table( Rama_Table_Type const type ) const {
		if ( type == flat_d_aa_ramatable ) return (extra_ram_probabil_.count(flat_l_aa_ramatable)!=0); //Special case: D-amino acids use the L-table, inverted, to save RAM.
		return (extra_ram_probabil_.count(type)!=0);
	}

	/// @brief Has a particular custom Rama count table been loaded?
	/// @details Custom Rama tables are lazily loaded.
	inline bool has_custom_rama_count_table( Rama_Table_Type const type ) const {
		if ( type == flat_d_aa_ramatable ) return (extra_ram_counts_.count(flat_l_aa_ramatable)!=0); //Special case: D-amino acids use the L-table, inverted, to save RAM.
		return (extra_ram_counts_.count(type)!=0);
	}

	/// @brief Has a particular custom Rama energy table been loaded?
	/// @details Custom Rama tables are lazily loaded.
	inline bool has_custom_rama_energy_table( Rama_Table_Type const type ) const {
		if ( type == flat_d_aa_ramatable ) return (extra_ram_energ_.count(flat_l_aa_ramatable)!=0); //Special case: D-amino acids use the L-table, inverted, to save RAM.
		return (extra_ram_energ_.count(type)!=0);
	}

	/// @brief Has a particular custom Rama cumulative distribution function been generated?
	/// @details Custom Rama tables are lazily loaded, and cumulative distribution functions are lazily generated.
	inline bool has_custom_rama_cdf( Rama_Table_Type const type ) const {
		if ( type == flat_d_aa_ramatable ) return (extra_cdf_.count(flat_l_aa_ramatable)!=0); //Special case: D-amino acids use the L-table, inverted, to save RAM.
		return (extra_cdf_.count(type)!=0);
	}

	/// @brief Generate a custom Rama cumulative distribution function from the corresponding energy table.
	/// @details This is generated from the ENERGY table, not from the COUNTS.  If the energy table has not been
	/// loaded, this function loads the energy table first.
	void generate_custom_rama_cdf( Rama_Table_Type const type );

	/// @brief Get a custom (extra) cumulative distribution function (CDF).
	/// @details This function is const.  The CDF must already have been generated.
	inline
	utility::vector1 < core::Real > const & custom_rama_cdf( Rama_Table_Type const type ) const {
		runtime_assert_string_msg( has_custom_rama_cdf(type), "Error in core::scoring::Ramachandran::custom_rama_cdf(): The custom CDF requested has not yet been generated." );
		return extra_cdf_.at(type);
	};

	/// @brief Get a custom (extra) rama probability table.
	/// @details This function is const.  The probability table must already have been loaded.  To keep this consistent with
	/// the older Rama probability tables, the phi- and psi-indices are zero-based (running from 0 to 35).
	inline
	core::Real const & custom_rama_probability_table( Rama_Table_Type const type, core::Size const phi, core::Size const psi ) const {
		runtime_assert_string_msg( has_custom_rama_probability_table( type ), "Error in core::scoring::Ramachandran::custom_rama_probability_table(): The custom probability table requested has not yet been loaded." );
		runtime_assert_string_msg( phi < static_cast<core::Size>(n_phi_), "Error in core::scoring::Ramachandran::custom_rama_probability_table(): phi index is out of range." );
		runtime_assert_string_msg( psi < static_cast<core::Size>(n_psi_), "Error in core::scoring::Ramachandran::custom_rama_probability_table(): psi index is out of range." );
		return extra_ram_probabil_.at( type )[phi+1][psi+1]; //The plus 1 is because the array is 1-based, but the interface must be zero-based for backwards consistency.  Ugh.
	};

private: // data

	/// @brief The Ramachandran probability tables.
	/// @details  This is a FORTRAN-style 4D array.  The first two dimensions are for
	/// phi and psi, respectively.  The third is for secondary structure (helix, sheet, loop), though
	/// only the loop values are used nowadays.  The fourth is for the amino acid type (20 entries).
	/// This is a horrible data structure that should be replaced by something much more intuitive.
	ObjexxFCL::FArray4D< Real > ram_probabil_;

	/// @brief Extra Ramachandran count tables.
	/// @details This is a map of (energy table type -> 2D array of probabilities by phi, psi).
	std::map< Rama_Table_Type, utility::vector1< utility::vector1 < core::Real > > > extra_ram_probabil_;

	/// @brief The Ramachandran count tables.
	/// @details  This is a FORTRAN-style 4D array.  The first two dimensions are for
	/// phi and psi, respectively.  The third is for secondary structure (helix, sheet, loop), though
	/// only the loop values are used nowadays.  The fourth is for the amino acid type (20 entries).
	/// This is another horrible data structure that should be replaced by something much more intuitive.
	ObjexxFCL::FArray4D_int ram_counts_;

	/// @brief Extra Ramachandran count tables.
	/// @details This is a map of (energy table type -> 2D array of counts by phi, psi).
	std::map< Rama_Table_Type, utility::vector1< utility::vector1 < core::Size > > > extra_ram_counts_;

	/// @brief The Ramachandran energy tables.
	/// @details  This is a FORTRAN-style 4D array.  The first two dimensions are for
	/// phi and psi, respectively.  The third is for secondary structure (helix, sheet, loop), though
	/// only the loop values are used nowadays.  The fourth is for the amino acid type (20 entries).
	/// This is also a horrible data structure that should be replaced by something much more intuitive.
	ObjexxFCL::FArray4D< Real > ram_energ_;

	/// @brief Extra Ramachandran energy tables.
	/// @details This is a map of (energy table type -> 2D array of energies by phi, psi).
	std::map< Rama_Table_Type, utility::vector1< utility::vector1 < core::Real > > > extra_ram_energ_;

	ObjexxFCL::FArray2D< Real > ram_entropy_;

	// loops only -- does not contain bicubic splines for the alpha or beta secondary structures
	// in fact, there is basically nothing in rosetta that uses the alpha or beta ramachandran tables
	utility::vector1< numeric::interpolation::spline::BicubicSpline > rama_energy_splines_;

	static int const n_phi_ = 36;
	static int const n_psi_ = 36;
	static Real const binw_; // 360 / n_phi_ = 10;
	static Real const rama_sampling_thold_;
	static Real const rama_sampling_factor_;
	static int const n_aa_ = 20; // Ramachandran score defined for the cananical AAs only.

	/// @brief The cumulative distribution functions (CDF) for all phi/psi bins for each amino acid.
	///
	utility::vector1< utility::vector1< Real > > cdf_;

	/// @brief Additional cumulative distribution functions (CDFs).
	/// @details Used for several things.  For example, for storing CDFs for flat L-alpha amino acid
	/// sampling.
	std::map< Rama_Table_Type, utility::vector1< Real > > extra_cdf_;

	// The CDF for all phi/psi bins given the left AA, the center AA, and the secondary-structure classification
	// (ABEGO or ABEGX?).  "Torsion bin" here means the ABEGO classification.  This classification is encoded
	// in core/conformation/util.cc  A phi/psi bin that doesn't land in the particular ss classification is
	// assigned the probability of the previous bin so that it will not be selected if you sample phi/psi bins
	// from a uniform probability distribution.
	ObjexxFCL::FArray2D< utility::vector1< Real > > cdf_by_torsion_bin_;

	// The count of the number of valid phi/psi bins for a particular torsion bin for each amino acid.
	ObjexxFCL::FArray2D< Size > n_valid_pp_bins_by_ppo_torbin_;

	// The phi/psi bin indexes with probability greater than the rama_sampling_thold_ value for each amino acid.
	utility::vector1< utility::vector1< Size > > phi_psi_bins_above_thold_;

};

}
}

#endif
