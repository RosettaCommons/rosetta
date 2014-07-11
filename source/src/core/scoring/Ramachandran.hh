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
// AUTO-REMOVED #include <core/scoring/ProteinTorsion.hh>

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
	//	Real & favorable, Real & allowed, Real & generous ) const;

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


private:

	void read_rama(
		std::string const & rama_map_filename,
		bool use_bicubic_interpolation);

	void read_rama_map_file ( utility::io::izstream * iunit );
	void initialize_rama_sampling_tables();

	void init_rama_sampling_table( conformation::ppo_torsion_bin torsion_bin );
	void init_uniform_sampling_table();

	void
	draw_random_phi_psi_from_cdf(
		utility::vector1< Real > const & cdf,
		Real & phi,
		Real & psi
	) const;

private: // data
	ObjexxFCL::FArray4D< Real > ram_probabil_;
	ObjexxFCL::FArray4D_int ram_counts_;
	ObjexxFCL::FArray4D< Real > ram_energ_;
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

	// The cumulative distribution functions (CDF) for all phi/psi bins for each amino acid.
	utility::vector1< utility::vector1< Real > > cdf_;
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
