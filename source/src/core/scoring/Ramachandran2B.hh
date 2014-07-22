// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/Ramachandran2B.hh
/// @brief  Neighbor-dependent Ramachandran potential class delcaration
/// @author Guoli Wang
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_scoring_Ramachandran2B_hh
#define INCLUDED_core_scoring_Ramachandran2B_hh

// Unit Headers
#include <core/scoring/Ramachandran2B.fwd.hh>
//#include <core/scoring/ProteinTorsion.hh>

// Project Headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/ppo_torsion_bin.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray4D.hh>

#include <utility/vector1.hh>

// C++ Headers
#include <map>


namespace core {
namespace scoring {


class Ramachandran2B : public utility::pointer::ReferenceCount
{
public:
	typedef pose::Pose Pose;
	typedef chemical::AA AA;

public:
	Ramachandran2B();
	virtual ~Ramachandran2B() ; // auto-removing definition from header{}

	Real
	eval_rama_score_residue(
		AA const res_aa,
		Real const phi,
		Real const psi
	) const;

	// Guoli Wang
	void
	eval_rama_score_residue(
		conformation::Residue const & res,
		chemical::AA const left_aa,
		chemical::AA const right_aa,
		Real & rama,
		Real & drama_dphi,
		Real & drama_dpsi
	) const;

	Real
	RamaE_Lower(
		conformation::Residue const &rsd,
		chemical::AA const &neighbor,
		Real &drama_dphi,
		Real &drama_dpsi
	) const;

	Real
	RamaE_Lower(
		conformation::Residue const &rsd,
		chemical::AA const &neighbor
	) const;

	Real
	RamaE_Upper(
		conformation::Residue const & rsd,
		chemical::AA const &neighbor,
		Real &drama_dphi,
		Real &drama_dpsi
	) const;

	Real
	RamaE_Upper(
		conformation::Residue const & rsd,
		chemical::AA const &neighbor
	) const;

	Real
	RamaE(
		conformation::Residue const & rsd,
		Real &drama_dphi,
		Real &drama_dpsi
	) const;

	Real
	RamaE(
		conformation::Residue const & rsd
	) const;
	// finished

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
	IdealizeRamaEnergy(
		Real const phi,
		Real const psi,
		Real & rama,
		Real & drama_dphi,
		Real & drama_dpsi,
		Real const entropy,
		ObjexxFCL::
		FArray2A< Real > const & rama_for_res
	) const;

	/// @brief pick a phi/psi pair chosen from the rama2b distribution so that the most
	/// commonly observed phi/psi pairs are chosen at higher probability than the least
	/// commonly observed pairs; use the phi/psi distribution for the central aa (the pos aa)
	/// based on its "left" (n-terminal) neighbor's amino acid type.
	void
	random_phipsi_from_rama_left(
		AA const left_aa,
		AA const pos_aa,
		Real & phi,
		Real & psi
	) const;

	/// @brief pick a phi/psi pair chosen from the rama2b distribution so that the most
	/// commonly observed phi/psi pairs are chosen at higher probability than the least
	/// commonly observed pairs; use the phi/psi distribution for the central aa (the pos aa)
	/// based on its "right" (c-terminal) neighbor's amino acid type.
	void
	random_phipsi_from_rama_right(
		AA const pos_aa,
		AA const right_aa,
		Real & phi,
		Real & psi
	) const;

	/// @brief function for torsion-bin specific but otherwise random phi/psi angles
	/// @author Amelie Stein
	void
	random_phipsi_from_rama_by_torsion_bin_left(
		AA const left_aa,
		AA const pos_aa,
		Real & phi,
		Real & psi,
		conformation::ppo_torsion_bin const torsion_bin
	) const;

	void
	random_phipsi_from_rama_by_torsion_bin_right(
		AA const pos_aa,
		AA const right_aa,
		Real & phi,
		Real & psi,
		conformation::ppo_torsion_bin const torsion_bin
	) const;


	// Undefined, commenting out to fix PyRosetta build  core::Size get_torsion_bin_index(conformation::ppo_torsion_bin torsion_bin) const;

	void
	get_entries_per_torsion_bin_left(
		AA const left_aa,
		AA const pos_aa,
		std::map< conformation::ppo_torsion_bin, core::Size > & tb_frequencies
	) const;

	void
	get_entries_per_torsion_bin_right(
		AA const pos_aa,
		AA const right_aa,
		std::map< conformation::ppo_torsion_bin, core::Size > & tb_frequencies
	) const;

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


	Size n_phi_bins() const;
	Size n_psi_bins() const;

	Real rama_bin_probability_left(  core::chemical::AA aa_left,   core::chemical::AA aa_center, Real phi, Real psi ) const;
	Real rama_bin_probability_right( core::chemical::AA aa_center, core::chemical::AA aa_right,  Real phi, Real psi ) const;
	Real minimum_sampling_probability() const;

	utility::vector1< Real > const & left_cdf(  core::chemical::AA aa_left,   core::chemical::AA aa_center ) const;
	utility::vector1< Real > const & right_cdf( core::chemical::AA aa_center, core::chemical::AA aa_right  ) const;

	utility::vector1< Real > const &
	left_cdf_for_torsion_bin(
		chemical::AA aa_left,
		chemical::AA aa_center,
		conformation::ppo_torsion_bin
	) const;

	utility::vector1< Real > const &
	right_cdf_for_torsion_bin(
		chemical::AA aa_center,
		chemical::AA aa_right,
		conformation::ppo_torsion_bin
	) const;

private:

	void read_rama();

	void initialize_rama_sampling_tables();

	void
	init_rama_sampling_table(
	   const conformation::ppo_torsion_bin torsion_bin,
		 ObjexxFCL::FArray4D< Real > const & ram_probability,
		 ObjexxFCL::FArray2D< utility::vector1< Real > > & cdf,
		 ObjexxFCL::FArray3D< utility::vector1< Real > > & cdf_by_torsion_bin,
		 ObjexxFCL::FArray3D< Size > & n_valid_pp_bins_by_ppo_torbin
	) const;

	void
	draw_random_phi_psi_from_cdf(
		utility::vector1< Real > const & cdf,
		Real & phi,
		Real & psi
	) const;

// Data
private:

	ObjexxFCL::FArray3D< Real > ram_energ_;
	ObjexxFCL::FArray1D< Real > ram_entropy_;
	ObjexxFCL::FArray4D< Real > ram_energ_left_;
	ObjexxFCL::FArray2D< Real > ram_entropy_left_;
	ObjexxFCL::FArray4D< Real > ram_energ_right_;
	ObjexxFCL::FArray2D< Real > ram_entropy_right_;

	static int const n_phi_ = 36;
	static int const n_psi_ = 36;
	static Real const binw_; // 360 / n_phi_ = 10;
	static int const n_aa_ = 20; // Ramachandran score defined for the cananical AAs only.
	static int const nullaa = 21; // Guoli Wang

	Real const rama_score_limit_;

	ObjexxFCL::FArray4D< Real > left_ram_probabil_; // AS: probability of each phi/psi combination given the left & actual  AA
	ObjexxFCL::FArray4D< Real > right_ram_probabil_; // AS: probability of each phi/psi combination given the actual & right  AA

	static Real const rama_sampling_thold_;

	// static Real const rama_sampling_factor_;

	// The cumulative distrubition function (CDF) for all phi/psi bins given the left AA and the center AA.
	ObjexxFCL::FArray2D< utility::vector1< Real > > left_cdf_;

	// The CDF for all phi/psi bins given the left AA, the center AA, and the secondary-structure classification (ABEGO or ABEGX?).
	// "Torsion bin" here means the ABEGO classifiction.  This classification is encoded in core/conformation/util.cc.
	// A phi/psi bin that doesn't land in the particular ss classification is assigned the probability of the previous
	// bin so that it will not be selected if you sample phi/psi bins from a uniform probability distribution
	ObjexxFCL::FArray3D< utility::vector1< Real > > left_cdf_by_torsion_bin_;

	// The count of the number of valid phi/psi bins given the left AA, the center AA, and the secondary structure classification.
	ObjexxFCL::FArray3D< Size > n_valid_left_pp_bins_by_ppo_torbin_;

	// the CDF for all phi/psi bins given the right AA, the center AA
	ObjexxFCL::FArray2D< utility::vector1< Real > > right_cdf_;

	// The CDF for all phi/psi bins given the left AA, the center AA, and the secondayr-structure classification (ABEGO).
	// See comments for left_cdf_by_torsion_bin_;
	ObjexxFCL::FArray3D< utility::vector1< Real > > right_cdf_by_torsion_bin_;

	// The count of the number of valid phi/psi bins given the right AA, the center AA, and the secondary structure classification.
	ObjexxFCL::FArray3D< Size > n_valid_right_pp_bins_by_ppo_torbin_;

};

}
}

#endif
