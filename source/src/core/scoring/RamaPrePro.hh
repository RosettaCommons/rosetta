// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/RamaPrePro.hh
/// @brief  RamaPrePro potential class delcaration
/// @author Frank DiMaio
/// @author Vikram K. Mulligan (vmullig@uw.edu) Feb 2016 -- made this compatible with canonical D-amino acids; returns 0 for noncanonicals.

#ifndef INCLUDED_core_scoring_RamaPrePro_hh
#define INCLUDED_core_scoring_RamaPrePro_hh

// Unit Headers
#include <core/scoring/RamaPrePro.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/izstream.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/mainchain_potential/MainchainScoreTable.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <numeric/interpolation/spline/BicubicSpline.hh>

#include <map>

namespace core {
namespace scoring {

class RamaPrePro : public utility::pointer::ReferenceCount
{
public:
	typedef pose::Pose Pose;
	typedef chemical::AA AA;

public:
	RamaPrePro();
	~RamaPrePro() override = default;

	/// @brief Evaluate the rama score for this residue (res1) given the identity of the next (res_aa2).
	/// @details This version only works for noncanonical or canonical residues with any number of mainchain
	/// torsions.  If the next residue's identity is pro or d-pro, a different score table is used.  Note:
	/// if return_derivs is true, the gradient vector is populated only.  If it is false, then only the
	/// score_rama value is populated.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void
	eval_rpp_rama_score(
		core::chemical::ResidueTypeCOP res1,
		AA const res_aa2,
		utility::vector1 < core::Real > mainchain_torsions, //Deliberately copied, not passed by reference
		Real & score_rama,
		utility::vector1 < core::Real > &gradient,
		bool const return_derivs
	) const;


	/// @brief Evaluate the rama score for this residue (res_aa1) given the identity of the next (res_aa2).
	/// @details This version only works for canonical L-amino acids, canonical D-amino acids, or glycine.  If the next
	/// residue's identity is pro or d-pro, a different score table is used.  Note:
	/// if return_derivs is true, the gradient vector is populated only.  If it is false, then only the
	/// score_rama value is populated.
	/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
	void
	eval_rpp_rama_score(
		AA const res_aa1,
		AA const res_aa2,
		Real const phi,
		Real const psi,
		Real & score_rama,
		Real & denergy_dphi,
		Real & denergy_dpsi,
		bool const return_derivs
	) const;


private: //Private methods.

	/// @brief Ensure that the RamaPrePro scoring tables for the 20 canonical amino acids are set up, and that we are storing
	/// pointers to them in a map of AA enum -> MainchainScoreTableCOP.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void read_canonical_rpp_tables();

	/// @brief Returns true if this aa is aa_pro or aa_dpr, false otherwise.
	///
	bool is_pro( core::chemical::AA const &aa) const;

private: //Private member variables.

	/// @brief Owning pointers to the MainchainScoreTables for the canonical amino acids.
	///
	std::map < core::chemical::AA, core::chemical::mainchain_potential::MainchainScoreTableCOP > canonical_score_tables_;

	/// @brief Owning pointers to the MainchainScoreTables for the canonical amino acids, pre-proline versions.
	///
	std::map < core::chemical::AA, core::chemical::mainchain_potential::MainchainScoreTableCOP > canonical_prepro_score_tables_;

};

}
}

#endif
