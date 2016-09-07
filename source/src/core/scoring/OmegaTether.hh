// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/OmegaTether.hh
/// @brief  OmegaTether potential class delcaration
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_scoring_OmegaTether_hh
#define INCLUDED_core_scoring_OmegaTether_hh

// Unit Headers
#include <core/scoring/OmegaTether.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/izstream.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.fwd.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <numeric/interpolation/spline/Bicubic_spline.hh>

namespace core {
namespace scoring {


class OmegaTether : public utility::pointer::ReferenceCount
{
public:
	typedef pose::Pose Pose;
	typedef chemical::AA AA;

public:
	OmegaTether();
	~OmegaTether() override = default;

	Real
	eval_omega_score_residue(
		AA const res_aa,
		Real const omega,
		Real const phi,
		Real const psi
	) const;

	void
	eval_omega_score_residue(
		conformation::Residue const & res,
		Real & energy,
		Real & denergy_domega,
		Real & denergy_dphi,
		Real & denergy_dpsi
	) const;

	void
	eval_omega_score_residue(
		AA const res_aa,
		Real const omega,
		Real const phi,
		Real const psi,
		Real & energy,
		Real & denergy_domega,
		Real & denergy_dphi,
		Real & denergy_dpsi
	) const;


	void
	eval_omega_score_all(
		Pose & pose,
		ScoreFunction const & scorefxn
	) const;

	/// @brief Returns the mainchain torsion index corresponding to "phi".
	/// @details Generally 1.  Set to 2 for beta-amino acids so that derivatives are calculated
	/// for the dihedral two spaces before the peptide bond.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	core::Size phi_index( core::conformation::Residue const &rsd ) const;

	/// @brief Returns the mainchain torsion index corresponding to "psi".
	/// @details Generally 2 (alpha-amino acids) or 3 (beta-amino acids).
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	core::Size psi_index( core::conformation::Residue const &rsd ) const;

	/// @brief Returns the mainchain torsion index corresponding to "omega".
	/// @details Should be 3 for alpha amino acids, 4 for beta amino acids.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	core::Size omega_index( core::conformation::Residue const &rsd ) const;


private:

	void read_omega_tables();
	void read_table_from_stream( utility::io::izstream &, ObjexxFCL::FArray2D< Real > &, ObjexxFCL::FArray2D< Real > &);
	void setup_interpolation( ObjexxFCL::FArray2D< Real > &, numeric::interpolation::spline::BicubicSpline  &);

	// phi-psi dependent only
	bool use_phipsi_dep_;
	utility::vector1< ObjexxFCL::FArray2D< core::Real > > omega_mus_all_, omega_sigmas_all_;
	utility::vector1< numeric::interpolation::spline::BicubicSpline > omega_mus_all_splines_,  omega_sigmas_all_splines_;
};

}
}

#endif
