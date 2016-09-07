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
#include <core/conformation/Residue.fwd.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <numeric/interpolation/spline/Bicubic_spline.hh>

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

	void
	eval_rpp_rama_score(
		AA const res_aa1,
		AA const res_aa2,
		Real const phi,
		Real const psi,
		Real & score_rama,
		Real & denergy_dphi,
		Real & denergy_dpsi
	) const;


protected:

	void read_rpp_tables();

	/// @brief adapted from Max's code, ramachandran.cc
	/// @details If symmetrize_gly is true, the plot for glycine is made symmetric.
	void read_rama_map_file_shapovalov ( std::string const & filename, utility::vector1<  ObjexxFCL::FArray2D< Real > > &data, bool const symmetrize_gly );

	void setup_interpolation( ObjexxFCL::FArray2D< Real > &, numeric::interpolation::spline::BicubicSpline  &, core::Real const &offset);

	/// @brief If the -symmetric_gly_tables option is used, symmetrize the aa_gly table.
	/// @details By default, the gly table is asymmetric because it is based on statistics from the PDB (which disproportionately put glycine
	/// in the D-amino acid region of Ramachandran space).  However, the intrinsic propensities of glycine make it equally inclined to favour
	/// right- or left-handed conformation.  (Glycine is achrial, and can't have a preference.)  Must be called AFTER gly table load, but prior
	/// to bicubic interpolation setup.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void symmetrize_gly_table( ObjexxFCL::FArray2D< core::Real > & data, core::Real &entropy ) const;

	// spline interpolated p(phi,psi)
	utility::vector1< numeric::interpolation::spline::BicubicSpline > rama_splines_, rama_pp_splines_;

};

}
}

#endif
