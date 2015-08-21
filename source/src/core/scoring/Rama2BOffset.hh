// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/Rama2BOffset.hh
/// @brief  Rama2BOffset potential class delcaration
/// @author

#ifndef INCLUDED_core_scoring_Rama2BOffset_hh
#define INCLUDED_core_scoring_Rama2BOffset_hh

// Unit Headers
#include <core/scoring/Rama2BOffset.fwd.hh>

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

// 4 different classes of AAs for these statistics:
//    G (glycine), P (proline), V (beta-branched), X (others)
// For trans omega we compute stats for all 16 combos
// For cis omega, only XP (preproline) and XX (all others)
enum Rama2BOffsetTables {
	TRANS_XX = 1,
	TRANS_XG,
	TRANS_XP,
	TRANS_XV,
	TRANS_GX,
	TRANS_GG,
	TRANS_GP,
	TRANS_GV,
	TRANS_PX,
	TRANS_PG,
	TRANS_PP,
	TRANS_PV,
	TRANS_VX,
	TRANS_VG,
	TRANS_VP,
	TRANS_VV,
	CIS_XP,
	CIS_XX,
	NRAMATABLES=CIS_XX
};

class Rama2BOffset : public utility::pointer::ReferenceCount
{
public:
	typedef pose::Pose Pose;
	typedef chemical::AA AA;

public:
	Rama2BOffset();
	~Rama2BOffset() {}

	void
	eval_r2bo_rama_score(
		AA const res_aa1,
		AA const res_aa2,
		Real const psi,
		Real const omega,
		Real const phi,
		Real & score_rama,
		Real & denergy_dpsi1,
		Real & denergy_domega2,
		Real & denergy_dphi2
	) const;

	void
	eval_r2bo_omega_score(
		AA const res_aa1,
		AA const res_aa2,
		Real const psi,
		Real const omega,
		Real const phi,
		Real & score_rama,
		Real & denergy_dpsi1,
		Real & denergy_domega2,
		Real & denergy_dphi2
	) const;

	void
	eval_p_aa_pp_score(
		AA const res_aa1,
		AA const res_aa2,
		Real const psi,
		Real const phi,
		Real & score_paapp,
		Real & denergy_dpsi1,
		Real & denergy_dphi2
	) const;

protected:

	Size
	aapair_to_table_index(
		chemical::AA const res_aa1,
		chemical::AA const res_aa2,
		bool cis
	) const;

	void read_r2bo_tables();
	void read_paapp_tables();

	void read_table_from_stream( utility::io::izstream &, ObjexxFCL::FArray2D< Real > &, ObjexxFCL::FArray2D< Real > &, ObjexxFCL::FArray2D< Real > &);
	void read_paapp_table_from_stream( utility::io::izstream &, utility::vector1 <ObjexxFCL::FArray2D< Real > > &);

	void setup_interpolation( ObjexxFCL::FArray2D< Real > &, numeric::interpolation::spline::BicubicSpline  &);

	// spline interpolated p(phi,psi)
	utility::vector1< numeric::interpolation::spline::BicubicSpline > rama_;

	// spline interpolated p(omega|phi,psi)
	utility::vector1< numeric::interpolation::spline::BicubicSpline > omega_mu_,  omega_sig_;

	// spline interpolated p_aa_pp_offset tables
	utility::vector1< numeric::interpolation::spline::BicubicSpline > paapp1_;
	utility::vector1< numeric::interpolation::spline::BicubicSpline > paapp2_;
};

}
}

#endif
