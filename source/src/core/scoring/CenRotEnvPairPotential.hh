// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/CenRotEnvPairPotential.cc
/// @brief  CenRot version of cen pair
/// @author Yuan Liu

#ifndef INCLUDED_core_scoring_CenRotEnvPairPotential_hh
#define INCLUDED_core_scoring_CenRotEnvPairPotential_hh

#include <core/scoring/CenRotEnvPairPotential.fwd.hh>

#include <core/scoring/SmoothEnvPairPotential.hh>
#include <numeric/interpolation/spline/Cubic_spline.hh>
#include <numeric/interpolation/spline/Bicubic_spline.hh>

#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/types.hh>
#include <core/chemical/AA.hh>

#include <numeric/constants.hh>
#include <numeric/types.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>


namespace core {
namespace scoring {

class CenRotEnvPairPotential : public utility::pointer::ReferenceCount {
public:
	CenRotEnvPairPotential();

	///pair
	void
	evaluate_cen_rot_pair_score(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		Real const cendist,
		Real & pair_contribution
	) const;

	void
	evaluate_cen_rot_pair_deriv(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		Real const cendist,
		Real & d_pair
	) const;

	void
	evaluate_cen_rot_pair_orientation_score(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		Real const cendist,
		Real & ang1_contribution,
		Real & ang2_contribution,
		Real & dih_contribution
	) const;

	void
	evaluate_cen_rot_pair_orientation_deriv(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		Real const cendist, Real const ang,
		Real & dE_dr,
		Real & dE_d_ang,
		Real & dE_d_dih
	) const;

	///env
	void
	evaluate_cen_rot_env_and_cbeta_score(
		pose::Pose const & pose,
		conformation::Residue const & rsd,
		Real & env_contribution,
		Real & cbeta6_contribution,
		Real & cbeta12_contribution
	) const ;

	void
	evaluate_cen_rot_env_and_cbeta_deriv(
		pose::Pose const & pose,
		conformation::Residue const & rsd,
		numeric::xyzVector<Real> & f2_cen_env,
		numeric::xyzVector<Real> & f2_cen_cb6,
		numeric::xyzVector<Real> & f2_cen_cb12,
		numeric::xyzVector<Real> & f2_cb_env,
		numeric::xyzVector<Real> & f2_cb_cb6,
		numeric::xyzVector<Real> & f2_cb_cb12
	) const ;

	//copy from smoothenvpairpotential
	void
	compute_centroid_environment(
		pose::Pose & pose
	) const;
	void
	compute_dcentroid_environment(
		pose::Pose & pose
	) const;
	void
	finalize( pose::Pose & pose ) const;

protected:
	Real cen_dist_cutoff_12_pad;

	SigmoidWeightedCenList< Real > const & cenlist_from_pose( pose::Pose const & ) const;
	SigmoidWeightedCenList< Real > & nonconst_cenlist_from_pose( pose::Pose & ) const;

	SigmoidWeightedCenList< numeric::xyzVector<Real> > const & dcenlist_from_pose( pose::Pose const & ) const;
	SigmoidWeightedCenList< numeric::xyzVector<Real> > & nonconst_dcenlist_from_pose( pose::Pose & ) const;

private:
	void
	fill_smooth_cenlist(
		SigmoidWeightedCenList< Real > & cenlist,
		Size const res1,
		Size const res2,
		Real const cendist
	) const;

	void
	fill_smooth_dcenlist(
		SigmoidWeightedCenList< numeric::xyzVector<Real> > & dcenlist,
		Size const res1,
		Size const res2,
		numeric::xyzVector<Real> const cenvec
	) const;

private:
	Real SIGMOID_SLOPE;

	utility::vector1< utility::vector1<
		numeric::interpolation::spline::CubicSpline > > pairsplines_;

	utility::vector1<numeric::interpolation::spline::CubicSpline> envsplines_;
	//numeric::interpolation::spline::CubicSpline cbeta6splines_;
	//numeric::interpolation::spline::CubicSpline cbeta12splines_;
	SmoothScoreTermCoeffs cbeta6_;
	SmoothScoreTermCoeffs cbeta12_;
	SmoothScoreTermCoeffs cenpack_;

	utility::vector1< utility::vector1<
		numeric::interpolation::spline::BicubicSpline > > angsplines_;

	//currently not implemented
	//utility::vector1< utility::vector1<
	// numeric::interpolation::spline::BicubicSpline > > dihsplines_;
};

} // ns scoring
} // ns core

#endif
