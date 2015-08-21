// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/SmoothEnvPairPotential.cc
/// @brief  Smooth, differentiable version of centroid env and pair terms
/// @author Frank DiMaio

#ifndef INCLUDED_core_scoring_SmoothEnvPairPotential_hh
#define INCLUDED_core_scoring_SmoothEnvPairPotential_hh

#include <core/scoring/SmoothEnvPairPotential.fwd.hh>

#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

#include <numeric/constants.hh>

#include <basic/datacache/CacheableData.hh>


#include <utility/vector1_bool.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>


namespace core {
namespace scoring {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Keep track of the cenlist information
/// stores both centroid counts (T = Real)
/// as well as d_centroid_counts (T = Vector)
template<class T>
class SigmoidWeightedCenList : public basic::datacache::CacheableData {
public:
	SigmoidWeightedCenList(): calculated_(false) {};
	SigmoidWeightedCenList( SigmoidWeightedCenList const & src ) :
		CacheableData() {
		fcen6_ = src.fcen6_;
		fcen10_ = src.fcen10_;
		fcen12_ = src.fcen12_;
		calculated_ = src.calculated_;
	}

	basic::datacache::CacheableDataOP clone() const {
		return basic::datacache::CacheableDataOP( new SigmoidWeightedCenList( *this ) );
	}

	Size size() const { return fcen6_.size();}

	T fcen6( Size const seqpos ) const { return fcen6_[ seqpos ];  }
	T & fcen6( Size const seqpos ) { return fcen6_[ seqpos ]; }

	T fcen10( Size const seqpos ) const { return fcen10_[ seqpos ]; }
	T & fcen10( Size const seqpos ) { return fcen10_[ seqpos ]; }

	T fcen12( Size const seqpos ) const { return fcen12_[ seqpos ]; }
	T & fcen12( Size const seqpos ) { return fcen12_[ seqpos ]; }

	bool calculated() const { return calculated_; }
	bool & calculated() { return calculated_; }

	void initialize( Size nres, T val ) {
		fcen6_.resize( nres );
		fcen10_.resize( nres );
		fcen12_.resize( nres );
		std::fill( fcen6_.begin(), fcen6_.end(), val );
		std::fill( fcen12_.begin(), fcen12_.end(), val );
		std::fill( fcen10_.begin(), fcen10_.end(), val );
	}

	// Setter functions
	void set_fcen6( Size const seqpos, T value ) {
		fcen6_[ seqpos ] = value;
	}

	void set_fcen10( Size const seqpos, T value ) {
		fcen10_[ seqpos ] = value;
	}

	void set_fcen12( Size const seqpos, T value ) {
		fcen12_[ seqpos ] = value;
	}

private:
	utility::vector1< T > fcen6_;
	utility::vector1< T > fcen10_;
	utility::vector1< T > fcen12_;
	bool calculated_;
};


///////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////
// fpd helper class holds some # of gaussian + sigmoid coeffs + a shift
class SmoothScoreTermCoeffs {
public:
	SmoothScoreTermCoeffs() : shift_(0.0) {}

	void add_sigmoid( numeric::xyzVector< Real > s_in ) {
		sigmoid_coeffs_.push_back( s_in );
	}
	void add_gaussian( numeric::xyzVector< Real > g_in ) {
		// normalize gaussian
		g_in[0] /= sqrt(2*numeric::constants::f::pi*g_in[2]*g_in[2]);
		gaussian_coeffs_.push_back( g_in );
	}
	void shift(Real s_in) { shift_=s_in; }

	Size nsigmoids() { return sigmoid_coeffs_.size(); }
	Size ngaussians() { return gaussian_coeffs_.size(); }

	numeric::xyzVector< Real > sigmoid( Size i ) {
		return sigmoid_coeffs_[i];
	}
	numeric::xyzVector< Real > gaussian( Size i ) {
		return gaussian_coeffs_[i];
	}
	Real shift() { return shift_; }

	void clear() {
		sigmoid_coeffs_.clear();
		gaussian_coeffs_.clear();
		shift_ = 0;
	}

	Real func( Real x ) const;
	Real dfunc( Real x ) const;

private:
	utility::vector1< numeric::xyzVector< Real > > sigmoid_coeffs_;
	utility::vector1< numeric::xyzVector< Real > > gaussian_coeffs_;
	Real shift_;
};

////////////////////////
////////////////////////

class SmoothEnvPairPotential : public utility::pointer::ReferenceCount {
public:
	SmoothEnvPairPotential();


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


	void
	evaluate_env_and_cbeta_scores(
		pose::Pose const & pose,
		conformation::Residue const & rsd,
		Real & env_score,
		Real & cb_score6,
		Real & cb_score12
	) const;


	void
	evaluate_pair_and_cenpack_score(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		Real const cendist,
		Real & pair_contribution,
		Real & cenpack_contribution
	) const;


	void
	evaluate_env_and_cbeta_deriv(
		pose::Pose const & pose,
		conformation::Residue const & rsd,
		numeric::xyzVector<Real> & d_env_score,
		numeric::xyzVector<Real> & d_cb_score6,
		numeric::xyzVector<Real> & d_cb_score12
	) const;


	void
	evaluate_pair_and_cenpack_deriv(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		Real const cendist,
		Real & d_pair,
		Real & d_cenpack
	) const;

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
		numeric::xyzVector<Real> const cendist
	) const;

	//fpd no need for this function anymore; functions are all defined (and well behaved) over all x
	//void truncate_cenlist_values( SigmoidWeightedCenList & cenlist ) const;

private: // data
	Real SIGMOID_SLOPE;

	SmoothScoreTermCoeffs cbeta6_;
	SmoothScoreTermCoeffs cbeta12_;
	SmoothScoreTermCoeffs cenpack_;
	utility::vector1< SmoothScoreTermCoeffs > env_;   // for each restype
	utility::vector1< utility::vector1< SmoothScoreTermCoeffs > > pair_;  // for each restype pair
};

} // ns scoring
} // ns core

#endif
