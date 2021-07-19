// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file source/src/protocols/electron_density/DockIntoDensityUtils.hh
/// @brief Utils for folding into density
/// @details
/// @author Frank DiMaio
/// @author Danny Farrell


#ifndef INCLUDED_protocols_electron_density_DockIntoDensityUtils_hh
#define INCLUDED_protocols_electron_density_DockIntoDensityUtils_hh

#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <protocols/electron_density/DensitySymmInfo.hh>

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <ObjexxFCL/FArray3D.hh>

#include <queue>


namespace protocols {
namespace electron_density {

/// CLASSES

// stores the intermediate result of spharm search
class RBfitResult {
public:
	core::Size pose_idx_;
	core::Real score_;
	numeric::xyzVector<core::Real> pre_trans_;
	numeric::xyzVector<core::Real> post_trans_;
	numeric::xyzMatrix<core::Real> rotation_;

	RBfitResult() {}

	RBfitResult(
		core::Size const pose_idx,
		core::Real const score,
		numeric::xyzMatrix<core::Real> const & rotation,
		numeric::xyzVector<core::Real> const & pre_trans,
		numeric::xyzVector<core::Real> const & post_trans
	) {
		pose_idx_ = pose_idx;
		score_ = score;
		pre_trans_ = pre_trans;
		post_trans_ = post_trans;
		rotation_ = rotation;
	}
};


// stores the result of refinement
class RefinementResult {
public:
	core::Real score_, prerefine_score_, spharm_score_;
	core::pose::PoseOP pose_;

	numeric::xyzVector< core::Real > center();

	RefinementResult() {}

	RefinementResult(
		core::Real const score,
		core::Real const prerefine_score,
		core::Real const spharm_score,
		core::pose::PoseOP const pose_in
	) {
		score_ = score;
		prerefine_score_ = prerefine_score;
		spharm_score_ = spharm_score;
		pose_ = pose_in;
	}
};


// comparators
class RBfitResultComparitor {
public:
	bool operator()(RBfitResult& t1, RBfitResult& t2) {
		return (t1.score_ > t2.score_);
	}

	bool operator()(RBfitResult const& t1, RBfitResult const& t2) {
		return (t1.score_ > t2.score_);
	}
};


class RefinementResultComparitor {
public:
	bool operator()(RefinementResult& t1, RefinementResult& t2) {
		return (t1.score_ > t2.score_);
	}

	bool operator()(RefinementResult const& t1, RefinementResult const& t2) {
		return (t1.score_ > t2.score_);
	}
};

class PointScoreComparator{
public:
	bool operator()(std::pair< numeric::xyzVector<core::Real>, core::Real > Pair1, std::pair< numeric::xyzVector<core::Real>, core::Real > Pair2){
		return (Pair1.second < Pair2.second);
	}

};

// database stores the top N results
template <class T,class Tcomp> class ResultDB  {
private:
	core::Size N_;
	std::priority_queue<T, std::vector<T>, Tcomp> queue_;

public:
	ResultDB(int N_in) { N_ = N_in; }

	void add_element( T elt ) {
		if ( queue_.size() <N_ ) {
			queue_.push( elt );
		} else if ( elt.score_ > queue_.top().score_ ) {
			queue_.pop();
			queue_.push( elt );
		}
	}

	// would we add an element with this score?
	bool to_add_element( core::Real score ) {
		return (queue_.size() <N_ || score > queue_.top().score_ );
	}

	T pop() {
		T retval = queue_.top();
		queue_.pop();
		return retval;
	}

	T top() {
		T retval = queue_.top();
		return retval;
	}

	core::Size size() { return queue_.size(); }
};


typedef ResultDB<RBfitResult,RBfitResultComparitor> RBfitResultDB;
typedef ResultDB<RefinementResult,RefinementResultComparitor> RefinementResultDB;


struct PoseSphericalSamplesOptions {
	core::Size B_, nRsteps_;
	core::Real delRsteps_, laplacian_offset_;
	bool center_on_middle_ca_;
};


struct SelectDensityPointsOptions {
	core::Size B_, nRsteps_, gridStep_, topNtrans_;
	core::Real delRsteps_, fragDens_, point_radius_, laplacian_offset_;
	bool center_on_middle_ca_, convolute_single_residue_;
	DensitySymmInfo symminfo_;
	core::pose::PoseOP native_;
	numeric::xyzVector< core::Real > native_com_;
};


struct DensityGridSearchOptions {
	core::Size B_, nRsteps_, max_rot_per_trans_;
	core::Real delRSteps_, cluster_radius_, laplacian_offset_;
	bool center_on_middle_ca_, convolute_single_residue_;
	DensitySymmInfo symminfo_;
	core::pose::PoseOP native_;
	numeric::xyzVector< core::Real > native_com_;
};


/// FUNCTIONS


core::Real
get_rms(core::pose::PoseOP const r1, core::pose::PoseOP const r2, DensitySymmInfo const &d);

// non-superposed RMS
core::Real
get_rms(RefinementResult const & r1, RefinementResult const & r2, DensitySymmInfo const & d );

void
apply_transform(
	core::pose::Pose & pose,
	RBfitResult const & transform);

core::Real get_rot_angle( numeric::xyzMatrix<core::Real> const & R );

/// @brief Step 0: map pose to spherically sampled density + mask
void
pose_spherical_samples(
	core::pose::Pose const &pose,
	ObjexxFCL::FArray3D< core::Real > &sigR,
	ObjexxFCL::FArray3D< core::Real > &epsR,
	PoseSphericalSamplesOptions const & params);

/// @brief Get max extent of pose
core::Real
get_radius(
	core::pose::Pose const & pose,
	numeric::xyzVector< core::Real > &com,
	bool const center_on_middle_ca );


void
map_from_spectrum(
	utility::vector1< core::Real > const& pose_1dspec,
	ObjexxFCL::FArray3D< core::Real > &rot,
	core::Real const delR );


void
dump_points_to_search_to_pdb(
	utility::vector1< numeric::xyzVector< core::Real > > const & points_to_search,
	std::string const & filename );

utility::vector1< numeric::xyzVector< core::Real > >
select_density_points( core::pose::Pose const & pose,
	SelectDensityPointsOptions const & params,
	core::Size & nRsteps );

/// @brief Get 1d power spectrum of a pose
core::Real
get_spectrum(
	core::pose::Pose const& pose,
	utility::vector1< core::Real > &pose_1dspec,
	core::Real const delR,
	core::Real const fragDens,
	bool const convolute_single_residue,
	bool const center_on_middle_ca);



void
do_filter( RefinementResultDB & results, DensitySymmInfo const & symminfo, core::Real const cluster_radius);


void
do_filter(
	utility::vector1< core::pose::PoseOP > const &poses,
	RBfitResultDB & results,
	bool const rescore,
	DensitySymmInfo const & symminfo,
	core::Real const cluster_radius
);


void
do_filter(
	RBfitResultDB & results,
	core::Size const delR,
	core::Size const nRsteps,
	core::Real const cluster_radius);


/// @brief do the main search over the map
void
density_grid_search(
	core::Size pose_idx,
	core::pose::Pose const & pose,
	RBfitResultDB & results,
	utility::vector1< numeric::xyzVector< core::Real > > const & points_to_search,
	DensityGridSearchOptions const & params
	//
);


}  // electron_density
}  // protocols
#endif
