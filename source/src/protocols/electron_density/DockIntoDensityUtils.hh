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
	friend std::ostream& operator<< (std::ostream& out, const RBfitResult& result);
};


// stores the result of refinement
class RefinementResult {
public:
	core::Real score_, prerefine_score_, spharm_score_;
	core::pose::PoseCOP pose_;

	numeric::xyzVector< core::Real > center();

	RefinementResult() {}

	RefinementResult(
		core::Real const score,
		core::pose::PoseOP const pose_in
	) {
		score_ = score;
		prerefine_score_ = 0.0;
		spharm_score_ = 0.0;
		pose_ = pose_in;
	}

	RefinementResult(
		core::Real const score,
		core::Real const prerefine_score,
		core::Real const spharm_score,
		core::pose::PoseCOP const pose_in
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
	bool operator()(RBfitResult const& t1, RBfitResult const& t2) const {
		return (t1.score_ > t2.score_);
	}
};


class RefinementResultComparitor {
public:
	bool operator()(RefinementResult const& t1, RefinementResult const& t2) const {
		return (t1.score_ > t2.score_);
	}
};


class RevRefinementResultComparitor {
public:
	bool operator()(RefinementResult const & t1, RefinementResult const & t2) const {
		return (t1.score_ < t2.score_);
	}
};

class PointScoreComparator{
public:
	bool operator()(
		std::pair< numeric::xyzVector<core::Real>, core::Real > const & Pair1,
		std::pair< numeric::xyzVector<core::Real>, core::Real > const & Pair2) const {
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


typedef ResultDB<RBfitResult, RBfitResultComparitor> RBfitResultDB;
typedef ResultDB<RefinementResult, RefinementResultComparitor> RefinementResultDB;
typedef ResultDB<RefinementResult, RevRefinementResultComparitor> RevRefinementResultDB;


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
	utility::vector1< core::pose::PoseCOP > natives_;
	utility::vector1< numeric::xyzVector< core::Real > > native_coms_;
	utility::vector1< numeric::xyzVector< core::Real > > native_middle_cas_;
};


struct DensityGridSearchOptions {
	core::Size B_, nRsteps_, max_rot_per_trans_, point_search_start_, point_search_end_;
	core::Real delRSteps_, cluster_radius_, laplacian_offset_, rms_cutoff_;
	bool center_on_middle_ca_, convolute_single_residue_, include_distance_during_fast_cluster_;
	std::string output_fn_;
	DensitySymmInfo symminfo_;
	utility::vector1< core::pose::PoseCOP > natives_;
	utility::vector1< numeric::xyzVector< core::Real > > native_coms_;
	utility::vector1< numeric::xyzVector< core::Real > > native_middle_cas_;
};

/// FUNCTIONS
template< typename T >
void
write_RBfitResultDB( RBfitResultDB fit_result_DB, T & outresults );

template< typename T >
void
dump_RefinementDB_to_silent(
	T resultDB,
	std::string const & outfile,
	std::string const & tag_prefix,
	std::string const & final_chain,
	bool const centroid_output,
	bool const append_to_outfile,
	utility::vector1< core::pose::PoseCOP > const & natives,
	DensitySymmInfo const & symminfo,
	bool const legacy_rms
);

// template
// void
// dump_RefinementDB_to_silent<RevRefinementResultDB>(
//   RevRefinementResultDB resultDB,
//   std::string const & outfile,
//   std::string const & tag_prefix,
//   std::string const & final_chain,
//   bool const centroid_output,
//   bool const append_to_outfile,
//   utility::vector1< core::pose::PoseCOP > const & natives,
//   DensitySymmInfo const & symminfo
// );
//
// template
// void
// dump_RefinementDB_to_silent<RefinementResultDB>(
//   RefinementResultDB resultDB,
//   std::string const & outfile,
//   std::string const & tag_prefix,
//   std::string const & final_chain,
//   bool const centroid_output,
//   bool const append_to_outfile,
//   utility::vector1< core::pose::PoseCOP > const & natives,
//   DensitySymmInfo const & symminfo
// );


void
compare_RBfitDB_to_native(
	RBfitResultDB resultDB,
	core::pose::Pose const & pose,
	core::pose::PoseCOPs const & natives,
	utility::vector1< numeric::xyzVector< core::Real > > const & native_coms,
	utility::vector1< numeric::xyzVector< core::Real > > const & native_middle_cas,
	DensitySymmInfo const & symminfo,
	bool const rot_middle_ca,
	core::Real const rms_cutoff );

// non-superposed RMS
core::Real
get_rms(core::pose::Pose const & r1, core::pose::Pose const & r2, DensitySymmInfo const & d);

core::Real
get_rms(RefinementResult const & r1, RefinementResult const & r2, DensitySymmInfo const & d );

// These respect the poses PDBInfo when doing rms
core::Real
get_rms(core::pose::Pose const & r1, core::pose::Pose const & r2, DensitySymmInfo const &d, bool const native);

core::Real
get_gdt(core::pose::Pose const & r1, core::pose::Pose const & r2, DensitySymmInfo const & d, bool const native);


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
dump_points_to_search_to_pdb_or_txt(
	utility::vector1< numeric::xyzVector< core::Real > > const & points_to_search,
	std::string const & pdb_filename,
	std::string const & txt_filename );

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

/// @brief cluster a RefinementResultDB based on a specific cluster_radius.  If target_size is unset or
//         set to zero then the database will not be filtered to a specific size.
void
cluster_RefinementDB( RefinementResultDB & results, DensitySymmInfo const & symminfo, core::Real const cluster_radius, core::Size const target_size = 0);

void
do_filter(
	utility::vector1< core::pose::PoseOP > const &poses,
	RBfitResultDB & results,
	bool const rescore,
	DensitySymmInfo const & symminfo,
	core::Real const cluster_radius
);


void
cluster_RBfitResultDB_fast(
	RBfitResultDB & results,
	core::Size const delR,
	core::Size const nRsteps,
	core::Real const cluster_radius,
	core::Size const max_results,
	bool const include_distance,
	core::scoring::electron_density::ElectronDensity const & dens);


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
