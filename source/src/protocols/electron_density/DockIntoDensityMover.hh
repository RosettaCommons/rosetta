// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief protocols for folding into density
/// @details
/// @author Frank DiMaio


#ifndef INCLUDED_protocols_electron_density_DockIntoDensityMover_hh
#define INCLUDED_protocols_electron_density_DockIntoDensityMover_hh

#include <basic/datacache/DataMap.fwd.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>

#include <ObjexxFCL/FArray3D.hh>

#include <protocols/loops/Loops.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

//// C++ headers
#include <string>
#include <queue>


namespace protocols {
namespace electron_density {

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
		core::Size pose_idx,
		core::Real score,
		numeric::xyzMatrix<core::Real> rotation,
		numeric::xyzVector<core::Real> pre_trans,
		numeric::xyzVector<core::Real> post_trans
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
	core::Real score_;
	core::pose::PoseOP pose_;

	numeric::xyzVector< core::Real > center() {
		core::Size centerres = (pose_->total_residue()+1)/2;
		return (pose_->residue(centerres).xyz(2));
	}

	RefinementResult() {}

	RefinementResult(
		core::Real score,
		core::pose::PoseOP pose_in
	) {
		score_ = score;
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


///
///  The workhorse
class DockIntoDensityMover : public protocols::moves::Mover {
public:
	DockIntoDensityMover() :
		topNtrans_(5000), topNfilter_(1000), topNfinal_(50), delR_(2),
		dens_wt_(20.0), cluster_radius_(2.0),fragDens_(0.7), B_(16), nRsteps_(0), gridStep_(2),
		center_on_middle_ca_(false), points_defined_(false), cluster_oversample_(2), max_rot_per_trans_(3),
		do_refine_(true), min_backbone_(true), ncyc_(1), normscores_(false), passthrough_(false), native_com_(0,0,0) {}

	// set options
	void setDelR( core::Real delR ) { delR_=delR; }
	void setNRsteps( core::Real nRsteps ) { nRsteps_=nRsteps; }
	void setB( core::Real B ) { B_=B; }
	void setTopN( core::Real topNtrans, core::Real topNfilter, core::Real topNfinal ) {
		topNtrans_=topNtrans;
		topNfilter_=topNfilter;
		topNfinal_=topNfinal;
	}
	void setGridStep( core::Size gridStep ) { gridStep_=gridStep; }
	void setDoRefine( bool do_refine ) { do_refine_=do_refine; }
	void setMinBackbone( bool min_backbone ) { min_backbone_=min_backbone; }
	void setNCyc( core::Size ncyc ) { ncyc_=ncyc; }
	void setOutputSilent( std::string silent_out ) { silent_ = silent_out; }
	void setClusterRadius( core::Real cluster_radius ) { cluster_radius_ = cluster_radius; }
	void setTag( std::string tag ) { tag_=tag; } // output tag
	void setFragDens( core::Real fragDens ) { fragDens_=fragDens; } // output tag
	void setPassThrough( bool passthrough ) { passthrough_ = passthrough; }
	void setNormScores(bool normscores) { normscores_ = normscores; }
	void setClusterOversamp( core::Size cluster_oversample ) { cluster_oversample_=cluster_oversample; }
	void setMaxRotPerTrans( core::Size max_rot_per_trans ) { max_rot_per_trans_=max_rot_per_trans; }

	// predefine search locations
	void predefine_search( utility::vector1< numeric::xyzVector<core::Real> > &pts_in );
	void setCenterOnMiddleCA( bool val ) { center_on_middle_ca_=val; }

	// set native pose (for reporting)
	void setNative( core::pose::PoseOP native );

	//
	void apply( core::pose::Pose & pose);

	// apply to multiple poses simultaneously
	void apply_multi( utility::vector1< core::pose::PoseOP > & poses);

	// get max extent of pose
	core::Real get_radius( core::pose::Pose const & pose, numeric::xyzVector< core::Real > &com );

	// apply an xform to a pose
	void apply_transform (
		core::pose::Pose & pose,
		RBfitResult const& transform
	);

	// get 1d power spectrum of a pose
	void
	get_spectrum( core::pose::Pose const& pose, utility::vector1< core::Real > &pose_1dspec );

	// convert 1d power spectrum into 3d map
	void
	map_from_spectrum( utility::vector1< core::Real > const& pose_1dspec, ObjexxFCL::FArray3D< double > &rot );

	/// @brief  step 0: map pose to spherically sampled density
	void
	poseSphericalSamples(
		core::pose::Pose const &pose,
		ObjexxFCL::FArray3D< double > &sigR);

	/// @brief  step 1: select points over which to search (saved in class variable)
	void
	select_points(
		core::pose::Pose & pose);

	/// @brief  step 2: perform the grid search, storing these results
	void
	density_grid_search (
		core::Size pose_idx,
		core::pose::Pose & pose,
		RBfitResultDB & results);

	/// @brief  step 3: local refinement of each hit
	///   empties the results_in DB
	void do_refinement (
		utility::vector1< core::pose::PoseOP > const &poses,
		RBfitResultDB & results_in,
		RefinementResultDB & results_out
	);

	/// @brief  step 4: filter similar hits (non-refined)
	void
	do_filter(
		utility::vector1< core::pose::PoseOP > const &poses,
		RBfitResultDB & results_in,
		bool rescore
	);

	/// @brief  step 4: (fast) filter hits using rotation only
	void
	do_filter(
		RBfitResultDB & results_in
	);

	/// @brief  step 4: filter similar hits (refined)
	void do_filter(
		RefinementResultDB & results_in
	);

	/// @brief   debugging: add some stats
	void print_best_rms( core::pose::Pose const &pose, RBfitResultDB const &results );


	virtual std::string get_name() const {
		return "DockIntoDensity";
	}

private:
	/// 3 stages, each with different filter numbers:
	///   1: select points
	///   2: select best by spharm rotation
	///   3: select best after refinement
	core::Size topNtrans_;
	core::Size topNfilter_;
	core::Size topNfinal_;

	// params of search
	core::Real delR_, dens_wt_, cluster_radius_, fragDens_;
	core::Size B_;
	core::Size nRsteps_,gridStep_;
	bool center_on_middle_ca_, points_defined_;

	// cluster params
	core::Size cluster_oversample_, max_rot_per_trans_;

	// refinement options
	bool do_refine_, min_backbone_;
	core::Size ncyc_;

	// silent file output
	std::string silent_, tag_;

	// points to search (can be input or chosen)
	utility::vector1< numeric::xyzVector<core::Real> > points_to_search_;
	bool normscores_;

	// don't actually search, pass through input as if it resulted from a search
	bool passthrough_;

	// native (for reporting)
	core::pose::PoseOP native_;
	numeric::xyzVector< core::Real > native_com_;
};

typedef utility::pointer::shared_ptr< DockIntoDensityMover > DockIntoDensityMoverOP;

}
}

#endif
