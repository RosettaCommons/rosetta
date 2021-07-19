// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief protocols for folding into density
/// @details
/// @author Frank DiMaio


#ifndef INCLUDED_protocols_electron_density_DockIntoDensityMover_hh
#define INCLUDED_protocols_electron_density_DockIntoDensityMover_hh


#include <core/conformation/Residue.hh> // only necessary here because of functions that ought to be in the .cc
#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>

#include <ObjexxFCL/FArray3D.fwd.hh>

#include <protocols/electron_density/DensitySymmInfo.hh>
#include <protocols/electron_density/DockIntoDensityUtils.hh>
#include <protocols/moves/Mover.hh>

#include <utility/vector1.hh>

//// C++ headers
#include <string>
#include <queue>


namespace protocols {
namespace electron_density {


///
///  The workhorse
class DockIntoDensityMover : public protocols::moves::Mover {
public:
	DockIntoDensityMover() :
		topNtrans_(5000), topNfilter_(1000), topNfinal_(50), delR_(2),
		dens_wt_(20.0), cluster_radius_(2.0),point_radius_(3),fragDens_(0.7), mindist_(3), laplacian_offset_(0), B_(16), nRsteps_(0), gridStep_(2),
		center_on_middle_ca_(false), points_defined_(false), convolute_single_residue_(false), cluster_oversample_(2), max_rot_per_trans_(3),
		do_refine_(true), min_backbone_(true), ncyc_(1), passthrough_(false), native_com_(0,0,0) {}

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
	void setMinDist( core::Real mindist ) { mindist_=mindist; }
	void setMinBackbone( bool min_backbone ) { min_backbone_=min_backbone; }
	void setNCyc( core::Size ncyc ) { ncyc_=ncyc; }
	void setOutputSilent( std::string silent_out ) { silent_ = silent_out; }
	void setClusterRadius( core::Real cluster_radius ) { cluster_radius_ = cluster_radius; }
	void setPointRadius( core::Real point_radius ) { point_radius_ = point_radius; }
	void setTag( std::string tag ) { tag_=tag; } // output tag
	void setFragDens( core::Real fragDens ) { fragDens_=fragDens; } // output tag
	void setPassThrough( bool passthrough ) { passthrough_ = passthrough; }
	void setClusterOversamp( core::Size cluster_oversample ) { cluster_oversample_=cluster_oversample; }
	void setMaxRotPerTrans( core::Size max_rot_per_trans ) { max_rot_per_trans_=max_rot_per_trans; }
	void setSymminfo( DensitySymmInfo const & symminfo ) { symminfo_=symminfo; }
	void setConvoluteSingleR( bool convolute_single_residue ) { convolute_single_residue_ = convolute_single_residue; }
	void setLaplacianOffset( core::Real laplacian_offset ) { laplacian_offset_ = laplacian_offset; }


	// predefine search locations
	void predefine_search( utility::vector1< numeric::xyzVector<core::Real> > &pts_in );
	void setCenterOnMiddleCA( bool val ) { center_on_middle_ca_=val; }

	// set native pose (for reporting)
	void setNative( core::pose::PoseOP native );

	//
	void apply( core::pose::Pose & pose) override;

	// apply to multiple poses simultaneously
	void apply_multi( utility::vector1< core::pose::PoseOP > & poses);

	// convert 1d power spectrum into 3d map
	void
	map_from_spectrum( utility::vector1< core::Real > const& pose_1dspec, ObjexxFCL::FArray3D< core::Real > &rot );

	/// @brief  step 3: local refinement of each hit
	///   empties the results_in DB
	void do_refinement (
		utility::vector1< core::pose::PoseOP > const &poses,
		RBfitResultDB & results_in,
		RefinementResultDB & results_out
	);

	/// @brief   debugging: add some stats
	void print_best_rms( core::pose::Pose const &pose, RBfitResultDB const &results );


	std::string get_name() const override {
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
	core::Real delR_, dens_wt_, cluster_radius_, point_radius_, fragDens_, mindist_, laplacian_offset_;
	core::Size B_;
	core::Size nRsteps_,gridStep_;
	bool center_on_middle_ca_, points_defined_, convolute_single_residue_;

	// symmetry
	DensitySymmInfo symminfo_;

	// cluster params
	core::Size cluster_oversample_, max_rot_per_trans_;

	// refinement options
	bool do_refine_, min_backbone_;
	core::Size ncyc_;

	// silent file output
	std::string silent_, tag_;

	// points to search (can be input or chosen)
	utility::vector1< numeric::xyzVector<core::Real> > points_to_search_;

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
