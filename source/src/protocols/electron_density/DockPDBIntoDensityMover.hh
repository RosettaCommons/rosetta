// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DockPDBIntoDensityMover.hh
/// @brief protocols for folding into density
/// @details Re-tooled DockFragmentsIntoDensityMover to support parallel execution, and more intermediate debugging steps
/// at the cost of only working with PDBs and not fragments.
/// @author Danny Farrell


#ifndef INCLUDED_protocols_electron_density_DockPDBIntoDensityMover_hh
#define INCLUDED_protocols_electron_density_DockPDBIntoDensityMover_hh

#include <basic/datacache/DataMap.fwd.hh>

#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>

#include <ObjexxFCL/FArray3D.hh>

#include <protocols/electron_density/DensitySymmInfo.hh>
#include <protocols/electron_density/DockIntoDensityUtils.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>

#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

#include <string>
#include <queue>


namespace protocols {
namespace electron_density {

struct MinimizePoseIntoDensityOptions {
	core::scoring::ScoreFunctionOP scorefxn_dens, scorefxn_refine, scorefxn_refine_rb;
	protocols::minimization_packing::MinMoverOP rbmin, bbmin;
	protocols::minimization_packing::PackRotamersMoverOP packer;
};

class DockPDBIntoDensityMover : public protocols::moves::Mover {
public:
	DockPDBIntoDensityMover() { // Not complex types so we use this method, even though it's less efficient.
		topNtrans_ = 5000;
		topNfilter_ = 1000;
		topNfinal_ = 50;
		delR_ = 2;
		dens_wt_ = 20.0;
		cluster_radius_ = 2.0;
		point_radius_ = 0;
		fragDens_ = 0.7;
		laplacian_offset_ = 0;
		B_ = 16;
		nRsteps_ = 0;
		gridStep_ = 2;
		rot_seq_center_ = false;
		rot_middle_ca_ = false;
		points_defined_ = false;
		cluster_oversample_ = 2;
		max_rot_per_trans_ = 3;
		do_refine_ = true;
		min_backbone_ = true;
		constrain_refinement_ = 0;
		cart_ref_ = false;
		convolute_single_residue_ = false;
		gaussian_blur_ = false;
		ncyc_ = 1;
		passthrough_ = false;
		centroid_silent_out_ = true;
		score_natives_ = false;
		cheat_native_com_ = false;
		cheat_native_mca_ = false;
		overwrite_ = false;
		points_to_search_fname_ = "";
		points_to_search_pdb_fname_ = "";
		start_model_ = "";
		dump_inter_ = false;
		dump_inter_silent_ = false;
		final_chain_ = "^";
		// TODO: initilize our natives
	}

	// Options

	// Search options
	void setDelR( core::Real delR ) { delR_=delR; }
	void setNRsteps( core::Real nRsteps ) { nRsteps_=nRsteps; }
	void set_nRsteps_from_pose( core::pose::Pose const & pose );
	void setB( core::Real B ) { B_=B; }
	void setTopN( core::Size const topNtrans, core::Size const topNfilter, core::Size const topNfinal ) {
		topNtrans_=topNtrans;
		topNfilter_=topNfilter;
		topNfinal_=topNfinal;
	}
	core::Size get_top_N_filter() const { return topNfilter_; }
	void setGridStep( core::Size gridStep ) { gridStep_=gridStep; }
	void setNCyc( core::Size ncyc ) { ncyc_=ncyc; }
	void setPointRadius( core::Real point_radius ) { point_radius_ = point_radius; }
	void setClusterOversamp( core::Size cluster_oversample ) { cluster_oversample_=cluster_oversample; }
	void setFragDens( core::Real fragDens ) { fragDens_=fragDens; } // output tag
	void setConvoluteSingleR( bool convolute_single_residue ) { convolute_single_residue_ = convolute_single_residue; } // how to make default true?
	void setGaussianBlur( bool gaussian_blur ) { gaussian_blur_ = gaussian_blur; }
	void setLaplacianOffset( core::Real laplacian_offset ) { laplacian_offset_ = laplacian_offset; }
	void setPointsToSearchFname( std::string const & pts_fname ) { points_to_search_fname_ = pts_fname; }
	void setPointsToSearchPDBFname( std::string const & pts_fname ) { points_to_search_pdb_fname_ = pts_fname; }
	void setRotateSeqCenter( bool rot_seq_center ) { rot_seq_center_ = rot_seq_center; }
	void setRotateMiddleCA( bool rot_middle_ca ) { rot_middle_ca_ = rot_middle_ca; }
	void setRotateOnSeqCenter( bool val ) { rot_seq_center_=val; }
	void setStartModel( std::string const & startmodel ) { start_model_ = startmodel; }

	// Refinement options
	void setConstrainRefinement( core::Real constrain_refinement ) { constrain_refinement_=constrain_refinement; }
	void setCartRef( core::Real cart_ref ) { cart_ref_=cart_ref; }
	void setDoRefine( bool do_refine ) { do_refine_=do_refine; }
	void setMinBackbone( bool min_backbone ) { min_backbone_=min_backbone; }

	// set native pose (for reporting)
	void setMultiNative( utility::vector1< core::pose::PoseOP > const & );
	void setScoreNatives( bool score_natives ) { score_natives_ = score_natives; }

	// General options
	void setCentroidSilentOut( bool centroid_silent_out ) { centroid_silent_out_ = centroid_silent_out; }
	void setClusterRadius( core::Real cluster_radius ) { cluster_radius_ = cluster_radius; }
	void setDumpInter( bool dump_inter ) { dump_inter_ = dump_inter; }
	void setDumpInterSilent( bool dump_inter_silent ) { dump_inter_silent_ = dump_inter_silent; }
	void setMaxRotPerTrans( core::Size max_rot_per_trans ) { max_rot_per_trans_=max_rot_per_trans; }
	void setOutputSilent( std::string silent_out ) { silent_ = silent_out; }
	void setOverwrite( bool overwrite ) { overwrite_ = overwrite; }
	void setPassThrough( bool passthrough ) { passthrough_ = passthrough; }
	void setSymminfo( DensitySymmInfo const & symminfo ) { symminfo_=symminfo; }
	void setTag( std::string tag ) { tag_=tag; } // output tag
	void setFinalChain( std::string const & final_chain ) { final_chain_=final_chain; }

	// Parallel options
	void setCoreIdx( utility::vector1< core::Size > core_idx ) { core_idx_ = core_idx; }
	void setRefineEnd( core::Size refine_end ) { refine_end_ = refine_end; }
	void setRefineStart( core::Size refine_start ) { refine_start_ = refine_start; }
	void setSearchEnd( core::Size point_search_end ) { point_search_end_ = point_search_end; }
	void setSearchStart( core::Size point_search_start ) { point_search_start_ = point_search_start; }
	void setPointSearchResultsFname( std::string const & point_search_results_fname ) { point_search_results_fname_ = point_search_results_fname; }
	void setCombinedSearchResultsFname( std::string const & combined_search_results_fname ) { combined_search_results_fname_ = combined_search_results_fname; }


	// predefine search locations
	void predefine_search( utility::vector1< numeric::xyzVector<core::Real> > &pts_in );

	// Cheat options
	void setCheatNativeCOM( bool cheat_native_com ) { cheat_native_com_ = cheat_native_com; }
	void setCheatNativeMCA( bool cheat_native_mca ) { cheat_native_mca_ = cheat_native_mca; }

	// get max extent of pose
	core::Real
	get_radius( core::pose::Pose const & );

	std::string get_name() const override {
		return "DockPDBIntoDensity";
	}

	void
	score_and_dump_natives();

	/// @brief Read in points to search
	void
	read_in_points_to_search();

	core::Real
	compare_and_align_poses( core::pose::Pose &, core::pose::Pose const & );

	RBfitResultDB
	read_in_partial_search_results( utility::vector1< std::string > const & ) const;

	void
	refine_partial_RBfitResultDB( utility::vector1< core::pose::PoseOP > const &,
		RBfitResultDB &,
		RevRefinementResultDB & );

	void
	set_search_responsibilities();

	void
	dump_RBfitDB_to_pdbs( RBfitResultDB & resultDB, core::pose::PoseOP const & ) const;

	// refinement functions
	void
	set_refinement_responsibilities( RBfitResultDB & );

	void
	check_for_existing_output_file( std::string const & ) const;

	void
	import_refinement_silent_files( utility::vector1< std::string > const &, RevRefinementResultDB & );

	void
	cluster_RevRefinementDB( RevRefinementResultDB &, core::Size ); // maybe rename?

	void
	minimize_poseOP_into_density(
		core::pose::PoseOP const,
		MinimizePoseIntoDensityOptions const & params,
		utility::vector1< core::Real > &);

	/// @brief  step 3: local refinement of each hit
	///   empties the results_in DB
	void refine_RBfitResultDB (
		core::pose::PoseOP const &,
		RBfitResultDB &,
		RevRefinementResultDB &
	);

	// Main callable functions
	utility::vector1< numeric::xyzVector<core::Real> >
	get_points_to_search( core::pose::PoseOP const );

	RBfitResultDB
	apply_search( core::pose::Pose & pose, core::Size const result_size );

	void
	combine_search( utility::vector1< std::string > const &, core::pose::Pose const &, RBfitResultDB & );

	void
	search_results_to_pdb( utility::vector1< std::string > const &, core::pose::PoseOP const & );

	void
	manual_refine_pdb( core::pose::PoseOP const & );

	RevRefinementResultDB
	apply_refinement( utility::vector1< std::string > const &, core::pose::PoseOP const, RBfitResultDB &, bool const );

	void
	combine_refinement( utility::vector1< std::string > const &, RevRefinementResultDB & refinement_results );

	void
	cluster_silent( utility::vector1< std::string > const & );

	void
	analyze_RefinementDB(RevRefinementResultDB);

	void
	run_aio( core::pose::PoseOP const);

	void
	apply( core::pose::Pose & ) override;

private:
	/// 3 stages, each with different filter numbers:
	///   1: select points
	///   2: select best by spharm rotation
	///   3: select best after refinement
	core::Size topNtrans_;
	core::Size topNfilter_;
	core::Size topNfinal_;
	// TODO: Document these params

	// params of search
	bool rot_seq_center_, rot_middle_ca_, points_defined_, convolute_single_residue_, gaussian_blur_;
	core::Real delR_, dens_wt_, cluster_radius_, point_radius_, fragDens_, laplacian_offset_;
	core::Size nRsteps_, gridStep_, B_;
	std::string points_to_search_fname_, points_to_search_pdb_fname_, start_model_;
	utility::vector1< numeric::xyzVector<core::Real> > points_to_search_;


	// symmetry
	DensitySymmInfo symminfo_;

	// cluster params
	core::Size cluster_oversample_, max_rot_per_trans_;

	// refinement options
	bool do_refine_, min_backbone_, centroid_silent_out_, dump_inter_, dump_inter_silent_, cart_ref_;
	core::Size ncyc_;
	core::Real constrain_refinement_;

	// silent file output
	std::string silent_, tag_, final_chain_, point_search_results_fname_, combined_search_results_fname_;

	// don't actually search, pass through input as if it resulted from a search
	bool passthrough_;

	// native (for reporting) (maybe have a native strut/class?)
	utility::vector1< core::pose::PoseOP > all_natives_;
	utility::vector1< numeric::xyzVector< core::Real > > all_native_com_;
	utility::vector1< numeric::xyzVector< core::Real > > all_native_mca_; // xyz coords of all CAs closest to COMs
	utility::vector1< numeric::xyzVector< core::Real > > all_native_seq_center_; // xyz coords of all CAs at sequence center
	bool cheat_native_com_; // dock only on native_com (cheat)
	bool cheat_native_mca_; // dock only on native_mca (cheat)
	bool score_natives_;

	// Parallel params
	core::Size point_search_start_, point_search_end_, refine_start_, refine_end_;
	utility::vector1< core::Size > core_idx_;

	// Misc.
	bool overwrite_;
};

typedef utility::pointer::shared_ptr< DockPDBIntoDensityMover > DockPDBIntoDensityMoverOP;

}
}

#endif
