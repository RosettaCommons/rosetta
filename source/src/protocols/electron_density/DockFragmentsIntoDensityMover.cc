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

#include <protocols/electron_density/DockFragmentsIntoDensityMover.hh>
#include <protocols/electron_density/DockIntoDensityUtils.hh>

#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/AtomInformation.hh>
#include <core/io/StructFileRep.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/util.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/fourier/SHT.hh>


#include <protocols/electron_density/DockFragmentsIntoDensityMover.hh>
#include <protocols/jd2/util.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <utility/file/file_sys_util.hh>
#include <utility/file/FileName.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


#include <basic/Tracer.hh>

#include <core/chemical/AtomTypeSet.hh> // AUTO IWYU For AtomTypeSet
#include <core/io/StructFileRepOptions.hh> // AUTO IWYU For StructFileRepOptions
#include <core/scoring/electron_density/xray_scattering.hh> // AUTO IWYU For OneGaussianScattering, get_A




namespace protocols {
namespace electron_density {

using core::Size;

static basic::Tracer TR( "protocols.electron_density.DockFragmentsIntoDensityMover" );


void
DockFragmentsIntoDensityMover::print_best_rms( core::pose::Pose const & pose, RBfitResultDB const & results ) {
	core::Real bestrms=1e4;
	core::Size bestrank=0;
	RBfitResultDB resultscopy = results;
	while ( resultscopy.size() > 0 ) {
		RBfitResult const sol_i = resultscopy.pop();
		core::pose::Pose posecopy(pose);
		apply_transform( posecopy, sol_i );
		core::pose::addVirtualResAsRoot( posecopy );
		core::Real const rms_i = get_rms(*native_, posecopy, symminfo_);
		if ( rms_i < bestrms ) {
			bestrms = rms_i; bestrank=resultscopy.size()+1;
		}
	}
	TR << "Best RMS = " << bestrms << " at rank " << bestrank << std::endl;
}


void
DockFragmentsIntoDensityMover::setNative( core::pose::PoseOP native ) {
	core::pose::addVirtualResAsRoot( *native );
	native_ = utility::pointer::make_shared< const core::pose::Pose >(*native);

	native_com_ = numeric::xyzVector< core::Real >(0,0,0);
	core::Size N=0;
	for ( int i=1; i<=(int)native_->size(); ++i ) {
		if ( !native_->residue(i).is_protein() ) continue;
		native_com_ += native_->residue(i).xyz(2);
		N++;
	}
	native_com_ /= N;


	native_middle_ca_ = numeric::xyzVector< core::Real >(0,0,0);
	core::Real closest_distance_to_COM(99999);
	for ( core::Size i=1; i <= native_->size(); ++i ) {
		if ( native_->residue(i).aa() == core::chemical::aa_vrt ) continue;

		core::Real const dist = native_->residue(i).atom(2).xyz().distance( native_com_ );

		if (  dist < closest_distance_to_COM ) {
			native_middle_ca_ = native_->residue(i).atom(2).xyz(); // set mca
			closest_distance_to_COM = dist;
		}
	}
}


void
DockFragmentsIntoDensityMover::map_from_spectrum( utility::vector1< core::Real > const& pose_1dspec, ObjexxFCL::FArray3D< core::Real > &rot ) {
	// ugly
	for ( int z=1; z<=(int)rot.u3(); ++z ) {
		for ( int y=1; y<=(int)rot.u2(); ++y ) {
			for ( int x=1; x<=(int)rot.u1(); ++x ) {
				numeric::xyzVector< core::Real > const idxX(
					x<=rot.u1()/2 ? x-1 : x-1-rot.u1(),
					y<=rot.u2()/2 ? y-1 : y-1-rot.u2(),
					z<=rot.u3()/2 ? z-1 : z-1-rot.u3()
				);

				numeric::xyzVector< core::Real > const cartX = [&]{
					numeric::xyzVector< core::Real > out;
					core::scoring::electron_density::getDensityMap().idxoffset2cart( idxX, out );
					return out;
				}();

				core::Real const d = cartX.length() / delR_;
				core::Real const fpart = d - std::floor(d);
				core::Size const dint = (core::Size) std::floor(d) + 1;
				if ( dint<pose_1dspec.size() ) { // last entry is always 0 so this check is valid
					rot(x,y,z) = (1-fpart)*pose_1dspec[dint] + (fpart)*pose_1dspec[dint+1];
				} else {
					rot(x,y,z) = 0.0;
				}
			}
		}
	}

	if ( basic::options::option[ basic::options::OptionKeys::edensity::debug ]() ) {
		core::scoring::electron_density::ElectronDensity(rot,1.0, numeric::xyzVector< core::Real >(0,0,0), false).writeMRC( "spectrum.mrc" );
	}
}


void
DockFragmentsIntoDensityMover::predefine_search( utility::vector1< numeric::xyzVector<core::Real> > &pts_in ) {
	points_defined_ = true;

	points_to_search_.clear();

	for ( int i=1; i<=(int)pts_in.size(); ++i ) {
		numeric::xyzVector< core::Real > x_idx;
		core::scoring::electron_density::getDensityMap().cart2idx( pts_in[i] , x_idx );
		points_to_search_.push_back( x_idx );
	}
}


void
DockFragmentsIntoDensityMover::do_refinement (
	utility::vector1< core::pose::PoseOP > const &poses,
	RBfitResultDB & results_in,
	RefinementResultDB & results_out
) {
	using namespace core::pack;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;

	core::scoring::ScoreFunctionOP scorefxn_dens( new core::scoring::ScoreFunction() );
	scorefxn_dens->set_weight( core::scoring::elec_dens_whole_structure_allatom, 1.0);

	core::scoring::ScoreFunctionOP scorefxn_refine = core::scoring::get_score_function();
	scorefxn_refine->set_weight( core::scoring::elec_dens_fast, dens_wt_);
	core::scoring::ScoreFunctionOP scorefxn_refine_rb( new core::scoring::ScoreFunction() );
	scorefxn_refine_rb->set_weight( core::scoring::elec_dens_fast, dens_wt_);

	core::kinematics::MoveMapOP bbmm( new core::kinematics::MoveMap );
	bbmm->set_bb( true ); bbmm->set_chi( true ); bbmm->set_jump( true );
	protocols::minimization_packing::MinMoverOP bbmin( new protocols::minimization_packing::MinMover(
		bbmm, scorefxn_refine, "lbfgs_armijo_nonmonotone", 1e-5, true ) );
	bbmin->max_iter(200); // make a parameter?

	// packer
	TaskFactoryOP tf( new TaskFactory() );
	tf->push_back( utility::pointer::make_shared< InitializeFromCommandline >()); // get extra rotamer flags from command line
	tf->push_back( utility::pointer::make_shared< operation::IncludeCurrent >()); // include current rotamer by default
	tf->push_back( utility::pointer::make_shared< RestrictToRepacking >()); // do not design
	protocols::minimization_packing::PackRotamersMoverOP packer( new protocols::minimization_packing::PackRotamersMover() );
	packer->task_factory( tf );
	packer->score_function( scorefxn_refine );

	protocols::simple_moves::SwitchResidueTypeSetMover to_all_atom( core::chemical::FA_STANDARD );

	// do refinement
	core::Size ntotal=results_in.size();
	while ( results_in.size() > 0 ) {
		RBfitResult sol_i = results_in.pop();
		//TR << "[" << ntotal-results_in.size() << "/" << ntotal << "]" << "\r" << std::flush;

		core::pose::PoseOP posecopy ( new core::pose::Pose( *(poses[ sol_i.pose_idx_ ]) ) );
		apply_transform( *posecopy, sol_i );
		core::pose::addVirtualResAsRoot( *posecopy );

		// Setup rigid-body movemap now! (we need to know root jump number)
		core::kinematics::MoveMapOP rbmm = utility::pointer::make_shared< core::kinematics::MoveMap >();
		rbmm->set_bb( false ); rbmm->set_chi( false ); rbmm->set_jump( false );
		int root = posecopy->fold_tree().root();
		utility::vector1< core::kinematics::Edge > root_edges = posecopy->fold_tree().get_outgoing_edges(root);
		for ( core::Size i=1; i<=root_edges.size(); ++i ) rbmm->set_jump ( root_edges[i].label() , true );

		// Setup rigid-body min now!
		protocols::minimization_packing::MinMoverOP rbmin(
			utility::pointer::make_shared< protocols::minimization_packing::MinMover >(
			rbmm, scorefxn_refine_rb, "lbfgs_armijo_nonmonotone", 1e-5, true ) );
		rbmin->max_iter(200); // make a parameter?

		core::Real const scoreb = (*scorefxn_dens)(*posecopy);
		core::Real scorei=0;

		// rbmin/pack/fullmin
		if ( do_refine_ ) {
			rbmin->apply( *posecopy );

			if ( posecopy->is_centroid() ) {
				to_all_atom.apply( *posecopy );
			}

			for ( core::Size i=1; i<=ncyc_; ++i ) {
				packer->apply( *posecopy );
				scorei = (*scorefxn_dens)(*posecopy);

				if ( min_backbone_ ) {
					//scorefxn_refine->show( *posecopy );
					bbmin->apply( *posecopy );
					//scorefxn_refine->show( *posecopy );
				} else {
					rbmin->apply( *posecopy );
				}
			}
		}

		// rescore (expensive sf)
		core::Real const scoref = (*scorefxn_dens)(*posecopy);

		core::Real rms=0.0;
		if ( native_ ) {
			rms = get_rms( *posecopy, *native_, symminfo_ );
		}

		TR << "[" << ntotal-results_in.size() << "/" << ntotal << "] " << scoreb << " : " << scorei << " : " << scoref << "  rms=" << rms << std::endl;

		// store
		results_out.add_element( RefinementResult( -scoref, -scoreb, sol_i.score_, posecopy ) );
	}
	TR << std::endl;
}


void
DockFragmentsIntoDensityMover::apply( core::pose::Pose & pose) {
	// call multipose mover
	utility::vector1< core::pose::PoseOP > temp;
	temp.push_back( utility::pointer::make_shared< core::pose::Pose >(pose) );
	apply_multi( temp );
	pose = *(temp[1]);
}


void
DockFragmentsIntoDensityMover::apply_multi( utility::vector1< core::pose::PoseOP > & poses) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose;

	std::string base_name = tag_;
	if ( base_name.size() == 0 ) {
		// get output name ... assumes jd2 ... anything that doesn't use jd2 must call setTag
		base_name = protocols::jd2::current_input_tag();
		utility::vector1< std::string > temp_out_names= utility::split( base_name );
		utility::file::FileName out_name = utility::file::combine_names( temp_out_names );
		base_name = out_name.base();
	}

	RBfitResultDB results_filtered(cluster_oversample_*topNfilter_);    // oversample for clustering
	RefinementResultDB results_refine(cluster_oversample_*topNfinal_);  // oversample for clustering

	if ( passthrough_ ) {
		poses.clear();
		poses.push_back( native_->clone() );

		RBfitResultDB results(poses.size());
		for ( int i=1; i<=(int)poses.size(); ++i ) {
			results.add_element(
				RBfitResult(1, -1000, numeric::xyzMatrix<core::Real>::identity(), numeric::xyzVector<core::Real>(0,0,0), numeric::xyzVector<core::Real>(0,0,0))
			);
			do_refinement ( poses, results, results_refine );
		}

		// force score to -1000 <---- HACK
		RefinementResult sol_i = results_refine.pop();
		sol_i.score_ = -1000;
		results_refine.add_element(sol_i);
	} else {
		// for each input fragment...
		for ( int i=1; i<=(int)poses.size(); ++i ) {
			core::pose::PoseOP pose_i = poses[i];

			// set nRsteps
			TR << " *** Input " << i << " of " << poses.size() << " ***" << std::endl;
			numeric::xyzVector< core::Real > com;

			// 1: get points to sample
			//  NOTE: this sets nRsteps_!
			points_to_search_ = select_density_points(
				*pose_i,
				SelectDensityPointsOptions{
				B_, nRsteps_, gridStep_, topNtrans_,
				delR_, fragDens_, point_radius_, laplacian_offset_,
				center_on_middle_ca_, convolute_single_residue_,
				symminfo_, get_multi_native(), get_multi_native_com(), get_multi_native_mca()
				},
				nRsteps_
			);
			TR << "Searching a max radius of " << nRsteps_*delR_ << " (" << nRsteps_ << " steps)" << std::endl;
			TR << "Selected " << points_to_search_.size() << " translations to search" << std::endl;

			// 2: do the grid search
			//    top N per position
			RBfitResultDB results(cluster_oversample_*topNfilter_);

			density_grid_search(
				i,*pose_i, results,
				points_to_search_,
				DensityGridSearchOptions{
				B_, nRsteps_, max_rot_per_trans_, 1 /*point_search_start*/, points_to_search_.size(),
				delR_, cluster_radius_, laplacian_offset_, 5.0 /*result cluster rms_cutoff*/,
				center_on_middle_ca_, convolute_single_residue_, false, /* <-include_distance */
				"" /*output_fn*/,
				symminfo_, get_multi_native(), get_multi_native_com(), get_multi_native_mca()
				}
			);

			TR << "Have " << results.size() << " solutions with score >= " << results.top().score_ << std::endl;
			if ( native_ ) print_best_rms( *pose_i, results );

			// 3: rescore
			if ( cluster_radius_>0.0 ) {
				do_filter( poses, results, /*rescore=*/ false, symminfo_, cluster_radius_ ); // rescore=false
			}

			while ( results.size() > topNfilter_ ) results.pop();

			TR << "Have " << results.size() << " solutions after clustering" << std::endl;
			if ( native_ ) print_best_rms( *pose_i, results );

			// 4: combine
			while ( results.size() > 0 ) {
				results_filtered.add_element( results.pop() );
			}
		}
	}

	TR << "Final output" << std::endl;

	// now cluster the aggregate pool (if necessary)
	if ( poses.size()>1 ) {
		TR << "Have " << results_filtered.size() << " total solutions before clustering" << std::endl;
		do_filter( poses, results_filtered, /*rescore=*/ false, symminfo_, cluster_radius_ ); // rescore=false
	}

	TR << "Have " << results_filtered.size() << " total solutions after clustering" << std::endl;

	// 5: refinement
	do_refinement ( poses, results_filtered, results_refine );
	TR << "Have " << results_refine.size() << " total solutions after refinement" << std::endl;

	// cluster a final time (in case structures minimized into the same energy well)
	// also apply local rescoring
	cluster_RefinementDB(results_refine, symminfo_, cluster_radius_, 0);

	// only dump the top N
	while ( results_refine.size() > topNfinal_ ) { results_refine.pop(); }

	if ( silent_.size() > 0 ) {
		dump_RefinementDB_to_silent(
			results_refine,
			silent_,
			base_name,
			"A",
			/*centroid_output*/ true,
			/*appent_to_outfile*/ false,
			get_multi_native(),
			symminfo_,
			/* legacy_rms */ true);
	} else {
		while ( results_refine.size() > 0 ) {
			// tag
			core::io::RemarkInfo remark;
			std::ostringstream oss;
			RefinementResult sol_i = results_refine.pop();
			core::pose::PoseOP const posecopy = sol_i.pose_->clone();

			if ( posecopy->pdb_info() ) {
				oss << "RANK = " << results_refine.size()+1;

				remark.num = 1; remark.value = oss.str();
				posecopy->pdb_info()->remarks().push_back( remark );

				oss.str(""); oss.clear();
				oss << "SCORE = " << sol_i.score_;
				remark.num = 1; remark.value = oss.str();
				posecopy->pdb_info()->remarks().push_back( remark );
			}

			std::string const outname = base_name+"_"+ObjexxFCL::right_string_of( results_refine.size()+1, 6, '0' )+".pdb";
			posecopy->dump_pdb( outname );
		}
	}
}

}
}
