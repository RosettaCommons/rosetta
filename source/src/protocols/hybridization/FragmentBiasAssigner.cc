// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief take out the "compute_frag_bias()" from CartesianSampler
/// @author Ray Wang wangyr@uw.edu
//
#include <protocols/hybridization/FragmentBiasAssigner.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/EnergyGraph.hh>

#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragSet.hh>

// symmetry
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/electron_density/util.hh>

#include <numeric/random/random.hh>
#include <numeric/random/WeightedSampler.hh>
#include <numeric/zscores.hh>

#include <utility/vector1.hh>

#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>

#include <basic/Tracer.hh>

namespace protocols {
namespace hybridization {

using namespace core;
using namespace core::scoring::electron_density;

static basic::Tracer fragbias_tr( "protocols.hybridization.FragmentBiasAssigner" );

FragmentBiasAssigner::
FragmentBiasAssigner(
	pose::Pose &pose
){
	init( pose );
}

// initialize values for fragmentProbs_
void
FragmentBiasAssigner::
init(
	pose::Pose &pose
){
	fragProbs_assigned_ = false;
	// get nres_, accounting for symmetry/vrt/ligands
	nres_ = pose.size();
	fragbias_tr.Trace << "init(): pose.size(): " << nres_ << std::endl;

	// symmetry
	symminfo_=nullptr;
	n_symm_subunit_=1;

	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		auto const & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation const &> ( pose.conformation()) );
		symminfo_ = SymmConf.Symmetry_Info();
		nres_ = symminfo_->num_independent_residues();
		n_symm_subunit_  = symminfo_->score_multiply_factor();
	}
	while ( !pose.residue(nres_).is_protein() ) nres_--;

	fragmentProbs_.resize(nres_, 1e-6);
	fragbias_tr.Trace << "init(): symmetrical_pose " << core::pose::symmetry::is_symmetric( pose )<< std::endl;
	fragbias_tr.Trace << "init(): n residues " << nres_ << std::endl;
	fragbias_tr.Trace << "init(): n_symm_subunit_ " << n_symm_subunit_ << std::endl;

	score_threshold_ = 9999;
	cumulative_ = false;

	wdw_to_freeze_ = 0;
	rsd_wdw_size_ = 3;
	cumulative_ = false;
}


// take the fragmentProbs_ assigned by methods selected, and compile them to frames
void
FragmentBiasAssigner::
compute_frag_bias(
	utility::vector1<numeric::random::WeightedSampler> &frag_bias, // output
	pose::Pose &pose,
	utility::vector1<core::fragment::FragSetOP> fragment_sets
){
	runtime_assert( fragProbs_assigned_ );

	// clean up the vector
	frag_bias.resize( fragment_sets.size(), 0.0 );

	// for each fragment size, smooth over the fragment window
	//   - handle mapping from frame->seqpos
	//   - don't allow any insertions that cross cuts
	for ( Size i_frag_set = 1; i_frag_set<=fragment_sets.size(); ++i_frag_set ) {
		utility::vector1< core::Real > frame_weights( pose.size(), 0.0 );

		for ( Size i_frame = 1; i_frame <= fragment_sets[i_frag_set]->nr_frames(); ++i_frame ) {
			core::fragment::ConstFrameIterator frame_it = fragment_sets[i_frag_set]->begin(); // first frame of the fragment library
			std::advance(frame_it, i_frame-1);  // point frame_it to the i_frame of the library
			core::Size seqpos_start = (*frame_it)->start();  // find starting and ending residue seqpos of the inserted fragment
			core::Size seqpos_end   = (*frame_it)->end();

			frame_weights[seqpos_start] = 0;
			for ( auto i_pos = (int)seqpos_start; i_pos<=(int)seqpos_end; ++i_pos ) {
				frame_weights[seqpos_start] += fragmentProbs_[i_pos];
			}
			frame_weights[seqpos_start] /= (seqpos_end-seqpos_start+1);

			// don't allow insertions at cut points
			for ( auto i_pos = (int)seqpos_start; i_pos<=(int)seqpos_end-1; ++i_pos ) {
				if ( pose.fold_tree().is_cutpoint(i_pos) ) {
					for ( int i_wdw = 0; i_wdw<wdw_to_freeze_; ++i_wdw ) {
						frame_weights[seqpos_start-i_wdw] = 0.0; // don't allow insertions at a frame where there are cut points downstream
						fragbias_tr.Trace << "chainbreak at: " <<  i_pos << " seqpos_start: " << seqpos_start << " residue_to_freeze: " << seqpos_start-i_wdw << " wdw_to_freeze_: " << wdw_to_freeze_ << std::endl;
					}
				}
			} // assign
		} // each frame in a fragment set
		frag_bias[i_frag_set].weights(frame_weights);

		//////////////////////////////////////////////
		/////////////// for debug only ///////////////
		// report the probability for each frame in a fragment set
		for ( Size i_frame = 1; i_frame <= fragment_sets[i_frag_set]->nr_frames(); ++i_frame ) {
			core::fragment::ConstFrameIterator frame_it = fragment_sets[i_frag_set]->begin(); // first frame of the fragment library
			std::advance(frame_it, i_frame-1);  // point frame_it to the i_frame of the library
			core::Size seqpos_start = (*frame_it)->start();  // find starting and ending residue seqpos of the inserted fragment
			core::Size seqpos_end   = (*frame_it)->end();

			for ( auto i_pos = (int)seqpos_start; i_pos<=(int)seqpos_end-1; ++i_pos ) { //
				fragbias_tr.Trace << "Prob_dens( " << i_frame << " " << seqpos_start << " ) = " << frame_weights[seqpos_start] << std::endl;
			}

			fragbias_tr.Trace << "i_frame: " << i_frame << " seqpos_start " << seqpos_start <<  " seqpos_end " << seqpos_end  <<  " frame_weights: " <<  frame_weights[seqpos_start] << std::endl;

		} // frame /////////////// for debug only ///////////////
	} // each fragment set

}


////////////////////////////////////////////////////
// a generic method to calculate per-residue score
// warning!! for hb terms it might not be accurate
void
FragmentBiasAssigner::
cal_perrsd_score(
	pose::Pose &pose,
	scoring::ScoreType const &score_type,
	utility::vector1<core::Real> &perrsd_score,
	Real weight
){

	runtime_assert( nres_ > 0 );
	runtime_assert( perrsd_score.size() == nres_ );

	core::scoring::ScoreFunctionOP myscore( new core::scoring::ScoreFunction() );
	myscore->set_weight( score_type, weight );

	if ( core::pose::symmetry::is_symmetric(pose) ) {
		myscore = core::scoring::symmetry::symmetrize_scorefunction(*myscore);
	}

	(*myscore)(pose);
	//myscore->show_line( fragbias_tr, pose ); fragbias_tr << std::endl;

	for ( int r=1; r<=(int)nres_; ++r ) {
		int rsrc = r;
		if ( symminfo_ && !symminfo_->bb_is_independent(r) ) { // for helical symmetry targets - the main chain for scoring might not be chain A
			rsrc = symminfo_->bb_follows(r);
		}
		perrsd_score[r] = weight*(pose.energies().residue_total_energies(rsrc)[ score_type ]/n_symm_subunit_);
	}
}


// rsd window size is controlled through set_rsd_wdw_to_assign_prob(rsd_wdw_size_)
void
FragmentBiasAssigner::
assign_prob_with_rsd_wdw(
	int rsn
){
	int offset=((int)rsd_wdw_size_-1)/2;
	runtime_assert( offset >= 0 );

	int start_rsn = rsn-offset;
	if ( start_rsn < 1 ) start_rsn = 1;

	Size end_rsn = rsn+offset;
	if ( end_rsn > nres_ ) end_rsn = nres_;

	for ( int i=start_rsn; i<=(int)end_rsn; ++i ) {
		fragmentProbs_[i] = 1.0;
		fragbias_tr.Trace << "rsn(rsn): " << rsn << " fragProb[" << i << "]: " << fragmentProbs_[i] << " rsd_window_size: " << rsd_wdw_size_ << std::endl;
	}

}


////////////////////////////////////////////////////
//
void
FragmentBiasAssigner::
automode(
	pose::Pose &pose,
	Real score_cut /*-0.5 is controlled through CartesianSampler*/
){
	fragmentProbs_.resize(nres_, 1e-6);

	// based on parameter scanning to predict regions to refine:
	// 0.45 0.05 0.15 0.35
	// density:9, density_nbrzscore:1, rama:3, geometry:7
	//std::map< Size, Real > per_rsd_dens, per_rsd_nbrdens, per_rsd_rama, per_rsd_geometry;
	std::map< Size, Real > zscore_dens, zscore_nbrdens, zscore_rama, zscore_geometry;
	std::map< Size, Real > scores;

	density_nbr( pose ); // perrsd_dens_ and perrsd_densnbr_
	rama( pose ); // perrsd_rama_
	geometry( pose ); // perrsd_geometry_

	numeric::calc_zscore( perrsd_dens_,     zscore_dens           );
	numeric::calc_zscore( perrsd_nbrdens_,  zscore_nbrdens        );
	numeric::calc_zscore( perrsd_rama_,     zscore_rama,     true );
	numeric::calc_zscore( perrsd_geometry_, zscore_geometry, true );

	//all_scores.resize(nres_, 0.0);
	for ( auto r : perrsd_dens_ ) {

		Real score =  0.45*zscore_dens[r.first]
			+ 0.05*zscore_nbrdens[r.first]
			+ 0.15*zscore_rama[r.first]
			+ 0.35*zscore_geometry[r.first];

		if ( score < score_cut ) assign_prob_with_rsd_wdw(r.first);

		//fragbias_tr << "rsn: " << r << " fragProb: " << fragmentProbs_[r] << " score: " << score << std::endl;
	}
}


////////////////////////////////////////////////////
//
void
FragmentBiasAssigner::
automode_scores(
	pose::Pose &pose,
	std::map< Size, Real > & scores
){
	// based on parameter scanning to predict regions to refine:
	// 0.45 0.05 0.15 0.35
	// density:9, density_nbrzscore:1, rama:3, geometry:7
	//std::map< Size, Real > per_rsd_dens, per_rsd_nbrdens, per_rsd_rama, per_rsd_geometry;
	std::map< Size, Real > zscore_dens, zscore_nbrdens, zscore_rama, zscore_geometry;

	density_nbr( pose ); // perrsd_dens_ and perrsd_densnbr_
	rama( pose ); // perrsd_rama_
	geometry( pose ); // perrsd_geometry_

	numeric::calc_zscore( perrsd_dens_,     zscore_dens           );
	numeric::calc_zscore( perrsd_nbrdens_,  zscore_nbrdens        );
	numeric::calc_zscore( perrsd_rama_,     zscore_rama,     true );
	numeric::calc_zscore( perrsd_geometry_, zscore_geometry, true );

	fragbias_tr.Trace << "Size: " << zscore_dens.size() << std::endl;
	//all_scores.resize(nres_, 0.0);
	for ( auto r : perrsd_dens_ ) {
		
		fragbias_tr.Trace << zscore_dens[r.first] << " " << zscore_nbrdens[r.first] << " " << zscore_rama[r.first] << " " << zscore_geometry[r.first] << std::endl;
		Real score =  0.45* zscore_dens[r.first]
			+ 0.05*zscore_nbrdens[r.first]
			+ 0.15*zscore_rama[r.first]
			+ 0.35*zscore_geometry[r.first];
		
		scores[r.first] = score;
	}
}

// assign 1 or 0
// add values into
void
FragmentBiasAssigner::assign_fragprobs(
	std::map< core::Size, core::Real > const & perrsd_score,
	Real threshold )
{
	for ( auto score_pair : perrsd_score) {
		if ( score_pair.second >= threshold ) {
			if ( cumulative_ ) {
				fragmentProbs_[ score_pair.first ] += 1.0; // no residue window size control here
			} else {
				assign_prob_with_rsd_wdw(score_pair.first);
			}
		}
		fragbias_tr.Trace << "rsn: " << score_pair.first << " fragProb: " << fragmentProbs_[ score_pair.first ] << " score: " << score_pair.second << " " << threshold << std::endl;
	}
}


void
FragmentBiasAssigner::
include_residues(
	std::set< core::Size > residues_to_include
	//int window_size
){
	for ( int r=1; r<=(int)nres_; ++r ) {
		if ( residues_to_include.find(r) != residues_to_include.end() ) {
			assign_prob_with_rsd_wdw(r);
		}
	}
}


void
FragmentBiasAssigner::
exclude_residues(
	std::set< core::Size > residues_to_exclude
	//int window_size
){
	for ( int r=1; r<=(int)nres_; ++r ) {
		if ( residues_to_exclude.find(r) != residues_to_exclude.end() ) {
			fragmentProbs_[r] = 0.0; // should I add window here? rw: probably no need
		}
	}
}


void
FragmentBiasAssigner::
uniform(
){
	fragProbs_assigned_=true;

	for ( int r=1; r<=(int)nres_; ++r ) {
		fragmentProbs_[r] = 1.0;
		fragbias_tr.Debug << "Prob_dens_uniform( " << r << " ) = " << 1.0 << std::endl;
	}
	fragProbs_assigned_=true;
}


void
FragmentBiasAssigner::
user(
	std::set<core::Size> user_pos,
	protocols::loops::LoopsOP loops
){
	fragProbs_assigned_ = true;

	// user defined segments to rebuild
	runtime_assert( user_pos.size()>0 || (loops && !loops->empty()) );

	for ( int r=1; r<=(int)nres_; ++r ) {
		//fragmentProbs_[r] = 0.0;
		if ( user_pos.find(r) != user_pos.end() ) assign_prob_with_rsd_wdw(r);
		if ( loops && loops->is_loop_residue(r) ) fragmentProbs_[r] = 1.0;
		fragbias_tr.Debug << "Prob_dens_user( " << r << " ) = " << fragmentProbs_[r] << std::endl;
	}
}

// fills in two vectors of member variable: perrsd_dens_ and perrsd_nbrdens_
void
FragmentBiasAssigner::
density_nbr(
	pose::Pose &pose
){
	using namespace core::scoring;
	// rescore pose
	fragProbs_assigned_=true;
	
	// clean and init the containers
	perrsd_dens_.clear();
	perrsd_nbrdens_.clear();
	for (Size i= 1; i <=nres_; ++i) {
		perrsd_dens_[i] = 0.0;
		perrsd_nbrdens_[i] = 0.0;
	}
	
	calculate_density_nbr( pose, perrsd_dens_, perrsd_nbrdens_, symminfo_ );
	
	for (auto score_pair : perrsd_nbrdens_){
		Real i_perrsd_dens_zscore = score_pair.second;
		Real i_nbrdens_zscore = perrsd_dens_[ score_pair.first ];
		
		if ( i_perrsd_dens_zscore < -1.0 ) { // the rscc for this residue is being considered worse overall
			if ( i_nbrdens_zscore < -2.0 ) {  // the rscc for this residue is worse than its neighbors
				if ( cumulative_ ) {
					fragmentProbs_[score_pair.first] += 1;
				} else {
					fragmentProbs_[score_pair.first] = 1;
				}
			} // i_perrsd_zscore
		} // i_nbr_zscore

		fragbias_tr.Trace << "fragmentProbs_: " << score_pair.first << " rsd: " << score_pair.first << " prob: " << fragmentProbs_[score_pair.first] << std::endl;
	} // for i in range(pose.size())
}


void
FragmentBiasAssigner::
density(
	pose::Pose &pose
){
	fragProbs_assigned_=true;

	// find segments with worst agreement to the density
	core::scoring::electron_density::ElectronDensity &edm = core::scoring::electron_density::getDensityMap();
	edm.setScoreWindowContext( true );
	edm.setWindow( 3 );  // smoother to use 3-res window

	// score the pose
	core::scoring::ScoreFunctionOP myscore( new core::scoring::ScoreFunction() );
	myscore->set_weight( core::scoring::elec_dens_window, 1.0 );

	// make sure interaction graph gets computed
	if ( pose.is_fullatom() ) {
		myscore->set_weight( core::scoring::fa_rep, 1.0 );
	} else {
		myscore->set_weight( core::scoring::vdw, 1.0 );
	}

	if ( core::pose::symmetry::is_symmetric(pose) ) {
		myscore = core::scoring::symmetry::symmetrize_scorefunction(*myscore);
	}

	(*myscore)(pose);

	perrsd_dens_.clear();
	for (Size i= 1; i <=nres_; ++i) perrsd_dens_[i] = 0.0;
	
	core::Real CCsum=0, CCsum2=0;

	// turn score_symm_complex flag on; otherwise it won't see symmetry, and will return 0
	for ( int r=1; r<=(int)nres_; ++r ) {

		perrsd_dens_[r] = core::scoring::electron_density::getDensityMap().matchRes( r, pose.residue(r), pose, symminfo_, false);
		CCsum += perrsd_dens_[r];
		CCsum2 += perrsd_dens_[r]*perrsd_dens_[r];
	}
	CCsum /= nres_;
	CCsum2 = sqrt( CCsum2/nres_-CCsum*CCsum );

	for ( auto r : perrsd_dens_) {
		if ( r.second < 0.6 ) {
			fragmentProbs_[r.first] = 1.0;
		} else if ( perrsd_dens_[r.first]<0.8 ) {
			fragmentProbs_[r.first] = 0.25;
		} else {
			fragmentProbs_[r.first] = 0.01;
		}
		fragbias_tr.Trace << "residue " << r.first << " rscc: " << r.second << " weight: " <<fragmentProbs_[r.first] << std::endl;
	}
}


void
FragmentBiasAssigner::
geometry(
	pose::Pose &pose,
	Real weight  /*1.0*/
){
	fragProbs_assigned_=true;
	if ( score_threshold_ == 9999 ) set_score_threshold( 0.6 );

	// clean the container
	perrsd_geometry_.clear();
	for (Size i= 1; i <=nres_; ++i) perrsd_geometry_[i] = 0.0;
	
	calculate_geometry( pose, perrsd_geometry_, n_symm_subunit_, weight);

	assign_fragprobs( perrsd_geometry_,score_threshold_ );
}


void
FragmentBiasAssigner::
rama(
	pose::Pose &pose,
	Real weight /*=0.2*/
){
	fragProbs_assigned_=true;
	if ( score_threshold_ == 9999 ) set_score_threshold( 0.7 );

	// clean the container
	perrsd_rama_.clear();
	for (Size i= 1; i <=nres_; ++i) perrsd_rama_[i] = 0.0;
	
	calculate_rama( pose, perrsd_rama_, n_symm_subunit_, weight);
	
	assign_fragprobs( perrsd_rama_, score_threshold_ );
}

void
FragmentBiasAssigner::
bfactors(
	pose::Pose &pose
){
	fragProbs_assigned_=true;

	// find segments with highest bfactors
	core::Real Btemp=25;  // no idea what value makes sense here
	// with Btemp = 25, a B=100 is ~54 times more likely to be sampled than B=0
	runtime_assert( pose.pdb_info() != nullptr );
	for ( int r=1; r<=(int)nres_; ++r ) {
		core::Real Bsum=0;
		core::Size nbb = pose.residue_type(r).last_backbone_atom();
		for ( core::Size atm=1; atm<=nbb; ++atm ) {
			Bsum += pose.pdb_info()->temperature( r, atm );
		}
		Bsum /= nbb;
		fragmentProbs_[r] = exp( Bsum/Btemp );
		fragbias_tr.Trace << "Prob_dens_bfact( " << r << " ) = " << fragmentProbs_[r] << " ; B=" << Bsum << std::endl;
	}
}


void
FragmentBiasAssigner::
chainbreak(
	pose::Pose &pose
){
	fragProbs_assigned_=true;

	for ( int r=1; r<(int)nres_; ++r ) {
		if ( !pose.residue_type(r).is_protein() ) continue;
		if ( pose.fold_tree().is_cutpoint(r+1) ) continue;

		numeric::xyzVector< core::Real > c0 , n1;
		c0 = pose.residue(r).atom(" C  ").xyz();
		n1 = pose.residue(r+1).atom(" N  ").xyz();
		core::Real d2 = c0.distance( n1 );
		if ( (d2-1.328685)*(d2-1.328685) > 0.1*0.1 ) {
			fragmentProbs_[r] = 1.0;
		} else {
			fragmentProbs_[r] = 0.001;
		}
		fragbias_tr.Debug << "Prob_dens_cb( " << r << " ) = " << fragmentProbs_[r] << std::endl;
	}
}


void
FragmentBiasAssigner::
fragbias_reporter(
	pose::Pose &/*pose*/
){
	for ( int r=1; r<=(int)nres_; ++r ) {
		//if (!pose.residue_type(r).is_protein()) continue;
		//if (pose.fold_tree().is_cutpoint(r+1)) continue;

		fragbias_tr.Trace << "rsn: " << r << " prob: " << fragmentProbs_[r] << std::endl;
	}
}


} // hybridization
} // protocol
