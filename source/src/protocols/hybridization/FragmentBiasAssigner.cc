// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

#include <numeric/random/random.hh>
#include <numeric/random/WeightedSampler.hh>

#include <utility/vector1.hh>

#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>

#include <basic/Tracer.hh>

namespace protocols {
namespace hybridization {

using namespace core;

static THREAD_LOCAL basic::Tracer fragbias_tr( "protocols.hybridization.FragmentBiasAssigner" );

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
	nres_ = pose.total_residue();
	fragbias_tr << "init(): pose.total_residue(): " << nres_ << std::endl;

	// symmetry
	symminfo_=NULL;
	n_symm_subunit_=1;

	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		core::conformation::symmetry::SymmetricConformation const & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation const &> ( pose.conformation()) );
		symminfo_ = SymmConf.Symmetry_Info();
		nres_ = symminfo_->num_independent_residues();
		n_symm_subunit_  = symminfo_->score_multiply_factor();
	}
	while ( !pose.residue(nres_).is_protein() ) nres_--;

	fragmentProbs_.resize(nres_, 0.0);
	fragbias_tr << "init(): symmetrical_pose " << core::pose::symmetry::is_symmetric( pose )<< std::endl;
	fragbias_tr << "init(): n residues " << nres_ << std::endl;
	fragbias_tr << "init(): n_symm_subunit_ " << n_symm_subunit_ << std::endl;
}


// take the fragmentProbs_ assigned by methods selected, and compile them to frames
void
FragmentBiasAssigner::
compute_frag_bias(
	utility::vector1<numeric::random::WeightedSampler> &frag_bias, // output
	pose::Pose &pose,
	utility::vector1<core::fragment::FragSetOP> fragment_sets
){
	fragbias_tr << "compute frag bias" << std::endl;
	runtime_assert( fragProbs_assigned_ );

	// clean up the vector
	frag_bias.resize( fragment_sets.size(), 0.0 );

	// for each fragment size, smooth over the fragment window
	//   - handle mapping from frame->seqpos
	//   - don't allow any insertions that cross cuts
	for ( Size i_frag_set = 1; i_frag_set<=fragment_sets.size(); ++i_frag_set ) {
		utility::vector1< core::Real > frame_weights( pose.total_residue(), 0.0 );

		for ( Size i_frame = 1; i_frame <= fragment_sets[i_frag_set]->nr_frames(); ++i_frame ) {
			core::fragment::ConstFrameIterator frame_it = fragment_sets[i_frag_set]->begin(); // first frame of the fragment library
			std::advance(frame_it, i_frame-1);  // point frame_it to the i_frame of the library
			core::Size seqpos_start = (*frame_it)->start();  // find starting and ending residue seqpos of the inserted fragment
			core::Size seqpos_end   = (*frame_it)->end();

			frame_weights[seqpos_start] = 0;
			for ( int i_pos = (int)seqpos_start; i_pos<=(int)seqpos_end; ++i_pos ) {
				frame_weights[seqpos_start] += fragmentProbs_[i_pos];
			}
			frame_weights[seqpos_start] /= (seqpos_end-seqpos_start+1);

			// don't allow insertions at cut points
			for ( int i_pos = (int)seqpos_start; i_pos<=(int)seqpos_end-1; ++i_pos ) {
				if ( pose.fold_tree().is_cutpoint(i_pos) ) {
					for ( int i_wdw = 0; i_wdw<=wdw_to_freeze_; ++i_wdw ) {
						frame_weights[seqpos_start-i_wdw] = 0.0; // don't allow insertions at a frame where there are cut points downstream
						fragbias_tr << "chainbreak at: " <<  i_pos << " seqpos_start: " << seqpos_start << " residue_to_freeze: " << seqpos_start-i_wdw << " wdw_to_freeze_: " << wdw_to_freeze_ << std::endl;
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

			for ( int i_pos = (int)seqpos_start; i_pos<=(int)seqpos_end-1; ++i_pos ) { //
				fragbias_tr << "Prob_dens( " << i_frame << " " << seqpos_start << " ) = " << frame_weights[seqpos_start] << std::endl;
			}

			fragbias_tr << "i_frame: " << i_frame << " seqpos_start " << seqpos_start <<  " seqpos_end " << seqpos_end  <<  " frame_weights: " <<  frame_weights[seqpos_start] << std::endl;

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
	myscore->show_line( fragbias_tr, pose ); fragbias_tr << std::endl;

	for ( int r=1; r<=(int)nres_; ++r ) {
		int rsrc = r;
		if ( symminfo_ && !symminfo_->bb_is_independent(r) ) { // for helical symmetry targets - the main chain for scoring might not be chain A
			rsrc = symminfo_->bb_follows(r);
		}
		perrsd_score[r] = weight*(pose.energies().residue_total_energies(rsrc)[ score_type ]/n_symm_subunit_);
	}
}


void
FragmentBiasAssigner::
cal_zscore(
	utility::vector1<core::Real> const &input_v,
	utility::vector1<core::Real> &zscore_v,
	bool negating
){
	//runtime_assert( input_v.size() == nres_ ); // do I need to have this here?
	Real sum=0, sq_sum=0;

	for ( Size i=1; i<=nres_; ++i ) {
		sum    += input_v[i];
		sq_sum += input_v[i]*input_v[i];
	}
	Real mean  = sum/nres_;
	Real stdev = std::sqrt( sq_sum/nres_ - mean*mean );

	for ( Size i=1; i<=nres_; ++i ) {

		Real i_zscore =  (input_v[i]-mean)/stdev;

		if ( negating ) i_zscore = -1*i_zscore;

		zscore_v[i] = i_zscore;
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
		fragbias_tr << "rsn(rsn): " << rsn << " fragProb[" << i << "]: " << fragmentProbs_[i] << " rsd_window_size: " << rsd_wdw_size_ << std::endl;
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
	fragbias_tr << "you are using automode!" << std::endl;
	fragmentProbs_.resize(nres_, 0.0);

	// based on parameter scanning to predict regions to refine:
	// 0.45 0.05 0.15 0.35
	// density:9, density_nbrzscore:1, rama:3, geometry:7
	utility::vector1<core::Real> zscore_dens, zscore_nbrdens, zscore_rama, zscore_geometry;

	// to catch some outliers
	zscore_dens.resize(nres_,0.0);
	zscore_nbrdens.resize(nres_,0.0);
	zscore_rama.resize(nres_,0.0);
	zscore_geometry.resize(nres_,0.0);

	density_nbr( pose ); // perrsd_dens_ and perrsd_densnbr_
	rama( pose ); // perrsd_rama_
	geometry( pose ); // perrsd_geometry_

	cal_zscore( perrsd_dens_,     zscore_dens           );
	cal_zscore( perrsd_nbrdens_,  zscore_nbrdens        );
	cal_zscore( perrsd_rama_,     zscore_rama,     true );
	cal_zscore( perrsd_geometry_, zscore_geometry, true );

	//all_scores.resize(nres_, 0.0);
	for ( int r=1; r<=(int)nres_; ++r ) {

		Real score =  0.45*zscore_dens[r]
			+ 0.05*zscore_nbrdens[r]
			+ 0.15*zscore_rama[r]
			+ 0.35*zscore_geometry[r];

		if ( score < score_cut ) assign_prob_with_rsd_wdw(r);

		//fragbias_tr << "rsn: " << r << " fragProb: " << fragmentProbs_[r] << " score: " << score << std::endl;
	}
}

// assign 1 or 0
// add values into
void
FragmentBiasAssigner::
assign_fragprobs(
	utility::vector1<core::Real> const &perrsd_score,
	Real threshold
){
	for ( int r=1; r<=(int)nres_; ++r ) {
		if ( perrsd_score[r] >= threshold ) {
			if ( cumulative_ ) {
				fragmentProbs_[r] += 1.0; // no residue window size control here
			} else {
				assign_prob_with_rsd_wdw(r);
			}
		}
		fragbias_tr << "rsn: " << r << " fragProb: " << fragmentProbs_[r] << " score: " << perrsd_score[r] << std::endl;
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
	fragbias_tr << "uniform method is selected" << std::endl;
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
	fragbias_tr << "user is chosen" << std::endl;
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
	fragbias_tr << "density_nbr is chosen" << std::endl;
	fragProbs_assigned_=true;

	core::scoring::electron_density::ElectronDensity &edm = core::scoring::electron_density::getDensityMap();
	edm.setScoreWindowContext( true );
	edm.setWindow( 3 );  // smoother to use 3-res window

	// score the pose
	core::scoring::ScoreFunctionOP myscore( new core::scoring::ScoreFunction() );
	myscore->set_weight( core::scoring::elec_dens_window, 1.0 );

	if ( pose.is_fullatom() ) {
		myscore->set_weight( core::scoring::fa_rep, 10e-30 );
	} else {
		myscore->set_weight( core::scoring::vdw, 10e-30 );
	}

	if ( core::pose::symmetry::is_symmetric(pose) ) {
		myscore = core::scoring::symmetry::symmetrize_scorefunction(*myscore);
	}

	fragbias_tr << "scoring the pose" << std::endl;
	(*myscore)(pose);
	myscore->show_line( fragbias_tr, pose ); fragbias_tr << std::endl;

	///////////////////////////////////////////////////////
	// get density correlation score from pose.total_residue(), rather than nres_
	// get zscore for real-space density correlation scores
	perrsd_dens_.resize(nres_, 0.0);
	core::Real rscc_sum=0, sq_rscc_sum=0;

	for ( Size r=1; r<=pose.total_residue(); ++r ) { // loop over the entire pose
		if ( pose.residue(r).aa() == core::chemical::aa_vrt ) continue;
		if ( symminfo_ && !symminfo_->bb_is_independent( r ) ) continue; // only the main chain gets selected

		Real dens_rscc = core::scoring::electron_density::getDensityMap().matchRes( r , pose.residue(r), pose, symminfo_ , false);
		Size asymm_num_r = ((r-1)%nres_)+1;
		perrsd_dens_[asymm_num_r] = dens_rscc;
		rscc_sum    += perrsd_dens_[asymm_num_r];
		sq_rscc_sum += perrsd_dens_[asymm_num_r]*perrsd_dens_[asymm_num_r];

		fragbias_tr << "res: " << asymm_num_r << " symmetric num: " << r << " dens_rscc: " << dens_rscc << std::endl;
	}
	// get mean and stdev for density rscc
	Real perrsd_dens_mean  = rscc_sum/nres_;
	Real perrsd_dens_stdev = std::sqrt( sq_rscc_sum/nres_ - perrsd_dens_mean*perrsd_dens_mean );


	///////////////////////////////////////////////////////
	// for each residue, get neighbors from energy graph,
	// and calculate density-zscore from neighbors
	EnergyGraph const & energy_graph( pose.energies().energy_graph() );
	perrsd_nbrdens_.resize(nres_, 0.0);

	for ( Size i=1; i<=pose.total_residue(); ++i ) {
		if ( pose.residue(i).aa() == core::chemical::aa_vrt ) continue; // prevent virtual
		if ( symminfo_ && !symminfo_->bb_is_independent( i ) ) continue; // only the main chain gets selected

		core::conformation::Residue const &rsd_i( pose.residue(i) );
		//if (rsd_i.name3()=="GLY") continue;

		Size asymm_num_i = ((i-1)%nres_)+1;
		Real i_dens_rscc = perrsd_dens_[asymm_num_i];
		Real sum    = i_dens_rscc;
		Real sq_sum = i_dens_rscc*i_dens_rscc;
		Size n_nbrs = 1;

		fragbias_tr << "rsd: " << i << " " << rsd_i.name3() << " i_dens_rscc: " << i_dens_rscc << std::endl;

		// get density score per i
		for ( graph::Graph::EdgeListConstIter
				iru  = energy_graph.get_node(i)->const_edge_list_begin(),
				irue = energy_graph.get_node(i)->const_edge_list_end();
				iru != irue; ++iru ) {

			EnergyEdge const * edge( static_cast< EnergyEdge const *> (*iru) );
			Size const j( edge->get_other_ind(i) );
			Size asymm_num_j = ((j-1)%nres_)+1;

			conformation::Residue const & rsd_j ( pose.residue(j) );
			//if (rsd_j.name3()=="GLY") continue;

			Real dist = std::pow( edge->square_distance(), 0.5 );
			if ( dist <= 10.0 ) {
				Real j_dens_rscc = perrsd_dens_[asymm_num_j];
				sum    += j_dens_rscc;
				sq_sum += j_dens_rscc*j_dens_rscc;
				n_nbrs ++;

				fragbias_tr << "computing energy graph: " << j
					//<< " dist: " <<  caled_dist
					<< " "                        << rsd_j.name3()
					<< " dist: "                 << dist
					<< " j_dens_rscc: "           << j_dens_rscc
					<< " sum: "                   << sum
					<< " sq_sum: "               << sq_sum
					<< " n_nbrs: "               << n_nbrs
					//<< " delta: "
					//<< dist-caled_dist
					<< std::endl;
			} else {
				continue;
			}
		} // res j

		fragbias_tr << " sum: " << sum
			<< " sq_sum: " << sq_sum
			<< " n_nrbs: " << n_nbrs
			<< std::endl;

		Real nbrdens_mean  = sum/n_nbrs;
		Real nbrdens_stdev = std::sqrt( sq_sum/n_nbrs - nbrdens_mean*nbrdens_mean );

		// z-score for rscc of residue i to rscc of its neighbors
		Real i_nbrdens_zscore = (i_dens_rscc - nbrdens_mean) / nbrdens_stdev;
		perrsd_nbrdens_[asymm_num_i] = i_nbrdens_zscore;

		// z-score for rscc of residue i to rscc of all residues
		Real i_perrsd_dens_zscore = (i_dens_rscc-perrsd_dens_mean)/perrsd_dens_stdev;

		fragbias_tr << "rsd: "              << asymm_num_i
			<< " rsn: "             << rsd_i.name3()
			<< " symm_rsd: "        << i
			<< " dens_rscc: "       << i_dens_rscc
			<< " nbrdens_mean: "    << nbrdens_mean
			<< " nbrdens_stdev: "   << nbrdens_stdev
			<< " i_nbr_zscore: "    << i_nbrdens_zscore
			<< " i_perrsd_zscore: " << i_perrsd_dens_zscore
			<< " nbrs: "            << n_nbrs
			<< std::endl;

		// assign nbr_zscore cutoff
		/*
		Real nbr_zscore_cutoff = -2;
		if ( n_nbrs > 20 ){   // the residue is buried
		nbr_zscore_cutoff = 0;
		} else if ( n_nbrs > 15 ){
		nbr_zscore_cutoff = -0.5;
		} else if ( n_nbrs > 10 ){
		nbr_zscore_cutoff = -1;
		}
		*/

		if ( i_perrsd_dens_zscore < -1.0 ) { // the rscc for this residue is being considered worse overall
			if ( i_nbrdens_zscore < -2.0 ) {  // the rscc for this residue is worse than its neighbors
				if ( cumulative_ ) {
					fragmentProbs_[asymm_num_i] += 1;
				} else {
					fragmentProbs_[asymm_num_i] = 1;
				}
			} // i_perrsd_zscore
		} // i_nbr_zscore

		fragbias_tr << "fragmentProbs_: " << asymm_num_i << " rsd: " << i << " prob: " << fragmentProbs_[asymm_num_i] << std::endl;
	} // for i in range(pose.total_residue())
}


void
FragmentBiasAssigner::
density(
	pose::Pose &pose
){
	fragbias_tr << "density is chose" << std::endl;
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

	perrsd_dens_.resize(nres_, 0.0);
	core::Real CCsum=0, CCsum2=0;

	// turn score_symm_complex flag on; otherwise it won't see symmetry, and will return 0
	for ( int r=1; r<=(int)nres_; ++r ) {
		int rsrc = r;
		if ( symminfo_ && !symminfo_->bb_is_independent(r) ) {
			rsrc = symminfo_->bb_follows(r);
		}

		perrsd_dens_[r] = core::scoring::electron_density::getDensityMap().matchRes( rsrc , pose.residue(rsrc), pose, symminfo_, false);
		CCsum += perrsd_dens_[r];
		CCsum2 += perrsd_dens_[r]*perrsd_dens_[r];
	}
	CCsum /= nres_;
	CCsum2 = sqrt( CCsum2/nres_-CCsum*CCsum );

	for ( int r=1; r<=(int)nres_; ++r ) {
		if ( perrsd_dens_[r]<0.6 ) {
			fragmentProbs_[r] = 1.0;
		} else if ( perrsd_dens_[r]<0.8 ) {
			fragmentProbs_[r] = 0.25;
		} else {
			fragmentProbs_[r] = 0.01;
		}
		//fragbias_tr.Debug << "residue " << r << ": " << " rscc=" << perrsd_dens_[r] << " weight=" <<fragmentProbs_[r] << std::endl;
		fragbias_tr << "residue " << r << " rscc: " << perrsd_dens_[r] << " weight: " <<fragmentProbs_[r] << std::endl;
	}
}


void
FragmentBiasAssigner::
geometry(
	pose::Pose &pose,
	Real weight  /*1.0*/
){
	fragbias_tr << "geometry is chose" << std::endl;
	fragProbs_assigned_=true;
	if ( score_threshold_ == 123456789 ) set_score_threshold( 0.6 );

	// clean the container
	perrsd_geometry_.resize(nres_, 0.0);
	cal_perrsd_score( pose,
		core::scoring::cart_bonded_angle,
		perrsd_geometry_,
		weight );

	assign_fragprobs( perrsd_geometry_,
		score_threshold_ );
}


void
FragmentBiasAssigner::
rama(
	pose::Pose &pose,
	Real weight /*=0.2*/
){
	fragProbs_assigned_=true;
	if ( score_threshold_ == 123456789 ) set_score_threshold( 0.7 );
	fragbias_tr << "rama is chosen, and the score_threshold is " << score_threshold_ << std::endl;

	// clean the container
	perrsd_rama_.resize(nres_, 0.0);
	cal_perrsd_score( pose,
		core::scoring::rama,
		perrsd_rama_,
		weight );

	assign_fragprobs( perrsd_rama_,
		score_threshold_ );
}


void
FragmentBiasAssigner::
old_rama(
	pose::Pose &pose
){
	using namespace core::scoring;
	fragProbs_assigned_=true;

	// find segments with worst rama score
	core::Real Rtemp=1;  // again, this is a guess

	// score the pose
	core::scoring::ScoreFunctionOP myscore( new core::scoring::ScoreFunction() );
	myscore->set_weight( core::scoring::rama, 1.0 );
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		myscore = core::scoring::symmetry::symmetrize_scorefunction(*myscore);
	}

	core::scoring::Energies & energies( pose.energies() );
	(*myscore)(pose);

	for ( int r=1; r<=(int)nres_; ++r ) {
		core::scoring::EnergyMap & emap( energies.onebody_energies( r ) );
		// i dont think this will work for symmetric systems where 1st subunit is not the scoring one
		core::Real ramaScore = emap[ core::scoring::rama ];
		fragmentProbs_[r] = exp( ramaScore / Rtemp );
		fragbias_tr << "Prob_dens_rama( " << r << " ) = " << fragmentProbs_[r] << " ; rama=" << ramaScore << std::endl;
	}
}


void
FragmentBiasAssigner::
bfactors(
	pose::Pose &pose
){
	fragbias_tr << "bfactors gets chosen" << std::endl;
	fragProbs_assigned_=true;

	// find segments with highest bfactors
	core::Real Btemp=25;  // no idea what value makes sense here
	// with Btemp = 25, a B=100 is ~54 times more likely to be sampled than B=0
	runtime_assert( pose.pdb_info() != 0 );
	for ( int r=1; r<=(int)nres_; ++r ) {
		core::Real Bsum=0;
		core::Size nbb = pose.residue_type(r).last_backbone_atom();
		for ( core::Size atm=1; atm<=nbb; ++atm ) {
			Bsum += pose.pdb_info()->temperature( r, atm );
		}
		Bsum /= nbb;
		fragmentProbs_[r] = exp( Bsum/Btemp );
		//fragbias_tr << "Prob_dens_bfact( " << r << " ) = " << fragmentProbs_[r] << " ; B=" << Bsum << std::endl;
	}
}


void
FragmentBiasAssigner::
chainbreak(
	pose::Pose &pose
){
	fragbias_tr << "chainbreak is chose" << std::endl;
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
	fragbias_tr << "report per-residue fragment bias" << std::endl;

	for ( int r=1; r<=(int)nres_; ++r ) {
		//if (!pose.residue_type(r).is_protein()) continue;
		//if (pose.fold_tree().is_cutpoint(r+1)) continue;

		std::cerr << "rsn: " << r << " prob: " << fragmentProbs_[r] << std::endl;
	}
}


} // hybridization
} // protocol
