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
/// @detailed
/// @author Frank DiMaio

#include <protocols/rbsegment_relax/OptimizeThreading.hh>
#include <protocols/rbsegment_relax/OptimizeThreadingCreator.hh>

#include <protocols/loops/loops_main.hh>
#include <protocols/simple_moves/FragmentMover.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Remarks.hh>

#include <protocols/moves/MonteCarlo.hh>

#include <core/fragment/picking_old/vall/util.hh>
#include <core/fragment/IndependentBBTorsionSRFD.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/util.hh>

#include <protocols/rbsegment_relax/util.hh>
#include <protocols/rbsegment_relax/RBSegmentMover.hh>
#include <protocols/rbsegment_relax/RBSegment.hh>
#include <protocols/hybridization/util.hh>
#include <protocols/loops/util.hh>
#include <protocols/moves/DsspMover.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/simple_moves/symmetry/SetupNCSMover.hh>

#include <protocols/moves/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <utility/tag/Tag.hh>
#include <numeric/random/random.hh>

#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>


using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace rbsegment_relax {

static basic::Tracer TR("protocols.rbsegment_relax.OptimizeThreading");

using namespace protocols;
using namespace core;


//////////////////
//////////////////
/// creator


std::string
OptimizeThreadingMoverCreator::keyname() const {
	return OptimizeThreadingMoverCreator::mover_name();
}

protocols::moves::MoverOP
OptimizeThreadingMoverCreator::create_mover() const {
	return new OptimizeThreadingMover;
}

std::string
OptimizeThreadingMoverCreator::mover_name() {
	return "OptimizeThreading";
}


//////////////////
//////////////////
/// mover

void OptimizeThreadingMover::apply( core::pose::Pose & pose ) {
	using namespace rbsegment_relax;

	core::Size nres = hybridization::get_num_residues_nonvirt( pose );
	while (!pose.residue_type(nres).is_protein()) --nres;
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		scorefxn_ = new core::scoring::symmetry::SymmetricScoreFunction( *scorefxn_ );
		scorefxn_sampling_ = new core::scoring::symmetry::SymmetricScoreFunction( *scorefxn_sampling_ );
	}

	// see if the pose has NCS
	simple_moves::symmetry::NCSResMappingOP ncs;
	if ( pose.data().has( core::pose::datacache::CacheableDataType::NCS_RESIDUE_MAPPING ) ) {
		ncs = ( static_cast< simple_moves::symmetry::NCSResMapping* >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::NCS_RESIDUE_MAPPING )() ));
	}

	// get sses & choose cutpoints
	core::scoring::dssp::Dssp dssp_obj( pose );
	dssp_obj.insert_ss_into_pose( pose );
	protocols::loops::Loops SSEs = protocols::loops::extract_secondary_structure_chunks( pose, "HE", 3, 6, 3, 4);
	SSEs.sequential_order();

	utility::vector1< RBResidueRange > segments;
	core::Size prevcut = 1;
	for (int j=2; j<=(int)SSEs.num_loop(); ++j) {
		core::Size j0stop = SSEs[j-1].stop();
		core::Size j1start = SSEs[j].start();

		// see if there is a cut between the two and use that instead
		bool found_cuts = false;
		for (int k=(int)j0stop; k<(int)j1start; ++k) {
			if (pose.fold_tree().is_cutpoint( k )) {
				segments.push_back( RBResidueRange( prevcut, k ) );
				//TR.Debug << "FC! Add segment: " << prevcut << " , " << k << std::endl;
				prevcut = k+1;
				found_cuts=true;
			}
		}

		if (!found_cuts) {
			core::Size cutpoint = (j0stop+j1start)/2;  // randomize?
			segments.push_back( RBResidueRange( prevcut, cutpoint-1 ) );
			//TR.Debug << "Add segment: " << prevcut << " , " << cutpoint-1 << std::endl;
			prevcut = cutpoint;
		}
	}
	segments.push_back( RBResidueRange( prevcut, nres ) );
	//TR.Debug << "F: Add segment: " << prevcut << " , " << nres << std::endl;

	core::Size nsegments=segments.size();
	for (int i=1; i<=(int)nsegments; ++i) {
		for (int j=i+1; j<=(int)nsegments; ++j) {
			bool crosses_cut = false;

			for (int k=i; k<j && !crosses_cut; ++k)
				crosses_cut = pose.fold_tree().is_cutpoint( segments[k].end() );

			if (!crosses_cut)
				segments.push_back( RBResidueRange( segments[i].start(), segments[j].end() ) );
		}
	}

	// mc loop
	SequenceShiftMover sshift( pose, RBSegment(), max_shift_ );

	core::Real best_score=1e30, acc_score=1e30;
	loops::LoopsOP best_loops = new loops::Loops();
	core::pose::Pose best_pose=pose, acc_pose;
	utility::vector1< int > offsets( nres, 0 );
	for (int i=1; i<=(int)nsteps_; ++i) {
		// random segment
		core::Size seg_i = numeric::random::random_range(1, segments.size() );

		if (i>1) {
			// apply movement & score pose
			utility::vector1 < RBResidueRange > ncs_segments ( 1, segments[seg_i] );
			if (ncs) {
				for (int j=1; j<=ncs->ngroups(); ++j ) {
					core::Size remap_start = ncs->get_equiv( j, segments[seg_i].start() );
					core::Size remap_stop = ncs->get_equiv( j, segments[seg_i].end() );
					for (int k=segments[seg_i].start(); k<=segments[seg_i].end() && remap_start==0; ++k)
						remap_start = ncs->get_equiv( j,k );
					if (remap_start==0) continue;  // undefined
					for (int k=segments[seg_i].end(); k>=segments[seg_i].start() && remap_stop==0; --k)
						remap_stop = ncs->get_equiv( j,k );
					//TR.Debug << "NCS: Add segment: " << remap_start << " , " << remap_stop << std::endl;
					ncs_segments.push_back( RBResidueRange(remap_start,remap_stop) );
				}
			}
			sshift.set_segment( RBSegment(ncs_segments) );
			sshift.apply( pose );
		}

		core::Real score_cst = (*scorefxn_)(pose);
		core::Real score_aln = sshift.score();

		// boltzmann
		core::Real score = score_cst + weight_*score_aln;
		core::Real boltz_factor = ( acc_score - score ) / temperature_;
		core::Real probability = std::exp( std::min (40.0, std::max(-40.0,boltz_factor)) );
		if ( probability < 1 ) {
			if ( numeric::random::uniform() >= probability ) {
				// reject
				pose = acc_pose;
			} else {
				// accept, not new best
				acc_pose = pose;
				acc_score = score;
				sshift.trigger_accept();
				TR << "Accept!    Step " << i << " score = " << score_cst << " + " << weight_ << " * " << score_aln << " = " << score << std::endl;
			}
		} else {
			// new best
			acc_pose = pose;
			acc_score = score;
			sshift.trigger_accept();
			if (score < best_score) {
				TR << "New best!  Step " << i << " score = " << score_cst << " + " << weight_ << " * " << score_aln << " = " << score << std::endl;
				best_pose = pose;
				best_loops = sshift.get_residues_to_rebuild();
				best_score = score;
			} else {
				TR << "Accept!    Step " << i << " score = " << score_cst << " + " << weight_ << " * " << score_aln << " = " << score << std::endl;
			}
		}
	}

	loops::LoopsOP loops = best_loops;
	TR << "Building:" << std::endl;
	TR << *loops << std::endl;

	pose = best_pose;
	if (rebuild_cycles_>0)
		rebuild_unaligned( pose, loops );
}

void
OptimizeThreadingMover::rebuild_unaligned(core::pose::Pose &pose, loops::LoopsOP loops) {
	if (loops->size() != 0) {
		// see if the pose has NCS
		simple_moves::symmetry::NCSResMappingOP ncs;
		if ( pose.data().has( core::pose::datacache::CacheableDataType::NCS_RESIDUE_MAPPING ) ) {
			ncs = ( static_cast< simple_moves::symmetry::NCSResMapping* >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::NCS_RESIDUE_MAPPING )() ));
		}

		core::kinematics::FoldTree f_in=pose.fold_tree(), f_new;

		// set foldtree + variants
		for( protocols::loops::Loops::iterator it=loops->v_begin(), it_end=loops->v_end(); it!=it_end; ++it ) {
			// if loop crosses a cut maintain that cut
			bool crosses_cut=false;
			Size i = 0;
			for( i=it->start(); i<=it->stop() && !crosses_cut; ++i ) {
				crosses_cut |= f_in.is_cutpoint(i);
			}
			if (crosses_cut) {
				it->set_cut(i-1); //?
			} else {
				it->set_cut( it->midpoint() );  // be deterministic so ncs copies match up
			}
			// to do? what if we cross two cuts?
		}


		protocols::loops::fold_tree_from_loops( pose, *loops, f_new);
		//std::cerr << f_new << std::endl;

		pose.fold_tree( f_new );
		protocols::loops::add_cutpoint_variants( pose );

		// set movemap
		core::kinematics::MoveMapOP mm_loop = new core::kinematics::MoveMap();
		for( protocols::loops::Loops::const_iterator it=loops->begin(), it_end=loops->end(); it!=it_end; ++it ) {
			for( Size i=it->start(); i<=it->stop(); ++i ) {
				mm_loop->set_bb(i, true);
				mm_loop->set_chi(i, true);
			}
		}

		// VERY FAST pick 3mers only in unaligned regions
		std::string tgt_seq = pose.sequence();
		core::fragment::FragSetOP frags3 = new core::fragment::ConstantLengthFragSet( 3 );
		for( protocols::loops::Loops::const_iterator it=loops->begin(), it_end=loops->end(); it!=it_end; ++it ) {
			for( Size i=it->start(); i+2<=it->stop(); ++i ) {
				core::fragment::FrameOP frame = new core::fragment::Frame( i, 3 );
				frame->add_fragment(
					core::fragment::picking_old::vall::pick_fragments_by_ss_plus_aa( "DDD", tgt_seq.substr( i-1, 3 ), 25, true, core::fragment::IndependentBBTorsionSRFD() ) );
				frags3->add( frame );
			}
		}

		//
		core::fragment::FragSetOP frags1 = new core::fragment::ConstantLengthFragSet( 1 );
		core::fragment::chop_fragments( *frags3, *frags1 );

		// make a vector of fragments for random access (why does FragSet not have this?!?)
		utility::vector1< core::fragment::FrameOP > frames1, frames3;
		for (core::fragment::FrameIterator it = frags1->nonconst_begin(); it != frags1->nonconst_end(); ++it) frames1.push_back( *it );
		for (core::fragment::FrameIterator it = frags3->nonconst_begin(); it != frags3->nonconst_end(); ++it) frames3.push_back( *it );

		// extend + idealize loops
		for( protocols::loops::Loops::const_iterator it=loops->begin(), it_end=loops->end(); it!=it_end; ++it ) {
			protocols::loops::Loop to_idealize( *it );
			protocols::loops::set_extended_torsions( pose, *it );
		}

		// setup MC
		core::Size nouterCyc=4, ninnerCyc=rebuild_cycles_;
		for (core::Size i=1; i<=nouterCyc; ++i) {
			if (i==nouterCyc-1) scorefxn_sampling_->set_weight( core::scoring::linear_chainbreak, 0.5 );
			if (i==nouterCyc)   scorefxn_sampling_->set_weight( core::scoring::linear_chainbreak, 2.0 );

			(*scorefxn_sampling_)(pose);
			protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( pose, *scorefxn_sampling_, 2.0 );

			for (Size n=1; n<=ninnerCyc; ++n) {
				utility::vector1< core::fragment::FrameOP > & working_frames = (n%2)?frames3:frames1;

				//frag3mover->apply( pose );
				int frame_idx = numeric::random::random_range(1,working_frames.size());
				core::fragment::FrameOP frame_i = working_frames[frame_idx];
				int frag_idx = numeric::random::random_range(1,frame_i->nr_frags());
				core::fragment::FragDataCOP to_insert = frame_i->fragment_ptr(frag_idx);
				to_insert->apply( pose, *frame_i );

				// now apply to ncs copies
				if (ncs) {
					for (int j=1; j<=ncs->ngroups(); ++j ) {
						bool all_are_mapped = true;
						for ( Size k=frame_i->start(); k<=frame_i->stop() && all_are_mapped; ++k )
							all_are_mapped &= (ncs->get_equiv( j,k )!=0);
						if (!all_are_mapped) continue;
						core::Size remap_start = ncs->get_equiv( j, frame_i->start() );
						core::Size remap_stop = ncs->get_equiv( j, frame_i->stop() );
						if (remap_stop-remap_start != frame_i->stop()-frame_i->start()) continue;
						to_insert->apply( pose, remap_start, remap_stop );
					}
				}

				(*scorefxn_sampling_)(pose);
				if (mc->boltzmann( pose , (n%2)?"frag3":"frag1" )) {
					;
					//std::ostringstream oss;
					//oss << "out" << i << "_" << n << ".pdb";
					//pose.dump_pdb( oss.str() );

					//std::cerr << "out " << i << "_" << n << ": " << frame_i->start() << std::endl;
					//for (int j=1; j<=ncs->ngroups(); ++j ) {
					//	bool all_are_mapped = true;
					//	for ( Size k=frame_i->start(); k<=frame_i->stop() && all_are_mapped; ++k ) all_are_mapped &= (ncs->get_equiv( j,k )!=0);
					//	if (!all_are_mapped) continue;
					//	core::Size remap_start = ncs->get_equiv( j, frame_i->start() );
					//	std::cerr << "  ---> " << remap_start << std::endl;
					//}
				}
			}
			mc->show_scores();
			mc->show_counters();
			mc->recover_low( pose );
		}

		// restore input ft
		protocols::loops::remove_cutpoint_variants( pose );
		pose.fold_tree( f_in );
	}
}


void OptimizeThreadingMover::parse_my_tag(
			utility::tag::TagPtr const tag,
			moves::DataMap & data,
			filters::Filters_map const & /*filters*/,
			moves::Movers_map const & /*movers*/,
			core::pose::Pose const & /*pose*/ )
{
	if( tag->hasOption( "scorefxn" ) ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "scorefxn" ) );
		scorefxn_ = (data.get< core::scoring::ScoreFunction * >( "scorefxns", scorefxn_name ))->clone();
	}
	if( tag->hasOption( "scorefxn_sampling" ) ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "scorefxn_sampling" ) );
		scorefxn_sampling_ = (data.get< core::scoring::ScoreFunction * >( "scorefxns", scorefxn_name ))->clone();
	}

	nsteps_ = tag->getOption<core::Size>( "nsteps", 1000 );
	rebuild_cycles_ = tag->getOption<core::Size>( "rebuild_cycles", 200 );
	weight_ = tag->getOption<core::Real>( "weight", 0.1 );
	temperature_ = tag->getOption<core::Real>( "temperature", 2.0 );
	max_shift_ = tag->getOption<core::Size>( "max_shift", 4 );
}

}
}
