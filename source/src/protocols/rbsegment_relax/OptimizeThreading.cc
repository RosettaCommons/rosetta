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

#include <protocols/rbsegment_relax/OptimizeThreading.hh>
#include <protocols/rbsegment_relax/OptimizeThreadingCreator.hh>

#include <protocols/loops/loops_main.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/PDBInfo.hh>

#include <protocols/moves/MonteCarlo.hh>

#include <core/fragment/picking_old/vall/util.hh>
#include <core/fragment/IndependentBBTorsionSRFD.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/import_pose/import_pose.hh>

#include <protocols/rbsegment_relax/RBSegmentMover.hh>
#include <protocols/rbsegment_relax/RBSegment.hh>
#include <protocols/hybridization/util.hh>
#include <protocols/loops/util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/symmetry/SetupNCSMover.hh>

#include <basic/datacache/DataMap.hh>

#include <utility/tag/Tag.hh>
#include <numeric/random/random.hh>

#include <basic/Tracer.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#include <core/chemical/ResidueType.hh> // AUTO IWYU For ResidueType
#include <core/kinematics/MoveMap.hh> // AUTO IWYU For MoveMap
#include <utility/string_util.hh> // AUTO IWYU For string_split
#include <protocols/loops/Loop.hh> // AUTO IWYU For Loop
#include <core/kinematics/FoldTree.hh> // AUTO IWYU For FoldTree, operator<<



namespace protocols {
namespace rbsegment_relax {

static basic::Tracer TR( "protocols.rbsegment_relax.OptimizeThreading" );

using namespace protocols;
using namespace core;


//////////////////
//////////////////
/// creator






//////////////////
//////////////////
/// mover

void OptimizeThreadingMover::apply( core::pose::Pose & pose ) {
	using namespace rbsegment_relax;

	core::Size nres = hybridization::get_num_residues_nonvirt( pose );
	while ( !pose.residue_type(nres).is_protein() ) --nres;

	// see if the pose has NCS
	symmetry::NCSResMappingOP ncs;
	if ( pose.data().has( core::pose::datacache::CacheableDataType::NCS_RESIDUE_MAPPING ) ) {
		ncs = ( utility::pointer::static_pointer_cast< symmetry::NCSResMapping > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::NCS_RESIDUE_MAPPING ) ));
	}

	// get sses & choose cutpoints
	protocols::loops::Loops SSEs;

	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		core::pose::Pose pose_asu;
		core::pose::symmetry::extract_asymmetric_unit(pose, pose_asu);
		core::scoring::dssp::Dssp dssp_obj( pose_asu );
		dssp_obj.insert_ss_into_pose( pose_asu );
		SSEs = protocols::loops::extract_secondary_structure_chunks( pose_asu, "HE", 3, 6, 3, 3, 4);
	} else {
		core::scoring::dssp::Dssp dssp_obj( pose );
		dssp_obj.insert_ss_into_pose( pose );
		SSEs = protocols::loops::extract_secondary_structure_chunks( pose, "HE", 3, 6, 3, 3, 4);
	}
	SSEs.sequential_order();


	// penalize insertions/deletions in the middle of SSEs
	utility::vector1< bool > sse_core(nres,false);
	for ( int j=1; j<=(int)SSEs.num_loop(); ++j ) {
		int jstart=SSEs[j].start(), jstop=SSEs[j].stop();
		for ( int k=jstart+4; k<=jstop-4; ++k ) {
			sse_core[ k ] = true;
		}
	}

	utility::vector1< RBResidueRange > segments;
	core::Size prevcut = 1;
	for ( int j=2; j<=(int)SSEs.num_loop(); ++j ) {
		core::Size j0stop = SSEs[j-1].stop();
		core::Size j1start = SSEs[j].start();

		std::string chain0 = pose.pdb_info()->chain(j0stop);
		std::string chain1 = pose.pdb_info()->chain(j1start);

		if ( chains_.size()>0 && std::find(chains_.begin(), chains_.end(), chain0) == chains_.end() ) continue;
		if ( chains_.size()>0 && std::find(chains_.begin(), chains_.end(), chain1) == chains_.end() ) continue;

		// see if there is a cut between the two and use that instead
		bool found_cuts = false;
		for ( auto k=(int)j0stop; k<(int)j1start; ++k ) {
			if ( pose.fold_tree().is_cutpoint( k ) ) {
				segments.push_back( RBResidueRange( prevcut, k ) );
				TR << "Add segment: " << prevcut << " , " << k << std::endl;
				prevcut = k+1;
				found_cuts=true;
			}
		}

		if ( !found_cuts ) {
			core::Size cutpoint = numeric::random::random_range(j0stop,j1start-1);
			segments.push_back( RBResidueRange( prevcut, cutpoint-1 ) );
			TR << "Add segment: " << prevcut << " , " << cutpoint-1 << std::endl;
			prevcut = cutpoint;
		}
	}
	segments.push_back( RBResidueRange( prevcut, nres ) );
	TR << "Add segment: " << prevcut << " , " << nres << std::endl;

	// combine segments
	core::Size nsegments=segments.size();
	for ( int i=1; i<=(int)nsegments; ++i ) {
		for ( int j=i+1; j<=(int)nsegments; ++j ) {
			bool crosses_cut = false;

			for ( int k=i; k<j && !crosses_cut; ++k ) {
				crosses_cut = pose.fold_tree().is_cutpoint( segments[k].end() );
			}

			if ( !crosses_cut ) {
				segments.push_back( RBResidueRange( segments[i].start(), segments[j].end() ) );
				TR << "Add segment: " << segments[i].start() << " , " << segments[j].end() << std::endl;
			}
		}
	}

	// mc loop
	SequenceShiftMover sshift( pose, RBSegment(), max_shift_ );
	sshift.set_extra_penalty( sse_core );

	core::Real best_score=1e30, acc_score=1e30;
	loops::LoopsOP best_loops( new loops::Loops() );
	core::pose::Pose best_pose=pose, acc_pose;
	utility::vector1< int > offsets( nres, 0 );

	if ( greedy_ ) {
		///
		/// greedy
		///
		core::Size steps_eff = nsteps_; //std::ceil( nsteps_ / (2.0*max_shift_*segments.size()));
		TR << "Greedy search for " << steps_eff << " steps" << std::endl;
		for ( int i=1; i<=(int)steps_eff; ++i ) {
			acc_pose = pose;
			core::Real score_cst = (*scorefxn_)(pose);
			core::Real score_aln = sshift.score();
			best_score = score_cst + weight_*score_aln;

			int bestSeg=1, bestK=0;
			for ( core::Size seg_i = 1; seg_i<=segments.size(); ++seg_i ) {
				utility::vector1 < RBResidueRange > ncs_segments ( 1, segments[seg_i] );
				if ( ncs ) {
					for ( core::uint j = 1; j <= ncs->ngroups(); ++j ) {
						core::Size remap_start = ncs->get_equiv( j, segments[seg_i].start() );
						core::Size remap_stop = ncs->get_equiv( j, segments[seg_i].end() );
						for ( core::uint k = segments[seg_i].start(); k <= segments[seg_i].end() && remap_start == 0; ++k ) {
							remap_start = ncs->get_equiv( j,k );
						}
						if ( remap_start==0 ) continue;  // undefined
						for ( core::uint k = segments[seg_i].end(); k >= segments[seg_i].start() && remap_stop == 0; --k ) {
							remap_stop = ncs->get_equiv( j,k );
						}
						//TR.Debug << "NCS: Add segment: " << remap_start << " , " << remap_stop << std::endl;
						ncs_segments.push_back( RBResidueRange(remap_start, remap_stop) );
					}
				}
				sshift.set_segment( RBSegment(ncs_segments) );

				int k_start = -1*(int)max_shift_;
				for ( int k=k_start; k<=(int)max_shift_; ++k ) {
					pose = acc_pose;
					if ( k!=0 ) {
						sshift.apply( pose, k );
						score_cst = (*scorefxn_)(pose);
						score_aln = sshift.score();

						// boltzmann
						core::Real score = score_cst + weight_*score_aln;
						if ( score < best_score ) {
							TR << "New best!  Step " << i << "." << seg_i << "[" << segments[seg_i].start() << "-" << segments[seg_i].end()
								<< "]." << k << " score = " << score_cst << " + " << weight_ << " * " << score_aln << " = " << score << std::endl;
							best_score = score;
							bestK=k;
							bestSeg=seg_i;
						}
					}
				}
			}

			// apply best shift we found
			{
				utility::vector1 < RBResidueRange > ncs_segments ( 1, segments[bestSeg] );
				if ( ncs ) {
					for ( core::uint j=1; j <= ncs->ngroups(); ++j ) {
						core::Size remap_start = ncs->get_equiv( j, segments[bestSeg].start() ),
							remap_stop = ncs->get_equiv( j, segments[bestSeg].end() );
						for ( core::uint k = segments[bestSeg].start(); k <= segments[bestSeg].end() && remap_start == 0; ++k ) {
							remap_start = ncs->get_equiv( j, k );
						}
						if ( remap_start==0 ) continue;
						for ( core::uint k = segments[bestSeg].end(); k >= segments[bestSeg].start() && remap_stop == 0; --k ) {
							remap_stop = ncs->get_equiv( j, k );
						}
						ncs_segments.push_back( RBResidueRange(remap_start,remap_stop) );
					}
				}
				if ( bestK==0 ) break;  // no move accepted last cycle


				pose = acc_pose;
				sshift.set_segment( RBSegment(ncs_segments) );
				sshift.apply( pose, bestK );
				sshift.trigger_accept();
			}
		}
		best_loops = sshift.get_residues_to_rebuild();
		*loops_ = *best_loops;
	} else {
		///
		/// non-greedy
		///
		for ( int i=1; i<=(int)nsteps_; ++i ) {
			// random segment
			core::Size seg_i = numeric::random::random_range(1, segments.size() );

			if ( i>1 ) {
				// apply movement & score pose
				utility::vector1 < RBResidueRange > ncs_segments ( 1, segments[seg_i] );
				if ( ncs ) {
					for ( core::uint j = 1; j <= ncs->ngroups(); ++j ) {
						core::Size remap_start = ncs->get_equiv( j, segments[seg_i].start() );
						core::Size remap_stop = ncs->get_equiv( j, segments[seg_i].end() );
						for ( core::uint k = segments[seg_i].start(); k <= segments[seg_i].end() && remap_start == 0; ++k ) {
							remap_start = ncs->get_equiv( j,k );
						}
						if ( remap_start==0 ) continue;  // undefined
						for ( core::uint k = segments[seg_i].end(); k >= segments[seg_i].start() && remap_stop == 0; --k ) {
							remap_stop = ncs->get_equiv( j,k );
						}
						//TR.Debug << "NCS: Add segment: " << remap_start << " , " << remap_stop << std::endl;
						ncs_segments.push_back( RBResidueRange(remap_start,remap_stop) );
					}
				}
				sshift.set_segment( RBSegment(ncs_segments) );
				sshift.apply( pose );
			}

			core::Real score_cst = (*scorefxn_)(pose);
			core::Real score_aln = sshift.score();

			if ( native_ ) {
				score_cst = core::scoring::CA_rmsd( pose, *native_ );
			}

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
					if ( !recover_low_ ) best_loops = sshift.get_residues_to_rebuild();
					TR << "Accept!    Step " << i << " score = " << score_cst << " + " << weight_ << " * " << score_aln << " = " << score << std::endl;
				}
			} else {
				// new best
				acc_pose = pose;
				acc_score = score;
				sshift.trigger_accept();
				if ( score < best_score ) {
					TR << "New best!  Step " << i << " score = " << score_cst << " + " << weight_ << " * " << score_aln << " = " << score << std::endl;
					best_pose = pose;
					best_loops = sshift.get_residues_to_rebuild();
					best_score = score;
				} else {
					TR << "Accept!    Step " << i << " score = " << score_cst << " + " << weight_ << " * " << score_aln << " = " << score << std::endl;
					if ( !recover_low_ ) best_loops = sshift.get_residues_to_rebuild();
				}
			}
			*loops_ = *best_loops;
		}
		if ( recover_low_ ) {
			pose = best_pose;
		}
	}

	TR << "Building:" << std::endl;
	TR << *loops_ << std::endl;

	rebuild_unaligned( pose );
}

void
OptimizeThreadingMover::rebuild_unaligned(core::pose::Pose &pose) {
	if ( loops_->size() != 0 && rebuild_cycles_!=0 ) {
		// see if the pose has NCS
		symmetry::NCSResMappingOP ncs;
		if ( pose.data().has( core::pose::datacache::CacheableDataType::NCS_RESIDUE_MAPPING ) ) {
			ncs = ( utility::pointer::static_pointer_cast< symmetry::NCSResMapping > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::NCS_RESIDUE_MAPPING ) ));
		}

		core::kinematics::FoldTree f_in=pose.fold_tree(), f_new;

		// set foldtree + variants
		for ( auto it=loops_->v_begin(), it_end=loops_->v_end(); it!=it_end; ++it ) {
			// if loop crosses a cut maintain that cut
			bool crosses_cut=false;
			core::Size i = 0;
			for ( i=it->start(); i<=it->stop() && !crosses_cut; ++i ) {
				crosses_cut |= f_in.is_cutpoint(i);
			}
			if ( crosses_cut ) {
				it->set_cut(i-1); //?
			} else {
				it->set_cut( it->midpoint() );  // be deterministic so ncs copies match up
			}
			// to do? what if we cross two cuts?
		}

		protocols::loops::fold_tree_from_loops( pose, *loops_, f_new);
		std::cerr << f_new << std::endl;

		pose.fold_tree( f_new );
		protocols::loops::add_cutpoint_variants( pose );

		// set movemap
		core::kinematics::MoveMapOP mm_loop( new core::kinematics::MoveMap() );
		for ( auto const & it : *loops_ ) {
			for ( core::Size i=it.start(); i<=it.stop(); ++i ) {
				mm_loop->set_bb(i, true);
				mm_loop->set_chi(i, true);
			}
		}

		// extend + idealize loops
		protocols::loops::set_extended_torsions_and_idealize_loops( pose, *loops_ );

		// pick 3mers only in unaligned regions
		std::string tgt_seq = pose.sequence();
		core::fragment::FragSetOP frags3( new core::fragment::ConstantLengthFragSet( 3 ) );
		for ( auto const & it : *loops_ ) {
			for ( core::Size i=it.start(); i+2<=it.stop(); ++i ) {
				core::fragment::FrameOP frame( new core::fragment::Frame( i, 3 ) );
				frame->add_fragment(
					core::fragment::picking_old::vall::pick_fragments_by_ss_plus_aa( "DDD", tgt_seq.substr( i-1, 3 ), 25, true, core::fragment::IndependentBBTorsionSRFD() ) );
				frags3->add( frame );
			}
		}

		core::fragment::FragSetOP frags1( new core::fragment::ConstantLengthFragSet( 1 ) );
		core::fragment::chop_fragments( *frags3, *frags1 );

		// make a vector of fragments for random access (why does FragSet not have this?!?)
		utility::vector1< core::fragment::FrameOP > frames1, frames3;
		for ( core::fragment::FrameIterator it = frags1->nonconst_begin(); it != frags1->nonconst_end(); ++it ) frames1.push_back( *it );
		for ( core::fragment::FrameIterator it = frags3->nonconst_begin(); it != frags3->nonconst_end(); ++it ) frames3.push_back( *it );

		// setup MC
		core::Size nouterCyc=4, ninnerCyc=rebuild_cycles_;
		for ( core::Size i=1; i<=nouterCyc; ++i ) {
			if ( i==nouterCyc-1 ) scorefxn_sampling_->set_weight( core::scoring::linear_chainbreak, 0.5 );
			if ( i==nouterCyc )   scorefxn_sampling_->set_weight( core::scoring::linear_chainbreak, 2.0 );

			(*scorefxn_sampling_)(pose);
			protocols::moves::MonteCarloOP mc( new protocols::moves::MonteCarlo( pose, *scorefxn_sampling_, 2.0 ) );

			for ( core::Size n=1; n<=ninnerCyc; ++n ) {
				utility::vector1< core::fragment::FrameOP > & working_frames = (n%2)?frames3:frames1;

				//frag3mover->apply( pose );
				int frame_idx = numeric::random::random_range(1,working_frames.size());
				core::fragment::FrameOP frame_i = working_frames[frame_idx];
				int frag_idx = numeric::random::random_range(1,frame_i->nr_frags());
				core::fragment::FragDataCOP to_insert = frame_i->fragment_ptr(frag_idx);
				to_insert->apply( pose, *frame_i );

				// now apply to ncs copies
				if ( ncs ) {
					for ( core::uint j = 1; j <= ncs->ngroups(); ++j ) {
						bool all_are_mapped = true;
						for ( core::Size k=frame_i->start(); k<=frame_i->stop() && all_are_mapped; ++k ) {
							all_are_mapped &= (ncs->get_equiv( j,k )!=0 && loops_->is_loop_residue(ncs->get_equiv( j,k )) );
						}
						if ( !all_are_mapped ) continue;
						core::Size remap_start = ncs->get_equiv( j, frame_i->start() );
						core::Size remap_stop = ncs->get_equiv( j, frame_i->stop() );
						if ( remap_stop-remap_start != frame_i->stop()-frame_i->start() ) continue;
						to_insert->apply( pose, remap_start, remap_stop );
					}
				}

				(*scorefxn_sampling_)(pose);
				if ( mc->boltzmann( pose , (n%2)?"frag3":"frag1" ) ) {
					;

					//std::ostringstream oss;
					//oss << "out" << i << "_" << n << ".pdb";
					//pose.dump_pdb( oss.str() );

					//std::cerr << "out " << i << "_" << n << ": " << frame_i->start() << std::endl;
					//for (int j=1; j<=ncs->ngroups(); ++j ) {
					// bool all_are_mapped = true;
					// for ( core::Size k=frame_i->start(); k<=frame_i->stop() && all_are_mapped; ++k ) all_are_mapped &= (ncs->get_equiv( j,k )!=0);
					// if (!all_are_mapped) continue;
					// core::Size remap_start = ncs->get_equiv( j, frame_i->start() );
					// std::cerr << "  ---> " << remap_start << std::endl;
					//}
				}
			}
			mc->show_scores();
			mc->show_counters();
			if ( recover_low_ ) mc->recover_low( pose );
		}

		// restore input ft
		protocols::loops::remove_cutpoint_variants( pose );
		pose.fold_tree( f_in );
	}
}


void OptimizeThreadingMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
)
{
	if ( tag->hasOption( "scorefxn" ) ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "scorefxn" ) );
		scorefxn_ = (data.get< core::scoring::ScoreFunction * >( "scorefxns", scorefxn_name ))->clone();
	}
	if ( tag->hasOption( "scorefxn_sampling" ) ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "scorefxn_sampling" ) );
		scorefxn_sampling_ = (data.get< core::scoring::ScoreFunction * >( "scorefxns", scorefxn_name ))->clone();
	}

	if ( tag->hasOption( "native" ) ) {
		std::string ref_model_pdb = tag->getOption<std::string>( "native" );
		native_ = utility::pointer::make_shared< core::pose::Pose >();
		core::import_pose::pose_from_file( *native_, ref_model_pdb , core::import_pose::PDB_file);
	}

	if ( tag->hasOption( "chains" ) ) {
		std::string chain_str = tag->getOption<std::string>( "chains" );
		utility::vector1< std::string > chns = utility::string_split(chain_str, ',');
		for ( int i=1; i<=(int)chns.size(); ++i ) chains_.push_back( chns[i] );
	}

	nsteps_ = tag->getOption<core::Size>( "nsteps", 5000 );
	step_penalty_ = tag->getOption<bool>( "step_penalty", false );
	recover_low_ = tag->getOption<bool>( "recover_low", true );
	greedy_ = tag->getOption<bool>( "greedy", false );
	rebuild_cycles_ = tag->getOption<core::Size>( "rebuild_cycles", 200 );
	weight_ = tag->getOption<core::Real>( "weight", 0.1 );
	temperature_ = tag->getOption<core::Real>( "temperature", 2.0 );
	max_shift_ = tag->getOption<core::Size>( "max_shift", 4 );


	// different nstep default for greedy
	if ( greedy_ && !tag->hasOption("nsteps") ) nsteps_=2;

	// add loopsOP to the DataMap
	if ( tag->hasOption( "loops_out" ) ) {
		std::string looptag = tag->getOption<std::string>( "loops_out" );
		data.add( "loops", looptag, loops_ );
	}
}

std::string OptimizeThreadingMover::get_name() const {
	return mover_name();
}

std::string OptimizeThreadingMover::mover_name() {
	return "OptimizeThreading";
}

void OptimizeThreadingMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "scorefxn", xs_string, "Scorefunction to use for threading" )
		+ XMLSchemaAttribute( "scorefxn_sampling", xs_string, "Scorefunction to use for Monte Carlo sampling" )
		+ XMLSchemaAttribute( "native", xs_string, "Path to PDB file containing native pose" )
		+ XMLSchemaAttribute( "chains", xsct_chain_cslist, "Specify chains within the pose to use for threading" )
		+ XMLSchemaAttribute::attribute_w_default( "nsteps", xsct_non_negative_integer, "Number of monte carlo steps", "5000" )
		+ XMLSchemaAttribute::attribute_w_default( "step_penalty", xsct_rosetta_bool, "This attribute is never actually used", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "recover_low", xsct_rosetta_bool, "Recover the lowest energy structure", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "greedy", xsct_rosetta_bool, "Perform a greedy alignment", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "rebuild_cycles", xsct_non_negative_integer, "Number of loop modeling cycles", "200" )
		+ XMLSchemaAttribute::attribute_w_default( "weight", xsct_real, "Weight for the sequence shift mover's score", "0.1" )
		+ XMLSchemaAttribute::attribute_w_default( "temperature", xsct_real, "Temperature for the Metropolis criterion", "2.0" )
		+ XMLSchemaAttribute::attribute_w_default( "max_shift", xsct_non_negative_integer, "Maximum number of residues to shift sequence", "4" )
		+ XMLSchemaAttribute( "loops_out", xs_string, "Add the LoopsOP to the datamap with this name" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Optimize threading of a sequence onto a pose", attlist );
}

std::string OptimizeThreadingMoverCreator::keyname() const {
	return OptimizeThreadingMover::mover_name();
}

protocols::moves::MoverOP
OptimizeThreadingMoverCreator::create_mover() const {
	return utility::pointer::make_shared< OptimizeThreadingMover >();
}

void OptimizeThreadingMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	OptimizeThreadingMover::provide_xml_schema( xsd );
}


}
}
