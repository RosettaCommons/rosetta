// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LoophashAssemblyMover.cc
///
/// @brief Derived from AssemblyMover, this mover generates an assembly and then uses Loophash segments
/// to connect any broken segments resulting from the assembly process
/// @author Tim Jacobs

// Unit Headers
#include <devel/sewing/sampling/LoophashAssemblyMover.hh>
#include <devel/sewing/sampling/LoophashAssemblyMoverCreator.hh>

//Protocol headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/Energies.hh>

#include <core/pose/util.hh>

#include <protocols/loops/Loops.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/util.hh>

#include <protocols/analysis/LoopAnalyzerMover.hh>

#include <protocols/simple_moves/MinMover.hh>

#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LoopHashMap.hh>
#include <protocols/loophash/LocalInserter.hh>
#include <protocols/sic_dock/loophash_util.hh>
#include <protocols/relax/AtomCoordinateCstMover.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

//Utility headers
#include <basic/Tracer.hh>

#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/keys/sewing.OptionKeys.gen.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

#include <devel/denovo_design/ConnectJumps.hh>

namespace devel {
namespace sewing {

static basic::Tracer TR( "devel.sewing.LoophashAssemblyMover" );

////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Boiler Plate Code   ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
protocols::moves::MoverOP
LoophashAssemblyMoverCreator::create_mover() const
{
	return new LoophashAssemblyMover;
}

std::string
LoophashAssemblyMoverCreator::keyname() const
{
	return LoophashAssemblyMoverCreator::mover_name();
}

std::string
LoophashAssemblyMoverCreator::mover_name()
{
	return "LoophashAssemblyMover";
}

protocols::moves::MoverOP
LoophashAssemblyMover::clone() const {
	return( protocols::moves::MoverOP( new LoophashAssemblyMover( *this ) ) );
}
protocols::moves::MoverOP
LoophashAssemblyMover::fresh_instance() const {
	return protocols::moves::MoverOP( new LoophashAssemblyMover );
}

std::string
LoophashAssemblyMover::get_name() const {
	return "LoophashAssemblyMover";
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  LoophashAssemblyMover functions   //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

LoophashAssemblyMover::LoophashAssemblyMover():
	AssemblyMover(),
	lh_library_(0),
	max_loop_segments_to_build_(100),
	max_loop_distance_(14.0),
	min_flanking_residues_(2),
	max_flanking_residues_(4),
	check_flanking_rms_(true),
	flanking_rms_cutoff_(20.0),
	max_insertion_rms_(0.3),
	loop_scorefxn_(0),
	loop_refine_scorefxn_(0)
{
	init();
}

void
LoophashAssemblyMover::init(){
	using namespace core;
	using namespace core::scoring;
	using namespace protocols::loophash;
	using namespace basic::resource_manager;
	using namespace basic::options;

	if ( ResourceManager::get_instance()->has_resource_with_description("LoopHashLibrary") ) {
		TR << "Retrieving lh library from resource manager." << std::endl;
		lh_library_ = get_resource<LoopHashLibrary>( "LoopHashLibrary" );
	}
	else{
		utility::vector1<core::Size> loop_sizes = option[ OptionKeys::lh::loopsizes ].value();
		TR << "Initializing lh library from command line" << std::endl;
		lh_library_ = new LoopHashLibrary( loop_sizes );
		lh_library_->load_mergeddb();
	}

	loop_scorefxn_ = new ScoreFunction();
	loop_scorefxn_->set_weight(coordinate_constraint,  0.5 );
	loop_scorefxn_->set_weight(cart_bonded,  0.5 );

	loop_refine_scorefxn_ = new ScoreFunction();
	//loop_refine_scorefxn_->set_weight(coordinate_constraint,  0.5 );
	loop_refine_scorefxn_->set_weight( cart_bonded,  0.5 );
	loop_refine_scorefxn_->set_weight( env			, 1.0);
	loop_refine_scorefxn_->set_weight( pair		 , 1.0);
	loop_refine_scorefxn_->set_weight( cbeta		, 1.0);
	loop_refine_scorefxn_->set_weight( vdw			, 1.0);
	loop_refine_scorefxn_->set_weight( rg			 , 3.0);
	loop_refine_scorefxn_->set_weight( cenpack	, 1.0);
	loop_refine_scorefxn_->set_weight( hs_pair	, 1.0);
	loop_refine_scorefxn_->set_weight( ss_pair	, 1.0);
	loop_refine_scorefxn_->set_weight( rsigma	 , 1.0);
	loop_refine_scorefxn_->set_weight( sheet		, 1.0);
}



bool
LoophashAssemblyMover::complete_assembly(
	AssemblyOP & assembly
){

	//First filter using base-class filters
	bool complete = AssemblyMover::complete_assembly(assembly);
	if(!complete) {
		return false;
	}

	//Now check for connectability and add loops
	if( ! basic::options::option[ basic::options::OptionKeys::sewing::skip_loop_generation ].value() ) {
		bool succesfully_rearranged = rearrange_assembly(assembly);
		if(!succesfully_rearranged) { return false; }

		core::pose::Pose cen_pose = assembly->to_pose(core::chemical::CENTROID);

		//devel::denovo_design::ConnectJumpsOP connect_jumps = new devel::denovo_design::ConnectJumps;
		//connect_jumps->jump1(1);
		//connect_jumps->jump2(2);
		//connect_jumps->apply(cen_pose);
		//if(basic::options::option[basic::options::OptionKeys::sewing::dump_pdbs]) {
		//	cen_pose.dump_pdb("reordered_assembly.pdb");
		//}

		built_loops_ = add_loophash_segments(assembly, cen_pose);
		if(built_loops_.empty()) {
			TR << "Failed to find loops for Assembly" << std::endl;
			return false;
		}
	}
	return true;
}




void
LoophashAssemblyMover::output_stats(
	AssemblyOP const & assembly,
	core::pose::Pose & pose
) {

	//Report all the base-class stats too
	AssemblyMover::output_stats(assembly, pose);

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	protocols::jd2::JobOP const job_me ( protocols::jd2::JobDistributor::get_instance()->current_job() );
	std::string const job_name ( protocols::jd2::JobDistributor::get_instance()->job_outputter()->output_name(job_me) );

	job_me->add_string_real_pair( "n_loops", built_loops_.num_loop());

	protocols::analysis::LoopAnalyzerMoverOP loop_analyzer = new protocols::analysis::LoopAnalyzerMover(built_loops_, true);
	//loop_analyzer->apply(pose);
	core::Real total_loops_score = loop_analyzer->get_total_score();
	job_me->add_string_real_pair( "loops_score", total_loops_score);
}



bool
LoophashAssemblyMover::rearrange_assembly(
	AssemblyOP & assembly
) const {

	utility::vector1< utility::vector1<core::Size> > orders = assembly->find_possible_orders(max_loop_distance_);
	numeric::random::random_permutation(orders, numeric::random::rg());
	TR << "Permuting " << orders.size() << " valid assembly orders" << std::endl;

	AssemblyOP reordered_assembly = assembly->clone();
	for(core::Size i=1; i<=orders.size(); ++i) {
		reordered_assembly = assembly->clone();
		reordered_assembly->reorder(orders[i]);
		core::pose::Pose assembly_pose = reordered_assembly->to_pose(core::chemical::CENTROID);
		//core::Size num_lh_frags = count_loophash_fragments(reordered_assembly, assembly_pose);
		//if((int)num_lh_frags >= basic::options::option[ basic::options::OptionKeys::sewing::min_lh_fragments].value()) {
			assembly = reordered_assembly;
			return true;
		//}
	}
	TR << "Assembly could not be rearranged in a connectable way." << std::endl;
	return false;
}



///@details count the number of loophash fragments for each
///unclosed loop in the current Assembly. Return the number
///of fragments for the jump with the least number of fragments
core::Size
LoophashAssemblyMover::count_loophash_fragments(
	AssemblyOP const assembly,
	core::pose::Pose const & pose
) const {
	using namespace protocols::loophash;
	using namespace basic::options;

	utility::vector1<core::Size> loop_sizes = option[ OptionKeys::lh::loopsizes ].value();
	core::Size max_radius = option[ OptionKeys::lh::max_radius ].value();

	utility::vector1<core::Size> loop_anchors = assembly->pose_loop_anchors();
	TR << "Checking lh fragments for " << loop_anchors.size() << " jumps." << std::endl;

	core::Size min_total = 100000000;
	for(core::Size i=1; i<=loop_anchors.size(); ++i ) {
		core::Size anchor_total = 0;

		for(core::Size n_overlapping_n = min_flanking_residues_; n_overlapping_n <= max_flanking_residues_; ++n_overlapping_n) {
			for(core::Size n_overlapping_c = min_flanking_residues_; n_overlapping_c <= max_flanking_residues_; ++n_overlapping_c) {

				numeric::geometry::hashing::Real6 loop_transform;
				if(!get_rt_over_leap_without_foldtree_bs( pose, loop_anchors[i]-n_overlapping_n, loop_anchors[i]+1+n_overlapping_c, loop_transform )){
					utility_exit_with_message("Unable to find rigid body transform over jump");
				}

				core::Real dist_sq = pose.residue(loop_anchors[i]).atom("CA").xyz().distance_squared(pose.residue(loop_anchors[i]+1).atom("CA").xyz());
				if(TR.Debug.visible()) {
					TR.Debug << "Loop anchor " << loop_anchors[i] << " distance sq " << dist_sq << std::endl;
				}

				for(core::Size j=1; j<=loop_sizes.size(); ++j) {
					//If this distance is too far, don't evaluate linker
					if( dist_sq >= (10.0*core::Real(loop_sizes[j])*core::Real(loop_sizes[j])) ) continue;
					LoopHashMap hashmap = lh_library_->gethash( loop_sizes[j] );
					anchor_total += hashmap.radial_count(max_radius, loop_transform);
				}
			}
		}
		TR << "Loop anchor " << loop_anchors[i] << " has " << anchor_total << " LH fragments" << std::endl;
		if(anchor_total == 0){
			return 0;
		}
		min_total = std::min(min_total, anchor_total);
	}
	return min_total;
}




protocols::loops::Loops
LoophashAssemblyMover::add_loophash_segments(
	AssemblyOP & assembly,
	core::pose::Pose & pose
) const{

	using namespace protocols::loophash;

	protocols::loops::Loops loophash_loops;

	utility::vector1<core::Size> loop_anchors = assembly->pose_loop_anchors();
	utility::vector1<core::Size> disconnected_segments = assembly->disconnected_segments();
	runtime_assert(loop_anchors.size() == disconnected_segments.size());

	for(core::Size i=1; i<=loop_anchors.size(); ++i) {
		TR.Debug << "building loop after residue " << loop_anchors[i] << std::endl;


		core::Size n_segment_start = assembly->pose_num(
			assembly->segments()[disconnected_segments[i]].model_id_, assembly->segments()[disconnected_segments[i]].residues_[1].resnum_);
		core::Size c_segment_end = assembly->pose_num(
			assembly->segments()[disconnected_segments[i]+1].model_id_, assembly->segments()[disconnected_segments[i]+1].residues_.back().resnum_);
		protocols::loops::Loop new_loop = add_single_loop(pose, loop_anchors[i], n_segment_start, c_segment_end);

		//If we failed to create this loop, then return an empty loops set
		if(new_loop.start() == 0 || new_loop.length() == 0) {
			loophash_loops.clear();
			return loophash_loops;
		}

		if(basic::options::option[basic::options::OptionKeys::sewing::dump_pdbs]) {
			pose.dump_pdb("best_loop_" + utility::to_string(i) + ".pdb");
		}

		assembly->add_loop_segment(pose, new_loop, disconnected_segments[i]);
		assembly->update_coords_from_pose(pose);
		for(core::Size j=i+1; j<=loop_anchors.size(); ++j) {
			loop_anchors[j]+=new_loop.length();
			++disconnected_segments[j];
		}
		loophash_loops.add_loop(new_loop);
	}

	return loophash_loops;
}



protocols::loops::Loop
LoophashAssemblyMover::add_single_loop(
	core::pose::Pose & pose,
	core::Size loop_anchor,
	core::Size n_segment_start,
	core::Size c_segment_end
) const {

	core::Real best_score = 1000000000.0;
	protocols::loops::Loop best_loop;
	core::pose::Pose best_pose;

	core::Size attempt = 0;
	for(core::Size n_overlapping_n = max_flanking_residues_; n_overlapping_n >= min_flanking_residues_; --n_overlapping_n) {
		for(core::Size n_overlapping_c = max_flanking_residues_; n_overlapping_c >= min_flanking_residues_; --n_overlapping_c) {

			if(TR.Debug.visible()) {
				TR.Debug << "Number of overlappign residues N-terminal of loop: " << n_overlapping_n << std::endl;
				TR.Debug << "Number of overlappign residues C-terminal of loop: " << n_overlapping_c << std::endl;
			}

			int loophash_fragment_start = (int)loop_anchor - (int)n_overlapping_n;
			int loophash_fragment_end = (int)loop_anchor + (int)n_overlapping_c;
			if(loophash_fragment_start < (int)n_segment_start || loophash_fragment_end > (int)c_segment_end){ continue; }

			utility::vector1< std::pair< protocols::loophash::BackboneSegment, std::string > > bb_segs =
				get_backbone_segments(pose, loophash_fragment_start, loophash_fragment_end);
			numeric::random::random_permutation(bb_segs, numeric::random::rg());

			if(check_flanking_rms_) {
				trim_bb_segs(pose, loop_anchor, n_overlapping_n, n_overlapping_c, bb_segs);
			}

			core::scoring::ScoreFunctionOP loop_scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "score4_smooth_cart" );
			loop_scorefxn->set_weight(core::scoring::cart_bonded, 1.0);

			for(core::Size i=1; i<=bb_segs.size(); ++i) {
				++attempt;
				if(attempt > max_loop_segments_to_build_) break;

				core::Size loop_size = bb_segs[i].first.length() - n_overlapping_n - n_overlapping_c - 1;

				std::pair<core::pose::Pose, core::Real> build_results =
					build_loop_pose(pose, loop_anchor, n_segment_start, c_segment_end, n_overlapping_n, n_overlapping_c, bb_segs[i].first, bb_segs[i].second);
				pose.constraint_set(NULL);

				core::pose::Pose loop_pose = build_results.first;
				core::Real insertion_rms = build_results.second;

				if(basic::options::option[basic::options::OptionKeys::sewing::dump_pdbs]) {
					loop_pose.dump_pdb("loop_pose_" + utility::to_string(attempt) + ".pdb");
				}

				loop_pose.energies().clear();
				core::Real loop_score = loop_refine_scorefxn_->score(loop_pose);

				protocols::loops::Loop loop = protocols::loops::Loop(loop_anchor+1, loop_anchor+loop_size);

				core::kinematics::FoldTree simple_ft = core::kinematics::FoldTree((int)loop_pose.total_residue());
				loop_pose.fold_tree(simple_ft);

				if(TR.Debug.visible()) {
					TR.Debug << "Loop score for loop (" << loop << ") build: " << attempt << ": " << loop_score << " RMS: " << insertion_rms << std::endl;
				}
				if(insertion_rms < max_insertion_rms_ && loop_score < best_score) {
					best_score = loop_score;
					best_loop = loop;
					best_pose = loop_pose;
				}
			}
		}
	}
	pose=best_pose;
	TR << "Best loop score for anchor: " << loop_anchor << " " << best_score << std::endl;

	return best_loop;
}


void
LoophashAssemblyMover::trim_bb_segs(
	core::pose::Pose const & pose,
	core::Size loop_anchor,
	core::Size n_overlapping_n,
	core::Size n_overlapping_c,
	BackboneSegments & bb_segs
) const {

	core::Size loophash_fragment_start = loop_anchor - n_overlapping_n;
	core::Size loophash_fragment_end = loop_anchor + n_overlapping_c;

	BackboneSegments::iterator it = bb_segs.begin();
	while(it != bb_segs.end()) {
		protocols::loophash::BackboneSegment backbone_seg = it->first;

		int loop_size = backbone_seg.length() - n_overlapping_n - n_overlapping_c - 1;
		if(loop_size < 1) {
			it = bb_segs.erase(it);
			continue;
		}

		//Check the rms
		std::vector<core::Real> phi = backbone_seg.phi();
		std::vector<core::Real> psi = backbone_seg.psi();
		std::vector<core::Real> omega = backbone_seg.omega();

		core::Real sumsqr = 0;
		core::Real diff = 0;
		core::Size count = 0;
		for(core::Size i=0; i<n_overlapping_n+1; ++i) {
			core::Real pose_phi = pose.phi(loophash_fragment_start+i);
			core::Real pose_psi = pose.psi(loophash_fragment_start+i);
			core::Real pose_omega = pose.omega(loophash_fragment_start+i);

			if(TR.Debug.visible()) {
				TR.Debug << loophash_fragment_start+i << ": " << pose_phi << " " << pose_psi << " " << pose_omega << std::endl;
			}

			diff = phi[i] - pose_phi;
			while( diff > 180  ) diff -= 360;
			while( diff < -180  ) diff += 360;
			sumsqr += diff*diff;
			++count;

			if(i < n_overlapping_n) {
				diff = psi[i] - pose_psi;
				while( diff > 180  ) diff -= 360;
				while( diff < -180  ) diff += 360;
				sumsqr += diff*diff;
				++count;

				diff = omega[i] - pose_omega;
				while( diff > 180  ) diff -= 360;
				while( diff < -180  ) diff += 360;
				sumsqr += diff*diff;
				++count;
			}
		}

		for(core::Size i=0; i<n_overlapping_c; ++i) {
			core::Real pose_phi = pose.phi(loophash_fragment_end-i);
			core::Real pose_psi = pose.psi(loophash_fragment_end-i);
			core::Real pose_omega = pose.omega(loophash_fragment_end-i);

			if(TR.Debug.visible()) {
				TR.Debug << loophash_fragment_end-i << ": " << pose_phi << " " << pose_psi << " " << pose_omega << std::endl;
			}

			if(i < n_overlapping_c - 1) {
				diff = phi[phi.size()-1-i] - pose_phi;
				while( diff > 180  ) diff -= 360;
				while( diff < -180  ) diff += 360;
				sumsqr += diff*diff;
				++count;
			}

			diff = psi[psi.size()-1-i] - pose_psi;
			while( diff > 180  ) diff -= 360;
			while( diff < -180  ) diff += 360;
			sumsqr += diff*diff;
			++count;

			diff = omega[omega.size()-1-i] - pose_omega;
			while( diff > 180  ) diff -= 360;
			while( diff < -180  ) diff += 360;
			sumsqr += diff*diff;
			++count;
		}

		core::Real rms = std::sqrt(sumsqr/(core::Real)count);
		if(TR.Debug.visible()) {
			TR.Debug << "RMS: " << rms << std::endl;
		}

		if(rms > flanking_rms_cutoff_) {
			it = bb_segs.erase(it);
		}
		else {
			++it;
		}
	}

}

///@details query the loophash library for backbone segments across the jump
///from loop_anchor to loop_anchor+1. Return a map from loop size to backbone segment
LoophashAssemblyMover::BackboneSegments
LoophashAssemblyMover::get_backbone_segments(
	core::pose::Pose & pose,
	core::Size loophash_fragment_start,
	core::Size loophash_fragment_end
) const {
	using namespace protocols::loophash;
	using namespace basic::options;

	utility::vector1<core::Size> loop_sizes = option[ OptionKeys::lh::loopsizes ].value();

	BBData bb_data;
	BBExtraData extra_data;

	//Get the loop transform
	numeric::geometry::hashing::Real6 loop_transform;
	if(!get_rt_over_leap_without_foldtree_bs( pose, loophash_fragment_start, loophash_fragment_end+1, loop_transform )){
		utility_exit_with_message("Unable to find rigid body transform over jump");
	}

	//Save all backbone segments from this transform
	BackboneSegments bb_segs;
	for(core::Size size_ind = 1; size_ind <= loop_sizes.size(); ++size_ind) {
		core::Size fragment_size = loop_sizes[size_ind];
		LoopHashMap hashmap = lh_library_->gethash( fragment_size );

		std::vector<core::Size> leap_index_bucket;
		hashmap.radial_lookup( core::Size(4/*max radius*/), loop_transform, leap_index_bucket);

		TR << "Leap index bucket size (" << fragment_size << " residues): " << leap_index_bucket.size() << std::endl;
		for(core::Size i=0; i<leap_index_bucket.size(); ++i) {
			BackboneSegment backbone_seg;
			LeapIndex cp = hashmap.get_peptide( leap_index_bucket[i] );
			lh_library_->backbone_database().get_backbone_segment( cp.index, cp.offset , fragment_size , backbone_seg );
			lh_library_->backbone_database().get_protein( cp.index, bb_data );
			lh_library_->backbone_database().get_extra_data(bb_data.extra_key, extra_data);

			runtime_assert(cp.offset % 3 == 0);
			core::Size seq_offset = (cp.offset/3);
			++seq_offset;//The first residue of the loop fragment is the same position as the anchor
			std::string sequence = extra_data.sequence;
			std::string loop_sequence = sequence.substr(seq_offset, fragment_size-1);
			bb_segs.push_back(std::make_pair(backbone_seg, loop_sequence));
		}
	}
	return bb_segs;
}



///@details create a loop size given the backbone segment from the loop library. Return the
///newly created pose.
std::pair<core::pose::Pose, core::Real>
LoophashAssemblyMover::build_loop_pose(
	core::pose::Pose const & pose,
	core::Size loop_anchor,
	core::Size /*n_segment_start*/,
	core::Size /*c_segment_end*/,
	core::Size n_overlapping_n,
	core::Size n_overlapping_c,
	protocols::loophash::BackboneSegment bb_seg,
	std::string loop_sequence
) const {

	using namespace protocols::loophash;

	core::Size fragment_size = bb_seg.length();
	core::Size num_new_residues = fragment_size - n_overlapping_n - n_overlapping_c - 1;
	//c_segment_end += num_new_residues;

	core::Size complete_loop_start = loop_anchor - n_overlapping_n;
	core::Size complete_loop_end = loop_anchor + num_new_residues + n_overlapping_c;

	core::pose::Pose combined_pose = pose;

	if( basic::options::option[ basic::options::OptionKeys::sewing::dump_pdbs ] ) {
		combined_pose.dump_pdb("start_pose.pdb");
	}

	//Append new loop residues to the pose
	core::chemical::ResidueTypeSetCOP rs = core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::CENTROID);
	for(core::Size i=0; i<num_new_residues; ++i) {
		core::conformation::ResidueOP new_rsd( NULL );
		core::chemical::ResidueTypeCOP type = rs->aa_map(core::chemical::aa_from_oneletter_code(loop_sequence[i+1+n_overlapping_n]))[1];
		if(TR.Debug.visible()) {
			TR.Debug << "Appending loophash residue " << type->name() << " (" << loop_sequence[i+1+n_overlapping_n] << ")" << std::endl;
		}
		new_rsd = core::conformation::ResidueFactory::create_residue( *type );
		combined_pose.conformation().safely_append_polymer_residue_after_seqpos(*new_rsd, loop_anchor+i, true);
	}
	core::pose::Pose unmodified_pose = combined_pose;

	if( basic::options::option[ basic::options::OptionKeys::sewing::dump_pdbs ] ) {
		combined_pose.dump_pdb("after_insertion.pdb");
	}

	TR << complete_loop_start << " " << complete_loop_end << " " << fragment_size << " " << n_overlapping_n << " " << n_overlapping_c << std::endl;
  LocalInserterOP inserter = new LocalInserter_SimpleMin();
	core::Real insertion_rmsd = inserter->make_local_bb_change_close_gaps(combined_pose, unmodified_pose, bb_seg, complete_loop_start);

	if( basic::options::option[ basic::options::OptionKeys::sewing::dump_pdbs ] ) {
		combined_pose.dump_pdb("after_apply.pdb");
	}

	//Build the pose
//
//	//Setup fold tree to create break
//	core::kinematics::FoldTree break_ft;
//	break_ft.add_edge(1, (int)complete_loop_end, core::kinematics::Edge::PEPTIDE);
//	break_ft.add_edge(1, (int)complete_loop_end+1, 1);
//	break_ft.add_edge((int)complete_loop_end+1, (int)combined_pose.total_residue(), core::kinematics::Edge::PEPTIDE);
//	combined_pose.fold_tree(break_ft);
//
//	//Add an extra residue such that torsions can be properly applied
//	TR << "Adding extra residue after position " << complete_loop_end << std::endl;
//	core::chemical::ResidueTypeCOP type = rs->aa_map(core::chemical::aa_from_oneletter_code('A'))[1];
//	core::conformation::ResidueOP extra_rsd = core::conformation::ResidueFactory::create_residue( *type );
//	combined_pose.conformation().safely_append_polymer_residue_after_seqpos(*extra_rsd, complete_loop_end, true);
//
//	//core::conformation::idealize_position(loop_anchor+num_new_residues, combined_pose.conformation());
//	for(core::Size i=complete_loop_start; i <= complete_loop_end; ++i) {
//		core::conformation::idealize_position(i, combined_pose.conformation());
//	}
//
//	//Apply loophash torsion angles
	std::vector<core::Real> phi   = bb_seg.phi();
	std::vector<core::Real> psi   = bb_seg.psi();
	std::vector<core::Real> omega = bb_seg.omega();
//	for ( core::Size i = 0; i < fragment_size; i++) {
//		core::Size pres = complete_loop_start+i;
//		if(TR.Debug.visible()) {
//			TR << "Updating res " << pres << " torsions:" << std::endl;
//			TR << "phi: " << combined_pose.phi(pres) << "->" << phi[i] << std::endl;
//			TR << "psi: " << combined_pose.psi(pres) << "->" << psi[i] << std::endl;
//			TR << "omega: " << combined_pose.omega(pres) << "->" << omega[i] << std::endl;
//		}
//		combined_pose.set_phi(pres, phi[i]);
//		combined_pose.set_psi(pres, psi[i]);
//		combined_pose.set_omega(pres, omega[i]);
//	}
//
//	//Now remove the extra residue
//	core::Size extra_resnum = complete_loop_end + 1;
//	TR << "Removing extra residue " << extra_resnum << std::endl;
//	combined_pose.conformation().delete_residue_slow(extra_resnum);
//
//	core::kinematics::FoldTree saved_ft = combined_pose.fold_tree();
//	core::kinematics::FoldTree simple_ft = core::kinematics::FoldTree((int)combined_pose.total_residue());
//	combined_pose.fold_tree(simple_ft);
//
//	//setup coordinate constraints
////	protocols::loops::Loops exclude_region;
////	exclude_region.add_loop( protocols::loops::Loop( complete_loop_start, res_pos + new_bs.length() ) );
////	transfer_phi_psi( original_pose, newpose );
////	add_coordinate_constraints_to_pose( newpose, original_pose, exclude_region );
//		protocols::relax::AtomCoordinateCstMover coord_cst_mover;
//		coord_cst_mover.apply(combined_pose);
//
//
//	core::kinematics::MoveMapOP move_map = new core::kinematics::MoveMap();
//	for(core::Size i=n_segment_start; i<=c_segment_end; ++i) {
//		TR << "Minimizing residue " << i << std::endl;
//		move_map->set_bb(i, true);
//		core::chemical::AtomIndices const & bb_indices = combined_pose.residue(i).mainchain_atoms();
//		for ( core::Size j = 1; j <= bb_indices.size(); ++j ) {
//			move_map->set( core::id::DOF_ID( core::id::AtomID( bb_indices[j], i ), core::id::D ), true );
//			move_map->set( core::id::DOF_ID( core::id::AtomID( bb_indices[j], i ), core::id::THETA ), true );
//			move_map->set( core::id::DOF_ID( core::id::AtomID( bb_indices[j], i ), core::id::PHI ), true );
//		}
//
//		if(i<pose.total_residue()-1){
//			TR.Debug << i << " bonded to " << i+1 << ": " << pose.residue(i).is_bonded(pose.residue(i+1)) << std::endl;
//		}
//	}
//
//	protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover();
//	min_mover->score_function(loop_scorefxn_);
//	min_mover->movemap(move_map);
//	min_mover->cartesian(true);

//	min_mover->score_function(loop_refine_scorefxn_);
//	min_mover->apply(combined_pose);

//	combined_pose.fold_tree(saved_ft);

	if(basic::options::option[basic::options::OptionKeys::sewing::dump_pdbs]) {
		Pose tmp;

		for(Size i = 1; i <= fragment_size+2; ++i) {
			core::conformation::ResidueOP new_rsd( NULL );
			new_rsd = core::conformation::ResidueFactory::create_residue( rs->name_map("ALA") );
			if(1==i) tmp.append_residue_by_jump( *new_rsd, 1 );
			else     tmp.append_residue_by_bond( *new_rsd, true );
			tmp.set_phi  ( tmp.n_residue(), 180.0 );
			tmp.set_psi  ( tmp.n_residue(), 180.0 );
			tmp.set_omega( tmp.n_residue(), 180.0 );
		}

		for ( Size i = 0; i < fragment_size; i++) {
			Size ires = 2+i;//due to the use of std:vector i has to start from 0, but positions offset by 1.
			tmp.set_phi  ( ires, phi[i]  );
			tmp.set_psi  ( ires, psi[i]  );
			tmp.set_omega( ires, omega[i]);
		}

		core::id::AtomID_Map< core::id::AtomID > atom_map( core::id::BOGUS_ATOM_ID );
		core::pose::initialize_atomid_map( atom_map, tmp, core::id::BOGUS_ATOM_ID );

//		atom_map.set( core::id::AtomID(tmp.residue( 2 ).atom_index("N"), 2), core::id::AtomID( pose.residue(complete_loop_start).atom_index("N"), complete_loop_start ) );
//		atom_map.set( core::id::AtomID(tmp.residue( 2 ).atom_index("CA"), 2 ), core::id::AtomID( pose.residue(complete_loop_start).atom_index("CA"), complete_loop_start ) );
//		atom_map.set( core::id::AtomID(tmp.residue( 2 ).atom_index("C"), 2 ), core::id::AtomID( pose.residue(complete_loop_start).atom_index("C"), complete_loop_start ) );

		for(core::Size i=0; i<n_overlapping_n+1; ++i) {
			core::Size frag_pos = 2+i;
			core::Size pose_pos = complete_loop_start+i;
			TR << "Aligning: " << frag_pos << " with " << pose_pos << std::endl;
			atom_map.set( core::id::AtomID(tmp.residue( frag_pos ).atom_index("N"), frag_pos), core::id::AtomID( pose.residue(pose_pos).atom_index("N"), pose_pos ) );
			atom_map.set( core::id::AtomID(tmp.residue( frag_pos ).atom_index("CA"), frag_pos ), core::id::AtomID( pose.residue(pose_pos).atom_index("CA"), pose_pos ) );
			atom_map.set( core::id::AtomID(tmp.residue( frag_pos ).atom_index("C"), frag_pos ), core::id::AtomID( pose.residue(pose_pos).atom_index("C"), pose_pos ) );
		}

		for(core::Size i=0; i<n_overlapping_c; ++i) {
			core::Size frag_pos = fragment_size+1-i;
			core::Size pose_pos = complete_loop_end-i;
			TR << "Aligning: " << frag_pos << " with " << pose_pos << std::endl;
			atom_map.set( core::id::AtomID(tmp.residue( frag_pos ).atom_index("N"), frag_pos ), core::id::AtomID( pose.residue(pose_pos).atom_index("N"), pose_pos ) );
			atom_map.set( core::id::AtomID(tmp.residue( frag_pos ).atom_index("CA"), frag_pos ), core::id::AtomID( pose.residue(pose_pos).atom_index("CA"), pose_pos ) );
			atom_map.set( core::id::AtomID(tmp.residue( frag_pos ).atom_index("C"), frag_pos ), core::id::AtomID( pose.residue(pose_pos).atom_index("C"), pose_pos ) );
		}

		core::Real rms = core::scoring::superimpose_pose(tmp, pose/*const*/, atom_map);
		TR << "Overlapping N: " << n_overlapping_n << std::endl;
		TR << "Overlapping C: " << n_overlapping_c << std::endl;
		TR << "LH Fragment size: " << fragment_size << std::endl;
		TR << "Num new residues: " << num_new_residues << std::endl;
		TR << "RMS of flanking region: " << rms << std::endl;
		tmp.dump_pdb("bb_seg.pdb");
	}

	return std::make_pair(combined_pose, insertion_rmsd);
}



void
LoophashAssemblyMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
){
	using namespace basic::options;

	AssemblyMover::parse_my_tag(tag, data, filters, movers, pose);
	
	if(tag->hasOption("max_loop_segments_to_build")) {
		max_loop_segments_to_build_ = tag->getOption<core::Size>("max_loop_segments_to_build");
	}

	if(tag->hasOption("max_loop_distance")) {
		max_loop_distance_ = tag->getOption<core::Real>("max_loop_distance");
	}

	if(tag->hasOption("max_insertion_rms")) {
		max_insertion_rms_ = tag->getOption<core::Real>("max_insertion_rms");
	}
}

} //sewing
} //devel
