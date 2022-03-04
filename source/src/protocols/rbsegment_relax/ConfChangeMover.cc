// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @file protocols/rbsegment_relax/ConfChangeMover.cc
/// @brief Sampling of new conformations from an input structure
/// @author Diego del Alamo (diego.delalamo@gmail.com) and Davide Sala (d.sala1388@gmail.com)

#include <protocols/rbsegment_relax/ConfChangeMover.hh>
#include <protocols/rbsegment_relax/ConfChangeMoverCreator.hh>

// RosettaScripts inputs
#include <basic/datacache/DataMap.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/Tracer.hh>
#include <boost/algorithm/string.hpp>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/IndependentBBTorsionSRFD.hh>
#include <core/fragment/picking_old/vall/util.hh>
#include <core/id/AtomID.hh>
#include <core/id/SequenceMapping.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/RT.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/func/FlatHarmonicFunc.hh>
#include <core/scoring/func/CircularSigmoidalFunc.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/types.hh>
#include <core/util/kinematics_util.hh>

#include <protocols/hybridization/CartesianHybridize.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rbsegment_relax/RBSegment.hh>
#include <protocols/rbsegment_relax/util.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/GunnCost.hh>

#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector1.hh>

// Numeric Headers
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <string>


namespace protocols {
namespace rbsegment_relax {

static basic::Tracer TR("protocols.rbsegment_relax.ConfChangeMover");

///////////////////////////////////////////////////////////////////////////////

/// @brief Constructor
ConfChangeMover::ConfChangeMover() : hybridization::CartesianHybridize() {
}

/// @brief Copy function - inherited
moves::MoverOP ConfChangeMover::clone() const {
	return utility::pointer::make_shared< ConfChangeMover >(*this);
}

///////////////////////////////////////////////////////////////////////////////
/// @brief  Pose manipulation for conformational changes
/// @details Sets up a "star" topology and obtains rigid body segments to
///      move around. A copy of the pose is then made and its rigid body
///      segments are moved in Cartesian space to observe the effects
///      of a rigid body transform. The loops are then closed
void ConfChangeMover::apply(core::pose::Pose &pose) {
	// MAKE SURE PROTEIN IS IN CENTROID MODE FOR THIS SAMPLING
	simple_moves::SwitchResidueTypeSetMover to_centroid("centroid");
	to_centroid.apply(pose);
	// MAKE A COPY OF THE ORIGINAL POSE FOR THE FOLLOWING REASONS:
	//  A. RETAINS ORIGINAL FOLD TREE FOR LOOP CLOSURE
	//  B. SOURCE OF FRAGMENTS FOR LOOP CLOSURE
	//CHECK OR CREATE FRAGS
	check_or_create_fragments(pose);
	auto original_pose = pose;
	auto original_foldtree = pose.fold_tree();

	// RIGID BODY MOVEMENT
	// VIRTUAL RESIDUE IS ESSENTIAL FOR "STAR TOPOLOGY"
	auto cst_set = pose.constraint_set()->clone();
	add_dihedral_csts_to_rb(pose, (template_) ? *template_ : original_pose);
	core::pose::addVirtualResAsRoot(pose);
	stage1(pose);
	core::pose::remove_virtual_residues(pose);
	pose.constraint_set(cst_set);

	core::util::remove_cutpoint_variants(pose);
	pose.fold_tree(original_foldtree);

	if ( stage2_moves_ > 0 && stage2_models_ == 1 ) {
		stage2(pose, original_pose);

	}

	if ( stage2_moves_ > 0 && stage2_models_ > 1 ) {

		TR.Warning << "Movers downstream to ConfChangeMover must be located within <MultiplePoseMover> ROSETTASCRIPTS subtags." << std::endl;

		for ( core::Size it = 1; it <= stage2_models_; ++it ) {
			store_poseOPs_.push_back(pose.clone());
			TR << std::endl;
			TR << "///////////////// Begin Stage2 iteration number " << it << " //////////////////////" << std::endl;
			TR << std::endl;
			stage2(*store_poseOPs_.back(), original_pose);
		}

		stage2_scorefxn_->set_weight(core::scoring::dihedral_constraint, 0.0);
		stage2_scorefxn_->set_weight(core::scoring::coordinate_constraint, 0.0);
		TR << "Total Number of expected poses in output:  " << store_poseOPs_.size() << std::endl;

		pose = *store_poseOPs_.back();
		TR << "Returning final pose" << std::endl;
		store_poseOPs_.pop_back();
	}
}

/// @brief necessary function to output multiple poses
core::pose::PoseOP ConfChangeMover::get_additional_output() {

	core::pose::PoseOP out_pose_OP;
	if ( store_poseOPs_.size() == 0 ) {
		return nullptr;
	}
	out_pose_OP = store_poseOPs_.back();
	store_poseOPs_.pop_back();

	TR << "Returning additional pose" << std::endl;

	return out_pose_OP;
}


/// @brief check if fragments (3 and 9) are presents or pick them
void ConfChangeMover::check_or_create_fragments(core::pose::Pose &pose) {

	std::string tgt_seq = pose.sequence();
	core::scoring::dssp::Dssp dssp(pose);
	dssp.insert_ss_into_pose(pose);
	std::string tgt_ss = dssp.get_dssp_secstruct();

	// pick from vall based on template SS and target sequence
	if ( !frags9_ ) {

		frags9_ = utility::pointer::make_shared<core::fragment::ConstantLengthFragSet>(9);

		for ( core::Size j = 1; j <= pose.total_residue() - 9 - 1; ++j ) {
			std::string ss_sub = tgt_ss.substr(j - 1, 9);
			std::string aa_sub = tgt_seq.substr(j - 1, 9);
			core::fragment::FrameOP frame(new core::fragment::Frame(j, 9));
			frame->add_fragment(core::fragment::picking_old::vall::pick_fragments_by_ss_plus_aa(ss_sub, aa_sub, 200, true, core::fragment::IndependentBBTorsionSRFD()));
			frags9_->add(frame);
		}
	}

	if ( !frags3_ ) {

		frags3_ = utility::pointer::make_shared<core::fragment::ConstantLengthFragSet>(3);
		for ( core::Size j = 1; j <= pose.total_residue() - 3 - 1; ++j ) {
			std::string ss_sub = tgt_ss.substr(j - 1, 3);
			std::string aa_sub = tgt_seq.substr(j - 1, 3);
			core::fragment::FrameOP frame(new core::fragment::Frame(j, 3));
			frame->add_fragment(core::fragment::picking_old::vall::pick_fragments_by_ss_plus_aa(ss_sub, aa_sub, 200, true, core::fragment::IndependentBBTorsionSRFD()));
			frags3_->add(frame);
		}
	}
}

/// @brief SETUP TOPOLOGY FOR RIGID BODY MOVES
void ConfChangeMover::stage1_pose_setup(core::pose::Pose &pose) {
	// IF AN RB FILE IS PROVIDED, USE IT
	core::kinematics::MoveMapOP mm1(new core::kinematics::MoveMap());
	if ( rb_file_ == "AUTO" ) {
		rbsegment_relax::guess_rbsegs_from_pose(pose, rbsegs_, rbsegs_2_, loops_);
	} else {
		rbsegment_relax::read_RBSegment_file(rbsegs_, loops_, rb_file_);
	}
	rbsegment_relax::setup_pose_rbsegs_keep_loops(pose, rbsegs_, loops_, mm1);

	// IF RESIDUE SELECTORS ARE PROVIDED, USE THEM
	if ( s1_res_ ) {
		// subset is a vector of size pose.size() with booleans
		auto const subset = s1_res_->apply(pose);
		utility::vector1<bool> use_rb(rbsegs_.size(), false);

		core::Size save_last = 0;

		if ( !seg_.empty() ) {

			for ( core::Size i = 1; i <= rbsegs_.size(); ++i ) {

				//setting new start and stop after matching seg_.front() with the i-sh segment
				if ( !seg_.empty() && i == seg_.front() ) {
					rbsegs_[i][1].set_start(newstart_.front());
					newstart_.erase(newstart_.begin());
					rbsegs_[i][1].set_end(newend_.front());
					newend_.erase(newend_.begin());
					seg_.erase(seg_.begin());

				}

				// if segment start is lower than previous segment end must be removed because of overlapping residues
				if ( rbsegs_[i][1].start() <= save_last ) {
					rbsegs_.erase(rbsegs_.begin() + i - 1);
					for ( core::Size it = 1; it <= seg_.size(); ++it ) {
						seg_[it] -= 1;
					}
					i = i-1;
				} else {
					save_last = rbsegs_[i][1].end();
				}
			}
		}

		//add new segments
		if ( !addstart_.empty() || !addend_.empty() ) {
			for ( core::Size i = rbsegs_.size() + 1; i <= (rbsegs_.size() + addstart_.size()); ++i ) {
				rbsegs_.push_back(rbsegs_[1]);
				rbsegs_[i][1].set_start(addstart_.front());
				rbsegs_[i][1].set_end(addend_.front());
				addstart_.erase(addstart_.begin());
				addend_.erase(addend_.begin());
			}
		}

		// set true for segment (residues) to move
		for ( core::Size i = 1; i <= rbsegs_.size(); ++i ) {
			//printing and setting residues to keep
			TR << "Segment " << i << ": ";
			for ( core::Size res = rbsegs_[i][1].start(); res <= rbsegs_[i][1].end(); ++res ) {
				TR << res << " ";
				if ( subset[res] ) {
					use_rb[i] = true;
				}
			}
			if ( use_rb[i] ) {
				TR << "TRUE";
			} else {
				TR << "FALSE";
			}
			TR << std::endl;
		}

		// remove segments out of the residue selector range
		for ( core::Size i = rbsegs_.size(); i >= 1; --i ) {
			if ( !use_rb[i] ) {
				rbsegs_.erase(rbsegs_.begin() + i - 1);
			}
		}
	}
	if ( rbsegs_.size() <= 2 ) {
		stage1_multi_freq_ = 0.0;
	}

	//store rbsegs to be restored after domain movement
	stored_rebsegs_ = rbsegs_;

	// Determine rigid bodies that should be moved together
	if ( s1_rigid_res_.size() ) {
		std::set<core::Size> all_segs;
		utility::vector1<std::set<core::Size> > rigid_segs;
		for ( auto const &selector : s1_rigid_res_ ) {
			TR << "Selector:" << std::endl;
			TR << "\t";
			auto const subset = selector->apply(pose);
			std::set<core::Size> segs_in_rb;
			for ( core::Size i = 1; i <= rbsegs_.size(); ++i ) {
				for ( auto const &res : get_residues_in_rbsegment(rbsegs_[i]) ) {
					if ( subset[res] ) {
						TR << i << ":" << res << " ";
						segs_in_rb.insert(i);
						all_segs.insert(i);
					}
				}
			}
			TR << std::endl;
			rigid_segs.push_back(segs_in_rb);
		}
		utility::vector1<rbsegment_relax::RBSegment> new_segs;
		for ( core::Size i = 1; i <= rbsegs_.size(); ++i ) {
			if ( all_segs.find(i) == all_segs.end() ) {
				new_segs.push_back(rbsegs_[i]);
			}
		}
		for ( core::Size i = 1; i <= rigid_segs.size(); ++i ) {
			utility::vector1<rbsegment_relax::RBSegment> new_segs_local;
			for ( auto const &seg : rigid_segs[i] ) {
				for ( core::Size j = 1; j <= rbsegs_[seg].nContinuousSegments(); ++j ) {
					new_segs_local.push_back(rbsegs_[seg][j]);
				}

			}
			new_segs.push_back(rbsegment_relax::RBSegment(new_segs_local));
		}
		rbsegs_ = new_segs;
	}

	TR << "Final segments to move:" << std::endl;
	for ( core::Size i = 1; i <= rbsegs_.size(); ++i ) {
		//printing and setting residues to keep
		TR << "Segment " << i << ": ";
		std::set<core::Size> residues_in_seg;
		recursive_residues_from_rbsegs(rbsegs_[i], residues_in_seg);
		for ( auto const &res : residues_in_seg ) {
			TR << res << " ";
		}
		TR << std::endl;
	}
}


/// @brief APPLY RIGID BODY MOVES
void ConfChangeMover::stage1(core::pose::Pose &pose) {

	// VARS, IN ORDER OF APPEARANCE
	core::id::SequenceMapping resmap;

	core::kinematics::MoveMapOP mm1(new core::kinematics::MoveMap());
	mm1->set_bb(true);

	if ( s1_res_ ) {
		TR.Warning << "Any selected residue belonging to a segment will activate such segment. Segments residue composition can be checked out in the stdout. " << std::endl;
		TR.Warning << "For each segment, a TRUE tag at the end of the residues range points out its activation." << std::endl;
		auto const subset = s1_res_->apply(pose);
		for ( core::Size it = 1; it <= subset.size(); ++it ) {
			if ( !subset[it] ) {
				mm1->set_bb(it, false);
			}
		}
	}

	simple_moves::SmoothFragmentMover smoothmover1(frags3_, mm1, utility::pointer::make_shared<simple_moves::GunnCost>());
	simple_moves::SmoothFragmentMover smoothmover2(frags9_, mm1, utility::pointer::make_shared<simple_moves::GunnCost>());
	core::optimization::MinimizerOptions options_5("linmin", 0.01, true, false, false);
	auto min_sfxn = core::scoring::ScoreFunctionFactory::create_score_function("score4_smooth_cart");
	auto vdw_round = 3 * stage1_moves_ / 4;

	// VAR SETUP
	smoothmover1.set_check_ss(true);
	smoothmover1.enable_end_bias_check(false);
	smoothmover2.set_check_ss(true);
	smoothmover2.enable_end_bias_check(false);
	options_5.max_iter(200);

	///////////////////
	// MOVER SETUP
	stage1_pose_setup(pose);

	// LET'S GO
	TR << "ConfChange Round 1: Rigid-body movement" << std::endl;
	auto mc = moves::MonteCarloOP(new moves::MonteCarlo(pose, *stage1_scorefxn_, stage1_temp_));
	mc->set_temperature(stage1_temp_);


	for ( core::Size iter = 1; iter <= stage1_moves_; ++iter ) {
		if ( iter == stage1_move_restore_ ) {
			rbsegs_ = stored_rebsegs_;
		}
		auto move_n = numeric::random::uniform();
		auto move = (move_n < stage1_frag_freq_) ? s1_frag_move_ : (1.0 - move_n > 1.0 - stage1_multi_freq_) ? s1_rb_multi_ : s1_rb_single_;
		if ( move == s1_frag_move_ ) {
			if ( numeric::random::uniform() > 0.5 ) {
				smoothmover1.apply(pose);
			} else {
				smoothmover2.apply(pose);
			}
		} else {
			auto seg = multi_rb(pose, (move == s1_rb_multi_));
			if ( seg.isHelix() && seg.isSimple() && numeric::random::uniform() < stage1_twist_freq_ ) {
				move = s1_twist_;
				helixmover_.set_segment(seg);
				helixmover_.apply(pose);
			} else {
				rbmover_.set_segment(seg);
				rbmover_.apply(pose);
			}
		}
		// CHECK GAPS. IF SSES ARE TOO FAR APART, THEN RESTORE PSOE
		if ( !check_gaps(pose) ) {
			pose = mc->last_accepted_pose();
			continue;
		}
		// CHECK MINIMIZATION
		if ( stage1_minimization_ && iter >= vdw_round ) {
			minimizer_.run(pose, *mm1, *min_sfxn, options_5);
		}
		// APPLY THE MOVE USING MCM
		mc->boltzmann(pose, move);
		if ( iter % 1000 == 0 ) {
			TR << "SSEs movement" << std::endl;
			mc->show_scores();
			mc->show_counters();
		}
	}
	stage1_scorefxn_->show(TR, pose);
}

/// @brief Checks that no unbridgeable gaps were introduced during sampling
bool ConfChangeMover::check_gaps(core::pose::Pose &pose) {
	// GAP FXN IS FROM BCL::SCORE PAPER
	for ( core::Size i = 1; i < rbsegs_.size(); ++i ) {
		auto last_res_i = *(get_residues_in_rbsegment(rbsegs_[i]).rbegin());
		auto first_res_j = *(get_residues_in_rbsegment(rbsegs_[i + 1]).begin());
		auto dist = pose.residue(last_res_i).atom("C").xyz().distance(pose.residue(first_res_j).atom("N").xyz());
		auto max_gap = 2.11 + 2.56 * (first_res_j - last_res_i);
		if ( dist > max_gap ) {
			return false;
		}
	}
	return true;
}


/// @brief Which SSEs to move if more than one is to be moved
rbsegment_relax::RBSegment ConfChangeMover::multi_rb(core::pose::Pose &pose, bool const &multiple_sses) {
	// FYI SSE MEANS SECONDARY STRUCTURAL ELEMENT
	std::set<int> sses{numeric::random::random_range(1, rbsegs_.size())};
	if ( multiple_sses && rbsegs_.size() > 2 ) {
		auto const contacts = sse_contact_strength(pose);
		auto const n_sses = numeric::random::random_range(2, rbsegs_.size() - 1);
		while ( int(sses.size()) < int(n_sses) ) {
			core::Real best_contact_str = 0.0;
			core::Size sse_w_best_contact = 0;
			for ( auto const &sse1 : sses ) {
				for ( core::Size i = 1; i <= contacts.at(sse1).size(); ++i ) {
					if ( sses.find(i) != sses.end() ) {
						continue;
					} else if ( contacts.at(sse1)[i] > best_contact_str ) {
						sse_w_best_contact = i;
					}
				}
			}
			if ( best_contact_str == 0 ) {
				break;
			}
			sses.insert(sse_w_best_contact);
		}
	}
	utility::vector1<rbsegment_relax::RBSegment> output;
	for ( auto const &it : sses ) {
		for ( auto it2 : rbsegs_[it].split() ) {
			output.push_back(it2);
		}
	}
	return rbsegment_relax::RBSegment(output);
}

utility::vector1<utility::vector1<core::Real>> ConfChangeMover::sse_contact_strength(core::pose::Pose &pose) const {
	utility::vector1<utility::vector1<core::Real>> output(rbsegs_.size(), utility::vector1<core::Real>(rbsegs_.size(), 0.0));

	pose.update_residue_neighbors();
	auto const &neighborgraph = pose.energies().tenA_neighbor_graph();
	for ( core::Size i = 1; i < rbsegs_.size(); ++i ) {
		auto res_in_i = get_residues_in_rbsegment(rbsegs_[i]);
		for ( core::Size j = i + 1; j <= rbsegs_.size(); ++j ) {
			auto contacts = 0;
			auto res_in_j = get_residues_in_rbsegment(rbsegs_[j]);
			auto const &smaller_sse = (res_in_i.size() < res_in_j.size()) ? res_in_i : res_in_j;
			auto const &larger_sse = (res_in_i.size() > res_in_j.size()) ? res_in_i : res_in_j;
			for ( auto const &res : smaller_sse ) {
				//core::Size res = pose_to_seg_id.at( res_pdb ); // Pesky pose numbering
				for ( auto edge_it = neighborgraph.get_node(res)->const_edge_list_begin(); edge_it != neighborgraph.get_node(res)->const_edge_list_end(); ++edge_it ) {
					auto const &node_index = (*edge_it)->get_other_node(res)->get_node_index();
					if ( pose.residue(node_index).is_virtual_residue() ) continue;
					if ( larger_sse.find(node_index) != larger_sse.end() ) {
						contacts += 1;
						break;
					}
				}
			}
			output[i][j] = core::Real(contacts) / core::Real(smaller_sse.size());
			output[j][i] = output[i][j];
		}
	}
	return output;
}

void ConfChangeMover::recursive_residues_from_rbsegs(rbsegment_relax::RBSegment const &rbseg, std::set<core::Size> &residues) const {
	if ( rbseg.nContinuousSegments() == 1 ) {
		for ( auto res = rbseg[1].start(); res <= rbseg[1].end(); ++res ) {
			residues.insert(res);
		}
	} else {
		for ( core::Size i = 1; i <= rbseg.nContinuousSegments(); ++i ) {
			recursive_residues_from_rbsegs(rbseg[i], residues);
		}
	}
}

std::set<core::Size> ConfChangeMover::get_residues_in_rbsegment(rbsegment_relax::RBSegment const &rbseg) const {

	std::set<core::Size> output;
	for ( auto x = 1; x <= int(rbseg.nContinuousSegments()); ++x ) {
		for ( auto res = rbseg[x].start(); res <= rbseg[x].end(); ++res ) {
			output.insert(res);
		}
	}
	return output;
}

/// @brief fragments insertion from either input structure or from picked
void ConfChangeMover::stage2(core::pose::Pose &pose, core::pose::Pose const &original_pose) {
	TR << "Stage2 activated" << std::endl;

	auto cst_set = pose.constraint_set()->clone();
	add_dihedral_csts(pose, original_pose);
	core::optimization::MinimizerOptions options_short("lbfgs_armijo_nonmonotone", 0.01, true, false, false);
	options_short.max_iter(5);
	core::optimization::MinimizerOptions options_lbfgs("lbfgs_armijo_nonmonotone", 0.01, true, false, false);
	options_lbfgs.max_iter(2000);
	auto min_sfxn = core::scoring::ScoreFunctionFactory::create_score_function("score4_smooth_cart");
	//auto min_sfxn = stage2_scorefxn_;
	core::Real dihedral_wt = stage2_scorefxn_->get_weight(core::scoring::dihedral_constraint);
	core::Real cart_wt = stage2_scorefxn_->get_weight(core::scoring::cart_bonded);
	if ( cart_wt == 0.0 ) {
		TR.Error << "Scoreterm cart_bonded must be nonzero! Setting to 0.05." << std::endl;
		cart_wt = 0.05;
	}


	if ( dihedral_wt == 0.0 ) {
		TR.Warning << "Scoreterm dihedral_constraint set to zero! Do not do this unless you know what you are doing!" << std::endl;
	}

	min_sfxn->set_weight(core::scoring::coordinate_constraint, 10.0);
	min_sfxn->set_weight(core::scoring::dihedral_constraint, dihedral_wt);
	min_sfxn->set_weight(core::scoring::cart_bonded, 0.5);

	core::kinematics::MoveMap mm;
	mm.set_bb(false);
	mm.set_chi(false);
	mm.set_jump(false);


	utility::vector1<loops::Loop> contigs;
	core::Size MIN_SIZE = 5;
	core::Size MAX_SIZE = 15;

	utility::vector1<core::Size> s2_size;
	core::Size start;
	core::Size stop;

	if ( s2_res_ ) {
		auto const subset = s2_res_->apply(pose);
		for ( core::Size it = 1; it <= subset.size(); ++it ) {
			if ( subset[it] ) {
				s2_size.push_back(it);
			}
		}
		start = s2_size.front();
		stop = s2_size.back();
	} else {
		start = 1;
		stop = pose.size();
	}

	for ( core::Size first = start; first <= stop - MIN_SIZE + 1; ++first ) {
		for ( core::Size last = first + MIN_SIZE - 1; last <= std::min(first + MAX_SIZE - 1, stop); ++last ) {
			contigs.push_back(loops::Loop(first, last));
		}
	}

	// INIT FRAMES
	for ( core::fragment::ConstFrameIterator it = frags9_->begin(); it != frags9_->end(); ++it ) {
		frames_[(*it)->start()] = **it;
	}

	auto mc = moves::MonteCarloOP(new moves::MonteCarlo(pose, *stage2_scorefxn_, stage2_temp_));

	core::Size r1 = stage2_moves_ / 2;

	stage2_scorefxn_->set_weight(core::scoring::cart_bonded, cart_wt / 10.0);
	stage2_scorefxn_->set_weight(core::scoring::coordinate_constraint, 1.0);
	min_sfxn->set_weight(core::scoring::dihedral_constraint, dihedral_wt / 10.0);

	TR << "ConfChange Round 2: Fragment hybridization" << std::endl;
	for ( core::Size iter = 1; iter <= stage2_moves_; ++iter ) {
		mm.set_bb(false);
		auto move = (numeric::random::uniform() < stage2_segment_freq_) ? s2_homolog_move_ : ((numeric::random::uniform() < stage2_targgaps_freq_) ? s2_frame_select_
			: s2_frame_random_);
		if ( iter == r1 ) {
			TR << "Resetting scores!" << std::endl;
			stage2_scorefxn_->set_weight(core::scoring::cart_bonded, cart_wt);
			stage2_scorefxn_->set_weight(core::scoring::dihedral_constraint, dihedral_wt);
			mc->recover_low(pose);
			mc->reset(pose);
			mc->reset_scorefxn(pose, *stage2_scorefxn_);
		}

		if ( move == s2_homolog_move_ || move == s2_frame_select_ ) {
			// Here we do the gap calculation
			utility::vector1<core::Real> biggest_gaps{0.0, 0.0, 0.0};
			utility::vector1<core::Size> res_w_gaps{0, 0, 0};


			for ( core::Size res = start; res < stop; ++res ) {
				core::Real dist = std::abs(pose.residue(res).atom("C").xyz().distance(pose.residue(res + 1).atom("N").xyz()) - 1.328685);
				if ( dist > biggest_gaps[1] ) {
					biggest_gaps[3] = biggest_gaps[2];
					biggest_gaps[2] = biggest_gaps[1];
					biggest_gaps[1] = dist;
					res_w_gaps[3] = res_w_gaps[2];
					res_w_gaps[2] = res_w_gaps[1];
					res_w_gaps[1] = res;
				} else if ( dist > biggest_gaps[2] ) {
					biggest_gaps[3] = biggest_gaps[2];
					biggest_gaps[2] = dist;
					res_w_gaps[3] = res_w_gaps[2];
					res_w_gaps[2] = res;
				} else if ( dist > biggest_gaps[3] ) {
					biggest_gaps[3] = dist;
					res_w_gaps[3] = res;
				}
			}
			numeric::random::random_permutation(res_w_gaps, numeric::random::rg());
			if ( move == s2_homolog_move_ ) {
				bool found_loop = false;
				for ( auto const &res_to_close : res_w_gaps ) {
					numeric::random::random_permutation(contigs, numeric::random::rg());
					for ( auto &loop : contigs ) {
						if ( loop.start() < res_to_close && loop.stop() > res_to_close ) {
							for ( core::Size r = std::max(1, int(loop.start() - 1)); r <= std::min(pose.size(), loop.stop() + 1); ++r ) {
								mm.set_bb(r, true);
							}
							auto temp_pose = original_pose;
							apply_frag(pose, temp_pose, loop);
							found_loop = true;
							break;
						}
					}
					if ( found_loop ) break;
				}

				if ( !found_loop ) {
					core::Size res = numeric::random::random_range(start, stop - 8);
					for ( core::Size r = std::max(1, int(res - 1)); r <= std::min(stop, res + 9); ++r ) {
						mm.set_bb(r, true);
					}
					apply_frame(pose, frames_[res]);
				}
			}
			if ( move == s2_frame_select_ ) {
				numeric::random::random_permutation(res_adj_, numeric::random::rg());
				core::Size res = std::max(start, std::min(stop - 8, res_w_gaps[1] - res_adj_[1]));
				for ( core::Size r = std::max(1, int(res - 1)); r <= std::min(stop, res + 9); ++r ) {
					mm.set_bb(r, true);
				}
				apply_frame(pose, frames_[res]);
			}
		} else {
			core::Size res = numeric::random::random_range(start, stop - 8);
			for ( core::Size r = std::max(1, int(res - 1)); r <= std::min(stop, res + 9); ++r ) {
				mm.set_bb(r, true);
			}
			apply_frame(pose, frames_[res]);
		}
		if ( iter >= r1 || move == s2_homolog_move_ ) {
			minimizer_.run(pose, mm, *min_sfxn, options_short);
		}
		mc->boltzmann(pose, move);
		if ( iter % 50 == 0 ) {
			mc->show_scores();
			mc->show_counters();
		}
	}

	mc->recover_low(pose);
	stage2_scorefxn_->show(TR, pose);

	if ( s2_res_ ) {
		for ( auto const &res : s2_size ) {
			mm.set_bb(res, true);
			mm.set_jump(res, true);
		}
	} else {
		mm.set_bb(true);
		mm.set_jump(true);
	}

	if ( stage2_models_ == 1 ) {
		stage2_scorefxn_->set_weight(core::scoring::dihedral_constraint, 0.0);
		stage2_scorefxn_->set_weight(core::scoring::coordinate_constraint, 0.0);
	}
	minimizer_.run(pose, mm, *min_sfxn, options_lbfgs);
	stage2_scorefxn_->show(TR, pose);
	pose.constraint_set(cst_set);


}

/// @brief add dihedral constraints to rigid-body moves
void ConfChangeMover::add_dihedral_csts_to_rb(core::pose::Pose &pose, core::pose::Pose const &original_pose) {
	TR << "add_dihedral_csts_to_rb activated" << std::endl;
	core::scoring::constraints::ConstraintSetOP cst_set = pose.constraint_set()->clone();

	std::set<core::Size> rbseg_residues;
	for ( auto const &rbseg : rbsegs_ ) {
		auto residues = get_residues_in_rbsegment(rbseg);
		for ( auto const &res : residues ) {
			rbseg_residues.insert(res);
		}
	}
	for ( core::Size res = 1; res <= pose.size(); ++res ) {
		if ( rbseg_residues.find(res) != rbseg_residues.end() ) {
			if ( res > 1 ) {
				core::scoring::func::FuncOP f_phi(new core::scoring::func::CircularSigmoidalFunc(original_pose.phi(res), numeric::constants::d::pi_over_3, 0.25));
				cst_set->add_constraint(core::scoring::constraints::ConstraintCOP(
					utility::pointer::make_shared<core::scoring::constraints::DihedralConstraint>(core::id::AtomID(pose.residue(res - 1).atom_index("C"), res - 1),
					core::id::AtomID(pose.residue(res).atom_index("N"), res),
					core::id::AtomID(pose.residue(res).atom_index("CA"), res),
					core::id::AtomID(pose.residue(res).atom_index("C"), res), f_phi)));
			}
			if ( res < pose.size() ) {
				core::scoring::func::FuncOP f_psi(new core::scoring::func::CircularSigmoidalFunc(original_pose.psi(res), numeric::constants::d::pi_over_3, 0.25));
				cst_set->add_constraint(core::scoring::constraints::ConstraintCOP(
					utility::pointer::make_shared<core::scoring::constraints::DihedralConstraint>(core::id::AtomID(pose.residue(res).atom_index("N"), res),
					core::id::AtomID(pose.residue(res).atom_index("CA"), res),
					core::id::AtomID(pose.residue(res).atom_index("C"), res),
					core::id::AtomID(pose.residue(res + 1).atom_index("N"), res + 1), f_psi)));
			}
		}
	}
	pose.constraint_set(cst_set);
}

/// @brief add dihedral constraints to stage2
void ConfChangeMover::add_dihedral_csts(core::pose::Pose &pose, core::pose::Pose const &original_pose) {
	TR << "add_dihedral_csts activated" << std::endl;
	core::scoring::constraints::ConstraintSetOP cst_set = pose.constraint_set()->clone();

	core::Size const &anchor = pose.fold_tree().root();
	std::set<core::Size> rbseg_residues;
	for ( auto const &rbseg : rbsegs_ ) {
		auto residues = get_residues_in_rbsegment(rbseg);
		for ( auto const &res : residues ) {
			rbseg_residues.insert(res);
		}
	}


	utility::vector1<bool> dihed(pose.size(), false);
	if ( dihed2_res_ ) {
		auto const subset = dihed2_res_->apply(pose);
		TR << "Dihedral constraints will be removed on the following residues: ";
		for ( core::Size it = 1; it <= pose.size(); ++it ) {
			if ( subset[it] ) {
				dihed[it] = true;
				TR << it << " ";
			}
		}
		TR << std::endl;
	}

	for ( core::Size res = 1; res <= pose.size(); ++res ) {
		if ( rbseg_residues.find(res) != rbseg_residues.end() ) {
			for ( core::Size atom = 1; atom <= pose.residue(res).last_backbone_atom(); ++atom ) {
				core::scoring::func::FuncOP f(new core::scoring::func::FlatHarmonicFunc(0.0, 1.0, 2.0));
				cst_set->add_constraint(core::scoring::constraints::ConstraintCOP(
					utility::pointer::make_shared<core::scoring::constraints::CoordinateConstraint>(core::id::AtomID(atom, res), core::id::AtomID(1, anchor), pose.residue(res).xyz(atom),
					f)));
			}
			if ( res > 1 && !dihed[res] ) {
				core::scoring::func::FuncOP f_phi(new core::scoring::func::CircularSigmoidalFunc(pose.phi(res), numeric::constants::d::pi_over_3, 0.25));
				cst_set->add_constraint(core::scoring::constraints::ConstraintCOP(
					utility::pointer::make_shared<core::scoring::constraints::DihedralConstraint>(core::id::AtomID(pose.residue(res - 1).atom_index("C"), res - 1),
					core::id::AtomID(pose.residue(res).atom_index("N"), res),
					core::id::AtomID(pose.residue(res).atom_index("CA"), res),
					core::id::AtomID(pose.residue(res).atom_index("C"), res), f_phi)));
			}
			if ( res < pose.size() && !dihed[res] ) {
				core::scoring::func::FuncOP f_psi(new core::scoring::func::CircularSigmoidalFunc(pose.psi(res), numeric::constants::d::pi_over_3, 0.25));
				cst_set->add_constraint(core::scoring::constraints::ConstraintCOP(
					utility::pointer::make_shared<core::scoring::constraints::DihedralConstraint>(core::id::AtomID(pose.residue(res).atom_index("N"), res),
					core::id::AtomID(pose.residue(res).atom_index("CA"), res),
					core::id::AtomID(pose.residue(res).atom_index("C"), res),
					core::id::AtomID(pose.residue(res + 1).atom_index("N"), res + 1), f_psi)));
			}
		} else {
			if ( res > 1 && !dihed[res] ) {
				core::scoring::func::FuncOP f_phi(new core::scoring::func::CircularSigmoidalFunc(original_pose.phi(res), numeric::constants::d::pi_over_3, 0.25));
				cst_set->add_constraint(core::scoring::constraints::ConstraintCOP(
					utility::pointer::make_shared<core::scoring::constraints::DihedralConstraint>(core::id::AtomID(pose.residue(res - 1).atom_index("C"), res - 1),
					core::id::AtomID(pose.residue(res).atom_index("N"), res),
					core::id::AtomID(pose.residue(res).atom_index("CA"), res),
					core::id::AtomID(pose.residue(res).atom_index("C"), res), f_phi)));
			}
			if ( res < pose.size() && !dihed[res] ) {
				core::scoring::func::FuncOP f_psi(new core::scoring::func::CircularSigmoidalFunc(original_pose.psi(res), numeric::constants::d::pi_over_3, 0.25));
				cst_set->add_constraint(core::scoring::constraints::ConstraintCOP(
					utility::pointer::make_shared<core::scoring::constraints::DihedralConstraint>(core::id::AtomID(pose.residue(res).atom_index("N"), res),
					core::id::AtomID(pose.residue(res).atom_index("CA"), res),
					core::id::AtomID(pose.residue(res).atom_index("C"), res),
					core::id::AtomID(pose.residue(res + 1).atom_index("N"), res + 1), f_psi)));
			}
		}
	}
	pose.constraint_set(cst_set);
}

void ConfChangeMover::parse_my_tag(utility::tag::TagCOP tag, basic::datacache::DataMap &data) {
	stage1_scorefxn_ = (tag->hasOption("stage1_scorefunction")) ? (data.get<core::scoring::ScoreFunction *>("scorefxns",
		tag->getOption<std::string>("stage1_scorefunction")))->clone()
		: core::scoring::ScoreFunctionFactory::create_score_function("score3");
	stage2_scorefxn_ = (tag->hasOption("stage2_scorefunction")) ? (data.get<core::scoring::ScoreFunction *>("scorefxns",
		tag->getOption<std::string>("stage2_scorefunction")))->clone()
		: core::scoring::ScoreFunctionFactory::create_score_function("score4_smooth_cart");
	stage1_moves_ = tag->getOption<core::Size>("stage1_moves", 10000);
	stage1_move_restore_ = tag->getOption<core::Size>("stage1_move_restoring_segs", 0);
	stage1_multi_freq_ = tag->getOption<core::Real>("stage1_multi_sse_freq", 0.5);
	stage1_frag_freq_ = tag->getOption<core::Real>("stage1_frag_freq", 0.5);
	stage1_twist_freq_ = tag->getOption<core::Real>("stage1_twist_freq", 0.2);
	stage1_temp_ = tag->getOption<core::Real>("stage1_temperature", 1.0);
	rot_ = tag->getOption<core::Real>("rotation_stdev", 10.0);
	transl_ = tag->getOption<core::Real>("translation_stdev", 1.0);

	std::string template_file = tag->getOption<std::string>("template_pose", "NONE");
	template_ = (template_file != "NONE") ? core::import_pose::pose_from_file(template_file) : core::pose::PoseOP(nullptr);
	if ( template_ ) {
		for ( auto &frags : {frags3_, frags9_} ) {
			for ( auto frame = frags->nonconst_begin(); frame != frags->nonconst_end(); ++frame ) {
				frame->steal(*template_);
			}
		}
	}


	stage1_minimization_ = tag->getOption<bool>("stage1_minimization", false);


	if ( tag->hasOption("modify_segments") ) {
		utility::vector1<std::string> segments = utility::string_split(tag->getOption<std::string>("modify_segments"), ',');
		for ( auto segment: segments ) {
			if ( segment.size() != 3 ) {
				TR.Error << "Format of each segment must be segIndex-firstRes-lastRes " << std::endl;
			}
			utility::vector1<std::string> index_start_end = utility::string_split(segment, '-');
			core::Size i(0), s(0), e(0);
			std::istringstream ss1(index_start_end[1]);
			std::istringstream ss2(index_start_end[2]);
			std::istringstream ss3(index_start_end[3]);
			ss1 >> i;
			ss2 >> s;
			ss3 >> e;
			seg_.push_back(i);
			newstart_.push_back(s);
			newend_.push_back(e);
		}
		std::sort(seg_.begin(), seg_.end());
		std::sort(newstart_.begin(), newstart_.end());
		std::sort(newend_.begin(), newend_.end());
	}

	if ( tag->hasOption("add_segments") ) {
		utility::vector1<std::string> add_segments = utility::string_split(tag->getOption<std::string>("add_segments"), ',');
		for ( auto add_segment: add_segments ) {
			if ( add_segment.size() != 2 ) {
				TR.Error << "Format of each added segment must be firstRes-lastRes " << std::endl;
			}
			utility::vector1<std::string> start_end = utility::string_split(add_segment, '-');
			core::Size j(0), k(0);
			std::istringstream ss1(start_end[1]);
			std::istringstream ss2(start_end[2]);
			ss1 >> j;
			ss2 >> k;
			addstart_.push_back(j);
			addend_.push_back(k);
		}
		std::sort(addstart_.begin(), addstart_.end());
		std::sort(addend_.begin(), addend_.end());
	}

	stage2_moves_ = tag->getOption<core::Size>("stage2_moves", 1000);
	stage2_segment_freq_ = tag->getOption<core::Real>("stage2_segment_freq", 0.25);
	stage2_targgaps_freq_ = tag->getOption<core::Real>("stage2_targgaps_freq", 0.25);
	stage2_temp_ = tag->getOption<core::Real>("stage2_temperature", 1.0);
	stage2_models_ = tag->getOption<core::Size>("stage2_models", 1);

	if ( tag->hasOption("frags3") ) {
		frags3_ = core::fragment::FragmentIO().read_data(tag->getOption<std::string>("frags3"));
	}
	if ( tag->hasOption("frags9") ) {
		frags9_ = core::fragment::FragmentIO().read_data(tag->getOption<std::string>("frags9"));
	}
	rbmover_ = rbsegment_relax::GaussianRBSegmentMover(rot_, transl_);
	helixmover_ = rbsegment_relax::HelicalGaussianMover(rot_, transl_, 0, 0);

	std::string s1_res = tag->getOption<std::string>("stage1_residues", "AUTO");
	if ( s1_res != "AUTO" ) {
		s1_res_ = core::select::residue_selector::get_residue_selector(s1_res, data);
	}

	std::string s2_res = tag->getOption<std::string>("stage2_residues", "AUTO");
	if ( s2_res != "AUTO" ) {
		s2_res_ = core::select::residue_selector::get_residue_selector(s2_res, data);
	}

	std::string dihed2_res = tag->getOption<std::string>("stage2_residues_no_dihedral_csts", "NONE");
	if ( dihed2_res != "NONE" ) {
		dihed2_res_ = core::select::residue_selector::get_residue_selector(dihed2_res, data);
	}

	if ( tag->hasOption("stage1_rigid_residues") ) {
		std::string s1_rigid_bodies = tag->getOption<std::string>("stage1_rigid_residues", "NONE");
		if ( s1_rigid_bodies != "NONE" ) {
			utility::vector1<std::string> selectors;
			boost::algorithm::split(selectors, s1_rigid_bodies, boost::is_any_of(", "));
			for ( auto const &selector : selectors ) {
				s1_rigid_res_.push_back(core::select::residue_selector::get_residue_selector(selector, data));
			}
		}
	}
}


void ConfChangeMover::provide_xml_schema(utility::tag::XMLSchemaDefinition &xsd) {
	utility::tag::AttributeList top_attlist;

	top_attlist + utility::tag::XMLSchemaAttribute::attribute_w_default("template_pose", utility::tag::xs_string, "Template pose for restraints generation", "NONE") +
		utility::tag::XMLSchemaAttribute("stage1_scorefunction", utility::tag::xs_string, "Stage1 scorefunction") +
		utility::tag::XMLSchemaAttribute("stage2_scorefunction", utility::tag::xs_string, "Stage 2 scorefunction") +
		utility::tag::XMLSchemaAttribute::attribute_w_default("stage1_moves", utility::tag::xsct_non_negative_integer, "Number of  moves to apply in stage 1",
		"10000") +
		utility::tag::XMLSchemaAttribute::attribute_w_default("stage1_move_restoring_segs", utility::tag::xsct_non_negative_integer,
		"At this stage1 move, the initial segments corresponding to DSSP SSEs will be restored (i.e. they can move independently)",
		"0")

		+ utility::tag::XMLSchemaAttribute("modify_segments", utility::tag::xs_string, "comma separated list of segments to modify in the format segIndex-firstRes-lastRes")

		+ utility::tag::XMLSchemaAttribute("add_segments", utility::tag::xs_string, "comma separated list of segments to add in the format firstRes-lastRes") +
		utility::tag::XMLSchemaAttribute::attribute_w_default("stage1_rigid_residues", utility::tag::xs_string,
		"comma separated list of residue selectors that should collectively be treated as rigid bodies during stage 1", "NONE")

		+ utility::tag::XMLSchemaAttribute::attribute_w_default("stage1_multi_sse_freq", utility::tag::xsct_real,
		"Frequency from 0 to 1 of simultaneous multi-SSE moves", "0.5") +
		utility::tag::XMLSchemaAttribute::attribute_w_default("stage1_frag_freq", utility::tag::xsct_real,
		"Frequency from 0 to 1 of fragment insertion moves in stage1",
		"0.5") +
		utility::tag::XMLSchemaAttribute::attribute_w_default("stage1_twist_freq", utility::tag::xsct_real, "Frequency from 0 to 1 of helical moves", "0.2") +
		utility::tag::XMLSchemaAttribute::attribute_w_default("stage1_temperature", utility::tag::xsct_real, "Temperature used for Monte Carlo Metropolis criterion", "1.0")

		+ utility::tag::XMLSchemaAttribute::attribute_w_default("stage1_residues", utility::tag::xs_string, "Residue selector to apply in stage1 (AUTO=detect automatically)", "AUTO")

		+ utility::tag::XMLSchemaAttribute::attribute_w_default("stage2_residues", utility::tag::xs_string, "Residue selector to apply in stage2 (AUTO=detect automatically)", "AUTO") +
		utility::tag::XMLSchemaAttribute::attribute_w_default("stage2_residues_no_dihedral_csts", utility::tag::xs_string,
		"Residue selector on which dihedral constraints won't be applied in stage2", "NONE")

		+ utility::tag::XMLSchemaAttribute::attribute_w_default("rotation_stdev", utility::tag::xsct_real, "Rotation angle (standard deviation) applied to SSEs", "10") +
		utility::tag::XMLSchemaAttribute::attribute_w_default("translation_stdev", utility::tag::xsct_real, "Translation distance (standard deviation) applied to SSEs", "1.0")

		+ utility::tag::XMLSchemaAttribute("stage1_minimization", utility::tag::xsct_rosetta_bool, "Minimize during stage 1, default false")
		+
		utility::tag::XMLSchemaAttribute::attribute_w_default("stage2_moves", utility::tag::xsct_non_negative_integer, "Number of moves to apply in stage 2", "1000") +
		utility::tag::XMLSchemaAttribute::attribute_w_default("stage2_segment_freq", utility::tag::xsct_real,
		"Frequency from 0 to 1 of segment insertion moves with segments derived from input conformation",

		"0.25") + utility::tag::XMLSchemaAttribute::attribute_w_default("stage2_targgaps_freq", utility::tag::xsct_real,
		"Frequency from 0 to 1 of targeting chain-breaks, higher freq means less diversification of loop regions",
		"0.5") +
		utility::tag::XMLSchemaAttribute::attribute_w_default("stage2_temperature", utility::tag::xsct_real, "Temperature used for Monte Carlo Metropolis criterion",
		"1.0") +
		utility::tag::XMLSchemaAttribute::attribute_w_default("stage2_models", utility::tag::xsct_non_negative_integer, "stage2 number of models that will be generated, models will different in stage2 sampling", "1")

		+ utility::tag::XMLSchemaAttribute("frags3", utility::tag::xs_string, "3-mer fragment file for fragments insertion") +
		utility::tag::XMLSchemaAttribute("frags9", utility::tag::xs_string, "9-mer fragment file for fragments insertion");
	std::string description = "The mover sample new protein conformations in two steps: 1) Rigid-body perturbation of secondary structure elements (SSEs) with optional dihedrals change through fragments insertion, 2) Fragment-based loops closure and modeling";
	moves::xsd_type_definition_w_attributes(xsd, mover_name(), description, top_attlist);
}

std::string ConfChangeMover::get_name() const {
	return mover_name();
}

std::string ConfChangeMover::mover_name() {
	return "ConfChangeMover";
}

std::string ConfChangeMoverCreator::keyname() const {
	return ConfChangeMover::mover_name();
}

moves::MoverOP ConfChangeMoverCreator::create_mover() const {
	return moves::MoverOP(new ConfChangeMover);
}

void ConfChangeMoverCreator::provide_xml_schema(utility::tag::XMLSchemaDefinition &xsd) const {
	ConfChangeMover::provide_xml_schema(xsd);
}

}
}
