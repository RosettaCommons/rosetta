// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/star/StarAbinitio.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <protocols/star/StarAbinitio.hh>

// C/C++ headers
#include <iostream>
#include <string>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <numeric/interpolate.hh>
#include <numeric/prob_util.hh>
#include <numeric/xyzVector.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/SecondaryStructure.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/ShortestPathInFoldTree.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/util/kinematics_util.hh>
#include <protocols/constraints_additional/MaxSeqSepConstraintSet.hh>
#include <protocols/comparative_modeling/ThreadingJob.hh>
#include <protocols/comparative_modeling/util.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/util.hh>
#include <protocols/medal/util.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/nonlocal/BiasedFragmentMover.hh>
#include <protocols/nonlocal/Policy.hh>
#include <protocols/nonlocal/PolicyFactory.hh>
#include <protocols/nonlocal/util.hh>
#include <protocols/simple_moves/SaneMinMover.hh>
#include <protocols/simple_moves/rational_mc/RationalMonteCarlo.hh>

// Package headers
#include <protocols/star/Extender.hh>
#include <protocols/star/util.hh>

namespace protocols {
namespace star {

static THREAD_LOCAL basic::Tracer TR( "protocols.star.StarAbinitio" );

typedef utility::vector1<double> Probabilities;

using core::Real;
using core::Size;
using core::kinematics::FoldTree;
using core::pose::Pose;
using core::scoring::ScoreFunctionOP;
using core::scoring::ScoreFunctionFactory;
using numeric::xyzVector;
using protocols::loops::Loop;
using protocols::loops::Loops;
using protocols::moves::MoverOP;
using protocols::simple_moves::rational_mc::RationalMonteCarlo;
using utility::vector1;

void compute_per_residue_probabilities(Size num_residues, Size fragment_len, const FoldTree& tree, Probabilities* probs) {
	probs->resize(num_residues, 1);
	protocols::medal::invalidate_residues_spanning_cuts(tree, fragment_len, probs);
	numeric::normalize(probs->begin(), probs->end());
	numeric::print_probabilities(*probs, TR);
}

/// @detail Regulates the application of constraints during folding based on
/// distance between residues in the fold tree. The MonteCarlo object should
/// be reset after calling this function.
void update_sequence_separation(Size distance, Pose* pose) {
	using protocols::constraints_additional::MaxSeqSepConstraintSet;
	using protocols::constraints_additional::MaxSeqSepConstraintSetOP;

	MaxSeqSepConstraintSetOP new_cst( new MaxSeqSepConstraintSet(*pose->constraint_set(), pose->fold_tree()) );
	new_cst->set_max_seq_sep(distance);
	pose->constraint_set(new_cst);
	TR << "max_seq_sep => " << distance << std::endl;
}

/// @detail Adds constraints between contacting residues in different
/// aligned chunks, as well as those specified on the command line via
/// -constraints:cst_file. Call this *after* kinematics are in place.
void setup_constraints(const Loops& aligned, Pose* pose) {
	using core::id::AtomID;
	using core::kinematics::ShortestPathInFoldTree;
	using core::scoring::constraints::AtomPairConstraint;
	using core::scoring::func::HarmonicFunc;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::scoring::constraints::add_constraints_from_cmdline_to_pose(*pose);

	Size n_csts = 0;
	for ( Size i = 1; i <= aligned.size(); ++i ) {
		const Loop& ci = aligned[i];

		for ( Size j = i + 1; j <= aligned.size(); ++j ) {
			const Loop& cj = aligned[j];

			for ( Size k = ci.start(); k <= ci.stop(); ++k ) {
				const AtomID ai(pose->conformation().residue(k).atom_index("CA"), k);
				const xyzVector<Real>& p = pose->xyz(ai);

				for ( Size l = cj.start(); l <= cj.stop(); ++l ) {
					const AtomID aj(pose->conformation().residue(l).atom_index("CA"), l);
					const xyzVector<Real>& q = pose->xyz(aj);

					Real distance = p.distance(q);
					if ( distance <= option[OptionKeys::abinitio::star::initial_dist_cutoff]() ) {
						pose->add_constraint(core::scoring::constraints::ConstraintCOP( core::scoring::constraints::ConstraintOP( new AtomPairConstraint(ai, aj, core::scoring::func::FuncOP(new HarmonicFunc(distance, 2))) ) ));
						TR << "AtomPair CA " << k << " CA " << l << " HARMONIC " << distance << " 2" << std::endl;
						++n_csts;
					}
				}
			}
		}
	}

	if ( !n_csts ) {
		TR.Error << "Failed to define constraints between chunks in the non-local pairing" << std::endl;
		TR.Error << "Check -abinitio:star:initial_dist_cutoff" << std::endl;
		utility_exit();
	}

	Size max_dist_ft = ShortestPathInFoldTree(pose->fold_tree()).max_dist();
	update_sequence_separation(max_dist_ft, pose);
}

void tear_down_constraints(Pose* pose) {
	pose->remove_constraints();
}

/// @detail Instantiates a score function by name and sets constraint terms to
/// the value specified by the options system (-constraints:cst_weight)
ScoreFunctionOP setup_score(const std::string& weights, Real cb) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	ScoreFunctionOP score = ScoreFunctionFactory::create_score_function(weights);
	score->set_weight(core::scoring::atom_pair_constraint, option[OptionKeys::constraints::cst_weight]());
	score->set_weight(core::scoring::linear_chainbreak, cb);
	return score;
}

/// @detail Updates RationalMonteCarlo instance.
/// Calling this method twice without scoring a pose in between will trigger a
/// runtime assertion in core/pose/Pose.cc.
void configure_rmc(MoverOP mover, ScoreFunctionOP score, Size num_cycles, Real temperature, bool recover_low, RationalMonteCarlo* rmc) {
	rmc->set_mover(mover);
	rmc->set_score_function(score);
	rmc->set_num_trials(num_cycles);
	rmc->set_temperature(temperature);
	rmc->set_recover_low(recover_low);
}

void StarAbinitio::setup_kinematics(const Loops& aligned, const vector1<unsigned>& interior_cuts, Pose & pose) const {
	assert(aligned.num_loop() >= 2);
	assert(interior_cuts.size() == (aligned.num_loop() - 1));

	const Size num_residues = pose.size();
	const Size vres = num_residues + 1;

	xyzVector<double> center;
	aligned.center_of_mass(pose, &center);
	core::pose::addVirtualResAsRoot(center, pose);

	vector1<std::pair<int, int> > jumps;
	for ( Size i = 1; i <= aligned.num_loop(); ++i ) {
		jumps.push_back(std::make_pair(vres, aligned[i].midpoint()));
	}

	vector1<int> cuts(interior_cuts);
	cuts.push_back(num_residues);

	ObjexxFCL::FArray2D_int ft_jumps(2, jumps.size());
	for ( Size i = 1; i <= jumps.size(); ++i ) {
		ft_jumps(1, i) = std::min(jumps[i].first, jumps[i].second);
		ft_jumps(2, i) = std::max(jumps[i].first, jumps[i].second);
	}

	ObjexxFCL::FArray1D_int ft_cuts(cuts.size());
	for ( Size i = 1; i <= cuts.size(); ++i ) {
		ft_cuts(i) = cuts[i];
	}

	FoldTree tree(vres);
	bool status = tree.tree_from_jumps_and_cuts(vres,          // nres_in
		jumps.size(),  // num_jump_in
		ft_jumps,      // jump_point
		ft_cuts,       // cuts
		vres);         // root

	runtime_assert(status);
	pose.fold_tree(tree);
	core::pose::correctly_add_cutpoint_variants(pose);

	TR << pose.fold_tree() << std::endl;
}

void StarAbinitio::apply(Pose& pose) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using protocols::comparative_modeling::ThreadingJob;
	using protocols::nonlocal::BiasedFragmentMover;
	using protocols::nonlocal::PolicyFactory;
	using protocols::nonlocal::PolicyOP;

	ThreadingJob const * const job = protocols::nonlocal::current_job();

	to_centroid(pose);
	emit_intermediate(pose, "star_initial.out");

	Extender extender(job->alignment().clone(), pose.size());
	extender.set_secondary_structure(pred_ss_);
	extender.extend_unaligned(&pose);
	emit_intermediate(pose, "star_extended.out");

	const Loops& aligned = *(extender.aligned());
	//const Loops& unaligned = *(extender.unaligned());
	TR << "Aligned: " << aligned << std::endl;

	const Size num_residues = pose.size();
	setup_kinematics(aligned, extender.cutpoints(), pose);
	setup_constraints(aligned, &pose);

	Probabilities probs_sm, probs_lg;
	compute_per_residue_probabilities(num_residues, fragments_sm_->max_frag_length(), pose.fold_tree(), &probs_sm);
	compute_per_residue_probabilities(num_residues, fragments_lg_->max_frag_length(), pose.fold_tree(), &probs_lg);

	MoverOP fragments_lg_uni( new BiasedFragmentMover(PolicyFactory::get_policy("uniform", fragments_lg_, 25), probs_lg) );
	MoverOP fragments_sm_uni( new BiasedFragmentMover(PolicyFactory::get_policy("uniform", fragments_sm_, 200), probs_sm) );
	MoverOP fragments_sm_smo( new BiasedFragmentMover(PolicyFactory::get_policy("smooth", fragments_sm_, 200), probs_sm) );

	// Simulation parameters
	const Real mult = option[OptionKeys::abinitio::increase_cycles]();
	const Real temperature = option[OptionKeys::abinitio::temperature]();
	const Real chainbreak = option[OptionKeys::jumps::increase_chainbreak]();

	ScoreFunctionOP score_stage1  = setup_score("score0", 0.1 * chainbreak);
	ScoreFunctionOP score_stage2  = setup_score("score1", 0.2 * chainbreak);
	ScoreFunctionOP score_stage3a = setup_score("score2", 0.4 * chainbreak);
	ScoreFunctionOP score_stage3b = setup_score("score5", 0.4 * chainbreak);
	ScoreFunctionOP score_stage4  = setup_score("score3", 0.6 * chainbreak);

	// Stage 1
	TR << "Stage 1" << std::endl;
	RationalMonteCarlo rmc(fragments_lg_uni, score_stage1, static_cast<Size>(mult * 2000), temperature, true);
	rmc.apply(pose);

	configure_rmc(fragments_sm_uni, score_stage1, static_cast<Size>(mult * 2000), temperature, true, &rmc);
	rmc.apply(pose);
	emit_intermediate(pose, "star_stage1.out");

	// Stage 2
	TR << "Stage 2" << std::endl;
	configure_rmc(fragments_lg_uni, score_stage2, static_cast<Size>(mult * 4000), temperature, true, &rmc);
	rmc.enable_autotemp(temperature);
	rmc.apply(pose);

	configure_rmc(fragments_sm_uni, score_stage2, static_cast<Size>(mult * 4000), temperature, true, &rmc);
	rmc.apply(pose);
	emit_intermediate(pose, "star_stage2.out");

	// Stage 3
	TR << "Stage 3" << std::endl;
	for ( Size i = 1; i <= 10; ++i ) {
		ScoreFunctionOP score = ((i % 2) == 0 && i <= 7) ? score_stage3a : score_stage3b;
		configure_rmc(fragments_lg_uni, score, static_cast<Size>(mult * 4000), temperature, true, &rmc);
		rmc.apply(pose);

		configure_rmc(fragments_sm_uni, score, static_cast<Size>(mult * 4000), temperature, true, &rmc);
		rmc.apply(pose);
	}
	emit_intermediate(pose, "star_stage3.out");

	// Stage 4
	TR << "Stage 4" << std::endl;
	configure_rmc(fragments_sm_smo, score_stage4, static_cast<Size>(mult * 8000), temperature, true, &rmc);
	for ( Size i = 1; i <= 3; ++i ) {
		rmc.apply(pose);
	}
	emit_intermediate(pose, "star_stage4.out");
}

void StarAbinitio::tear_down_kinematics(Pose & pose) const {
	pose.conformation().delete_residue_slow(pose.size());
	core::util::remove_cutpoint_variants(pose);
	simple_fold_tree(pose);
}

StarAbinitio::StarAbinitio() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using core::fragment::FragmentIO;
	using core::fragment::SecondaryStructure;
	using core::optimization::MinimizerOptions;
	using core::optimization::MinimizerOptionsOP;
	using core::kinematics::MoveMap;
	using core::kinematics::MoveMapOP;
	using protocols::simple_moves::SaneMinMover;

	FragmentIO io;
	fragments_lg_ = io.read_data(option[in::file::frag9]());
	fragments_sm_ = io.read_data(option[in::file::frag3]());

	// Approximate secondary structure from fragments when psipred isn't available
	if ( option[in::file::psipred_ss2].user() ) {
		pred_ss_ = core::fragment::SecondaryStructureOP( new SecondaryStructure() );
		pred_ss_->read_psipred_ss2(option[in::file::psipred_ss2]());
	} else {
		pred_ss_ = core::fragment::SecondaryStructureOP( new SecondaryStructure(*fragments_sm_) );
	}

	// Configure the minimizer
	ScoreFunctionOP min_score = ScoreFunctionFactory::create_score_function("score4_smooth_cart");
	min_score->set_weight(core::scoring::atom_pair_constraint, option[OptionKeys::constraints::cst_weight]());

	MinimizerOptionsOP min_options;
	min_options = MinimizerOptionsOP( new MinimizerOptions("lbfgs_armijo_nonmonotone", 0.01, true, false, false) );
	min_options->max_iter(100);

	MoveMapOP mm( new MoveMap() );
	mm->set_bb(true);
	mm->set_chi(true);
	mm->set_jump(true);

	minimizer_ = protocols::simple_moves::SaneMinMoverOP( new SaneMinMover(mm, min_score, min_options, true) );
}

std::string StarAbinitio::get_name() const {
	return "StarAbinitio";
}

protocols::moves::MoverOP StarAbinitio::clone() const {
	return protocols::moves::MoverOP( new StarAbinitio(*this) );
}

protocols::moves::MoverOP StarAbinitio::fresh_instance() const {
	return protocols::moves::MoverOP( new StarAbinitio() );
}

}  // namespace star
}  // namespace protocols
