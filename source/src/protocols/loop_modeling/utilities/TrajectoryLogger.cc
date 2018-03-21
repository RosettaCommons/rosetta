// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/utilities/TrajectoryLogger.hh>
#include <protocols/loop_modeling/LoopMover.hh>

// Protocol headers
#include <protocols/loops/Loops.hh>
#include <protocols/simple_filters/BuriedUnsatHbondFilter.hh>

// Core headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <boost/format.hpp>

// C++ headers
#include <iostream>

using namespace std;
using core::scoring::ScoreFunctionCOP;

namespace protocols {
namespace loop_modeling {
namespace utilities {

TrajectoryLogger::TrajectoryLogger(string const & prefix) {
	using namespace basic::options;

	reset_timer();
	set_prefix(prefix);

	have_native_pose_ = option[OptionKeys::in::file::native].user();
	if ( have_native_pose_ ) {
		pose_from_file(native_pose_, option[OptionKeys::in::file::native](), core::import_pose::PDB_file);
	}
}

// Use a raw pointer to be compatible with `this`.
TrajectoryLogger::TrajectoryLogger(
	LoopMover const * mover,
	string const & prefix)

: TrajectoryLogger(prefix) {

	init(mover);
}

TrajectoryLogger::TrajectoryLogger(
	LoopsCOP loops,
	ScoreFunctionCOP scorefxn,
	string const & prefix)

: TrajectoryLogger(prefix) {

	init(loops, scorefxn);
}

// Use a raw pointer to be compatible with `this`.
void TrajectoryLogger::init(LoopMover const * mover) {
	init(
		mover->get_tool<LoopsCOP>(ToolboxKeys::LOOPS),
		mover->get_tool<ScoreFunctionCOP>(ToolboxKeys::SCOREFXN));
}

void TrajectoryLogger::init(LoopsCOP loops, ScoreFunctionCOP scorefxn) {
	set_loops(loops);
	set_score_function(scorefxn);
	reset_timer();
}

void TrajectoryLogger::record_move(
	basic::Tracer & tr,
	Pose & pose,
	Size const i,
	Size const j,
	Size const k,
	bool const proposed,
	bool const accepted) const {

	using namespace basic::options;

	if ( tr.Trace.visible() ) {
		tr.Trace << "Propose move:         ";
		tr.Trace << boost::format("iteration: %d,%d,%d; ") % i % j % k;
		tr.Trace << boost::format("proposed: %d; ") % proposed;
		tr.Trace << boost::format("accepted: %d; ") % accepted;
		tr.Trace << boost::format("score: %.3f REU; ") % calc_score(pose);
		if ( have_native_pose_ ) {
			tr.Trace << boost::format(u8"RMSD: %.3f Å; ") % calc_rmsd_to_native(pose);
		}
		if ( ! option[ OptionKeys::run::no_prof_info_in_silentout ] ) {
			tr.Trace << boost::format("time: %d s; ") % get_timer();
		}
		tr.Trace << endl;
	}
}

void TrajectoryLogger::record_move(
	basic::Tracer & tr,
	Pose & pose,
	Size const i,
	bool const proposed,
	bool const accepted) const {

	using namespace basic::options;

	if ( tr.Trace.visible() ) {
		tr.Trace << "Propose move:         ";
		tr.Trace << boost::format("iteration: %d; ") % i;
		tr.Trace << boost::format("proposed: %d; ") % proposed;
		tr.Trace << boost::format("accepted: %d; ") % accepted;
		tr.Trace << boost::format("score: %.3f REU; ") % calc_score(pose);
		if ( have_native_pose_ ) {
			tr.Trace << boost::format(u8"RMSD: %.3f Å; ") % calc_rmsd_to_native(pose);
		}
		if ( ! option[ OptionKeys::run::no_prof_info_in_silentout ] ) {
			tr.Trace << boost::format("time: %d s; ") % get_timer();
		}
		tr.Trace << endl;
	}
}

void TrajectoryLogger::record_new_pose(
	basic::Tracer & tr,
	Pose & pose) const {

	tr.Trace << "Recover best pose:    ";
	tr.Trace << boost::format("score: %.3f REU; ") % calc_score(pose);
	if ( have_native_pose_ ) {
		tr.Trace << boost::format(u8"RMSD: %.3f Å; ") % calc_rmsd_to_native(pose);
	}
	tr.Trace << endl;
}

void TrajectoryLogger::record_new_score_function(
	basic::Tracer & tr,
	Pose & pose) const {

	using core::scoring::fa_rep;
	using core::scoring::rama;
	using core::scoring::rama2b;
	using core::scoring::chainbreak;

	tr.Trace << "Ramp score function:  ";
	tr.Trace << boost::format("score: %.3f REU; ") % calc_score(pose);
	tr.Trace << boost::format("chainbreak: %.2f; ") % scorefxn_->get_weight(chainbreak);
	tr.Trace << boost::format("fa_rep: %.2f; ") % scorefxn_->get_weight(fa_rep);
	tr.Trace << boost::format("rama: %.2f; ") % scorefxn_->get_weight(rama);
	tr.Trace << boost::format("rama2b: %.2f; ") % scorefxn_->get_weight(rama2b);
	tr.Trace << endl;
}

void TrajectoryLogger::record_new_temperature(
	basic::Tracer & tr,
	Real temperature) const {

	tr.Trace << "Ramp temperature:     ";
	tr.Trace << boost::format("temperature: %.3f; ") % temperature;
	tr.Trace << endl;
}

void TrajectoryLogger::record_endpoint(
	basic::Tracer & tr,
	Pose & pose) const {

	using namespace basic::options;
	using core::scoring::chainbreak;
	using core::pose::setPoseExtraScore;
	using protocols::simple_filters::BuriedUnsatHbondFilter;

	// Report and save the final score.

	Real score = calc_score(pose);
	Real chainbreak_score = pose.energies().total_energies()[chainbreak];

	tr.Info << boost::format("Final Score: %.3f REU; ") % score;
	tr.Info << boost::format("Chainbreak: %.3f REU") % chainbreak_score << endl;
	setPoseExtraScore(pose, prefix_ + "score", score);
	setPoseExtraScore(pose, prefix_ + "chainbreak", chainbreak_score);

	// Report and save the final RMSD-to-native.

	if ( have_native_pose_ ) {
		Real rmsd = calc_rmsd_to_native(pose);
		tr.Info << boost::format(u8"Final RMSD: %.3f Å (all loops, backbone heavy atom, no superposition)") % rmsd << endl;
		setPoseExtraScore(pose, prefix_ + "rmsd", rmsd);
	}

	// Report and save the number of buried unsatisfied H-bonds.

	if ( have_native_pose_ and pose.is_fullatom() ) {
		protocols::simple_filters::BuriedUnsatHbondFilter unsats;
		Real delta_unsats, delta_bb_unsats, delta_sc_unsats;

		unsats.set_report_bb_heavy_atom_unsats(true);
		delta_bb_unsats = unsats.compute(pose) - unsats.compute(native_pose_);
		unsats.set_report_bb_heavy_atom_unsats(false);

		unsats.set_report_sc_heavy_atom_unsats(true);
		delta_sc_unsats = unsats.compute(pose) - unsats.compute(native_pose_);
		unsats.set_report_sc_heavy_atom_unsats(false);

		delta_unsats = delta_bb_unsats + delta_sc_unsats;

		tr.Info << boost::format(u8"Δ Buried Unsats (Total): %.3f") % delta_unsats << endl;
		tr.Info << boost::format(u8"Δ Buried Unsats (Backbone): %.3f") % delta_bb_unsats << endl;
		tr.Info << boost::format(u8"Δ Buried Unsats (Sidechain): %.3f") % delta_sc_unsats << endl;

		setPoseExtraScore(pose, prefix_ + "delta_unsats", delta_unsats);
		setPoseExtraScore(pose, prefix_ + "delta_bb_unsats", delta_bb_unsats);
		setPoseExtraScore(pose, prefix_ + "delta_sc_unsats", delta_sc_unsats);
	}

	// Report how much time passed.

	if ( ! option[ OptionKeys::run::no_prof_info_in_silentout ] ) {
		long elapsed_time = get_timer();
		tr.Info << boost::format("Elapsed Time: %d sec") % elapsed_time << endl;
		setPoseExtraScore(pose, prefix_ + "time", elapsed_time);
	}
}

Real TrajectoryLogger::calc_score(Pose & pose) const {
	return scorefxn_->score(pose);
}

Real TrajectoryLogger::calc_rmsd_to_native(Pose & pose) const {
	using core::scoring::rmsd_no_super_subset;
	using core::scoring::is_protein_backbone;

	debug_assert(have_native_pose_);

	ObjexxFCL::FArray1D<bool> loop_residues (pose.size(), false);
	for ( Size i = 1; i <= pose.size(); i++ ) {
		loop_residues[i-1] = loops_->is_loop_residue(i);
	}

	return rmsd_no_super_subset(
		native_pose_, pose, loop_residues, is_protein_backbone);
}

string TrajectoryLogger::get_prefix() const {
	return prefix_;
}

void TrajectoryLogger::set_prefix(string prefix) {
	prefix_ = prefix;
}

Pose const & TrajectoryLogger::get_native_pose() const {
	return native_pose_;
}

void TrajectoryLogger::set_native_pose(Pose const & pose) {
	native_pose_ = pose;
	have_native_pose_ = true;
}

void TrajectoryLogger::unset_native_pose() {
	have_native_pose_ = false;
}

LoopsCOP TrajectoryLogger::get_loops() const {
	return loops_;
}

void TrajectoryLogger::set_loops(LoopsCOP loops) {
	loops_ = loops;
}

ScoreFunctionCOP TrajectoryLogger::get_score_function() const {
	return scorefxn_;
}

void TrajectoryLogger::set_score_function(ScoreFunctionCOP scorefxn) {
	scorefxn_ = scorefxn;
}

long TrajectoryLogger::get_timer() const {
	return time(nullptr) - start_time_;
}

void TrajectoryLogger::reset_timer() {
	start_time_ = time(nullptr);
}


} // namespace utilities
} // namespace loop_modeling
} // namespace protocols
