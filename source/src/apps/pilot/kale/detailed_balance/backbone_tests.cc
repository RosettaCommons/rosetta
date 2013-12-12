// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Headers {{{1

// Core headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

// Protocol headers
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/kinematic_closure/ClosureProblem.hh>
#include <protocols/kinematic_closure/ClosureSolution.hh>
#include <protocols/kinematic_closure/samplers/BalancedKicSampler.hh>
#include <protocols/kinematic_closure/perturbers/Perturber.hh>
#include <protocols/kinematic_closure/perturbers/RamaPerturber.hh>
#include <protocols/kinematic_closure/perturbers/UniformPerturber.hh>
#include <protocols/kinematic_closure/perturbers/VicinityPerturber.hh>
#include <protocols/kinematic_closure/perturbers/WalkingPerturber.hh>
#include <protocols/kinematic_closure/pivot_pickers/PivotPicker.hh>
#include <protocols/kinematic_closure/pivot_pickers/FixedPivots.hh>
#include <protocols/kinematic_closure/solution_pickers/FilteredSolutions.hh>
#include <protocols/kinematic_closure/solution_pickers/RandomSolutions.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>

// Utility headers
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/backrub.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>
#include <numeric/random/random.hh>
#include <boost/foreach.hpp>

#define foreach BOOST_FOREACH

// C++ headers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

// Namespaces {{{1
using namespace std;
using namespace basic::options;

using core::Size;
using core::Real;
using core::pose::Pose;
using core::scoring::ScoreFunction;
using core::scoring::ScoreFunctionOP;
using core::scoring::getScoreFunction;

using protocols::moves::Mover;
using protocols::moves::MoverOP;
using protocols::moves::MonteCarlo;
using protocols::moves::MonteCarloOP;
using protocols::loops::Loop;
using protocols::kinematic_closure::ClosureProblem;
using protocols::kinematic_closure::ClosureProblemOP;
using protocols::kinematic_closure::ClosureSolution;
using protocols::kinematic_closure::ClosureSolutionCOP;
using protocols::kinematic_closure::SolutionList;
using protocols::kinematic_closure::samplers::BalancedKicSampler;
using protocols::kinematic_closure::perturbers::PerturberOP;
using protocols::kinematic_closure::perturbers::RamaPerturber;
using protocols::kinematic_closure::perturbers::UniformPerturber;
using protocols::kinematic_closure::perturbers::VicinityPerturber;
using protocols::kinematic_closure::perturbers::WalkingPerturber;
using protocols::kinematic_closure::pivot_pickers::PivotPickerOP;
using protocols::kinematic_closure::pivot_pickers::FixedPivots;
using protocols::kinematic_closure::solution_pickers::FilteredSolutions;
using protocols::kinematic_closure::solution_pickers::RandomSolutions;
using protocols::simple_moves::sidechain_moves::SidechainMover;
using protocols::simple_moves::sidechain_moves::SidechainMoverOP;

basic::Tracer TR("apps.pilot.kale.monte_carlo");

// Options {{{1
OPT_2GRP_KEY(File, kale, in, pdb)
OPT_2GRP_KEY(File, kale, out, pdb)
OPT_2GRP_KEY(Boolean, kale, out, quiet)
OPT_2GRP_KEY(Integer, kale, mc, iterations)
OPT_2GRP_KEY(Integer, kale, mc, temperature)
OPT_2GRP_KEY(IntegerVector, kale, kic, pivots)
OPT_2GRP_KEY(String, kale, kic, closure_move)
OPT_2GRP_KEY(String, kale, kic, sampling_move)
OPT_2GRP_KEY(String, kale, kic, breadth_move)
OPT_2GRP_KEY(RealVector, kale, kic, move_weights)
OPT_2GRP_KEY(String, kale, kic, score_function)
OPT_2GRP_KEY(Integer, kale, kic, dump_pose)
OPT_2GRP_KEY(Integer, kale, kic, dump_stats)

// Typedefs {{{1
class ClosureMover;
class BreadthMover;

typedef utility::pointer::owning_ptr<ClosureMover> ClosureMoverOP;
typedef utility::pointer::owning_ptr<ClosureMover const> ClosureMoverCOP;

typedef utility::pointer::owning_ptr<BreadthMover> BreadthMoverOP;
typedef utility::pointer::owning_ptr<BreadthMover const> BreadthMoverCOP;
// }}}1

class ClosureMover : public Mover { // {{{1

public:

	ClosureMover(
			Pose const & pose,
			Size first, Size last, Size cut,
			string closure_move, string sampling_move) {

		loop_ = Loop(first, last, cut);
		pivot_picker_ = new FixedPivots(first, last, cut);
		closure_move_ = closure_move;

		if (sampling_move == "vicinity") perturber_ = new VicinityPerturber(pose);
		else if (sampling_move == "rama") perturber_ = new RamaPerturber;
		else if (sampling_move == "uniform") perturber_ = new UniformPerturber;
		else if (sampling_move == "walking") perturber_ = new WalkingPerturber;
		else utility_exit_with_message("Unknown sampling move: " + sampling_move);
	}

	void apply(Pose & pose) {
		if (closure_move_ == "naive") naive_apply(pose);
		else if (closure_move_ == "balanced") balanced_apply(pose);
		else utility_exit_with_message("Unknown closure move: " + closure_move_);
	}

	string get_name() const { return "KicSampler"; }

protected:

	void naive_apply(Pose & pose) {
		ClosureProblemOP problem = new ClosureProblem;
		SolutionList solutions;

		// Setup and solve the closure problem.
		problem->frame(pose, loop_, pivot_picker_);
		perturber_->perturb(pose, problem);
		problem->solve(solutions);

		// Randomly pick a solution to apply.
		if (not solutions.empty()) {
			Size index = numeric::random::random_range(1, solutions.size());
			solutions[index]->apply(pose);
		}
	}

	void balanced_apply(Pose & pose) {
		ClosureProblemOP problem = new ClosureProblem;
		ClosureSolutionCOP solution;
		SolutionList unperturbed_solutions, perturbed_solutions;

		// Solve the unperturbed problem.
		problem->frame(pose, loop_, pivot_picker_);
		problem->solve(unperturbed_solutions);

		// Solve the perturbed problem.
		perturber_->perturb_with_balance(pose, problem);
		problem->solve(perturbed_solutions);

		// Pick a solution to apply.
		solution = BalancedKicSampler::pick_solution(
				unperturbed_solutions, perturbed_solutions);
		solution->apply(pose);
	}

	void uniform_perturb();

	void rama_perturb();

	void biased_rama_perturb();

private:
	Loop loop_;
	PivotPickerOP pivot_picker_;
	PerturberOP perturber_;
	string closure_move_;
	bool check_rama_;
};

class BreadthMover : public Mover { // {{{1

public:

	BreadthMover(
			Pose const & pose,
			Size first, Size last, Size cut,
			string breadth_move, string sampling_move) {

		first_ = first;
		last_ = last;
		breadth_move_ = breadth_move;
		sampling_move_ = sampling_move;
	}

	void apply(Pose & pose) {
		if (breadth_move_ == "uniform") uniform_apply(pose);
		else if (breadth_move_ == "omega") omega_apply(pose);
		else if (breadth_move_ == "rama") rama_apply(pose);
		else utility_exit_with_message("Unknown breadth move: " + breadth_move_);
	}

	string get_name() const { return "BreadthMove"; }

private:

	void uniform_apply(Pose & pose) {
		Size index = numeric::random::random_range(first_, last_);
		Size which_angle = numeric::random::random_range(1, 3);
		Real angle_value = 360 * numeric::random::uniform();

		switch (which_angle) {
			case 1: pose.set_phi(index, angle_value); break;
			case 2: pose.set_psi(index, angle_value); break;
			case 3: pose.set_omega(index, angle_value); break;
		}
	}

	void omega_apply(Pose & pose) {
		Size index = numeric::random::random_range(first_, last_);
		Real angle = 360 * numeric::random::uniform();
		pose.set_omega(index, angle);
	}

	void rama_apply(Pose & pose) {
		using core::scoring::Ramachandran;
		using core::scoring::ScoringManager;

		Real phi, psi;
		Ramachandran const & rama =
			ScoringManager::get_instance()->get_Ramachandran();

		Size index = numeric::random::random_range(first_, last_);
		rama.uniform_phipsi_from_allowed_rama(pose.aa(index), phi, psi);

		pose.set_phi(index, phi);
		pose.set_psi(index, psi);
	}

private:
	Size first_, last_;
	string breadth_move_, sampling_move_;

};
// }}}1

class SamplingManager { // {{{1

	friend class OutputManager;

	public:

		SamplingManager() { // {{{2
			if (option[OptionKeys::kale::in::pdb].active() == false) {
				utility_exit_with_message("No input PDB file specified.");
			}

			if (option[OptionKeys::kale::kic::closure_move].active() == false) {
				utility_exit_with_message("No closure move specified.");
			}

			if (option[OptionKeys::kale::kic::breadth_move].active() == false) {
				utility_exit_with_message("No breadth move specified.");
			}

			closure_move = option[OptionKeys::kale::kic::closure_move]();
			sampling_move = option[OptionKeys::kale::kic::sampling_move]();
			breadth_move = option[OptionKeys::kale::kic::breadth_move]();
			weights = option[OptionKeys::kale::kic::move_weights]();
			score_name = option[OptionKeys::kale::kic::score_function]();
			pdb_path = option[OptionKeys::kale::in::pdb]();
			first_index = option[OptionKeys::kale::kic::pivots]()[1];
			cut_index = option[OptionKeys::kale::kic::pivots]()[2];
			last_index = option[OptionKeys::kale::kic::pivots]()[3];
			temperature = option[OptionKeys::kale::mc::temperature]();
			iterations = option[OptionKeys::kale::mc::iterations]();

			if (temperature < 0) {
				temperature = numeric_limits<Real>::infinity();
			}

			if (weights[1] == 0) closure_move = "off";
			if (weights[2] == 0) breadth_move = "off";

			core::import_pose::pose_from_pdb(pose, pdb_path);

			score_function = new ScoreFunction;

			if (score_name == "rama") {
				score_function->set_weight(core::scoring::rama, 1);
			}
			else if (score_name != "off") {
				utility_exit_with_message("Unknown score function: " + score_name);
			}

			monte_carlo = new MonteCarlo(pose, *score_function, temperature);
		}

		void setup() { // {{{2
			Size first = first_index;
			Size cut = cut_index;
			Size last = last_index;

			closure_mover = new ClosureMover(
					pose, first, last, cut, closure_move, sampling_move);
			breadth_mover = new BreadthMover(
					pose, first, last, cut, breadth_move, sampling_move);

			sidechain_mover = new SidechainMover;
			sidechain_mover->set_preserve_detailed_balance(true);

			using core::pack::task::TaskFactory;
			using core::pack::task::TaskFactoryOP;
			using core::pack::task::operation::InitializeFromCommandline;
			using core::pack::task::operation::IncludeCurrent;

			TaskFactoryOP task_factory = new TaskFactory;
			task_factory->push_back(new InitializeFromCommandline);
			task_factory->push_back(new IncludeCurrent);

			sidechain_mover->set_task_factory(task_factory);
		}

		void update() { // {{{2
			Size total_weight = weights[1] + weights[2] + weights[3];
			Real closure_threshold = weights[1] / total_weight;
			Real breadth_threshold = weights[2] / total_weight + closure_threshold;
			Real move_choice = numeric::random::uniform();
			Real proposal_ratio = 1;
			string move_name;

			cout << "total_weight:      " << total_weight << endl;
			cout << "closure_threshold: " << closure_threshold << endl;
			cout << "breadth_threshold: " << breadth_threshold << endl;
			cout << "move_choice:       " << move_choice << endl;

			if (total_weight == 0) {
				utility_exit_with_message("No active movers.");
			}

			if (move_choice < closure_threshold) {
				closure_mover->apply(pose);
				move_name = "closure";
				cout << "  Picking closure." << endl;
			}
			else if (move_choice < breadth_threshold) {
				breadth_mover->apply(pose);
				move_name = "breadth";
				cout << "  Picking breadth." << endl;
			}
			else {
				sidechain_mover->apply(pose);
				proposal_ratio *= sidechain_mover->last_proposal_density_ratio();
				move_name = "sidechain";
				cout << "  Picking sidechain." << endl;
			}
			cout << endl;

			monte_carlo->boltzmann(pose, move_name, proposal_ratio);
		}

		int get_iterations() { // {{{2
			return iterations;
		}
		// }}}2

	private:

		// Sampling helpers {{{2
		Pose pose;
		ClosureMoverOP closure_mover;
		BreadthMoverOP breadth_mover;
		SidechainMoverOP sidechain_mover;
		ScoreFunctionOP score_function;
		MonteCarloOP monte_carlo;

		// Command line options {{{2
		string pdb_path;
		string closure_move;
		string sampling_move;
		string breadth_move;
		string score_name;
		Size first_index, cut_index, last_index;
		utility::vector1<Real> weights;
		Size iterations;
		Real temperature;
		// }}}2

};

class OutputManager { // {{{1

	public:

		OutputManager(SamplingManager &master) : sampler(master) { // {{{2
			log_coordinates.open("coordinates.dat");
			log_solutions.open("solutions.dat");
			log_pivots.open("pivots.dat");

			quiet = option[OptionKeys::kale::out::quiet]();
			dump_pose_freq = option[OptionKeys::kale::kic::dump_pose]();
			dump_stats_freq = option[OptionKeys::kale::kic::dump_stats]();
		}

		~OutputManager() { // {{{2
			log_coordinates.close();
			log_solutions.close();
			log_pivots.close();
		}
		// }}}2

		void write_header() { // {{{2
			cout << "Closure Move:   " << sampler.closure_move << endl;
			cout << "Sampling Move:  " << sampler.sampling_move << endl;
			cout << "Breadth Move:   " << sampler.breadth_move << endl;
			cout << "Move Weights:   " << sampler.weights[1] << "/"
			                           << sampler.weights[2] << "/"
			                           << sampler.weights[3] << endl;
			cout << "Score Function: " << sampler.score_name << endl;
			cout << "Peptide:        " << sampler.pdb_path << endl;
			cout << "Pivots:         " << sampler.first_index << "/"
			                           << sampler.cut_index << "/"
			                           << sampler.last_index << endl;
			cout << "Random Seed:    " << numeric::random::RG.get_seed() << endl;
			cout << "Iterations:     " << sampler.iterations << endl;
			cout << "Temperature:    " << sampler.temperature << endl;

			log_pivots << sampler.first_index << " ";
			log_pivots << sampler.cut_index << " ";
			log_pivots << sampler.last_index << endl;
		}

		void write_progress(Size iterations_so_far) { // {{{2
			if (quiet == true) return;
			cerr << "\r[" << iterations_so_far << "/" << sampler.iterations << "]";
		}

		void write_footer() { // {{{2
			cerr << endl;
			sampler.monte_carlo->show_counters();
		}
		// }}}2

		void dump_pose(int iteration, bool force=false) { // {{{2
			if (force == false && dump_pose_freq <= 0) return;
			if (force == false && iteration % dump_pose_freq != 0) return;

			static int total_dumps = sampler.iterations + 1;
			static int total_digits = (int) ceil(log10(total_dumps));

			stringstream stream;
			stream << setw(total_digits) << setfill('0') << iteration;
			string counter = stream.str();

			sampler.pose.dump_pdb("trajectory/" + counter + ".pdb");
		}

		void dump_statistics(int iteration) { // {{{2
			if (dump_stats_freq <= 0)
				return;

			if (iteration % dump_stats_freq != 0)
				return;

			// Log torsion information.

			log_coordinates << endl << "iteration  " << iteration;

			//if (sampler.closure_mover != NULL) {
			//	Size solution_count = sampler.closure_mover->solutions_to_last_move();
			//	log_coordinates << endl << "solutions  " << solution_count;
			//}

			log_coordinates << endl << "phi        ";
			for (int i = sampler.first_index; i <= sampler.last_index; i++) {
				log_coordinates << setw(10) << sampler.pose.phi(i) << " ";
			}

			log_coordinates << endl << "psi        ";
			for (int i = sampler.first_index; i <= sampler.last_index; i++) {
				log_coordinates << setw(10) << sampler.pose.psi(i) << " ";
			}

			log_coordinates << endl << "omega      ";
			for (int i = sampler.first_index; i <= sampler.last_index; i++) {
				log_coordinates << setw(10) << sampler.pose.omega(i) << " ";
			}

			// Log solution count information.

			log_coordinates << endl;
		}
		// }}}2

	private:
		SamplingManager &sampler;
		ofstream log_coordinates, log_solutions, log_pivots;

		bool quiet;
		int dump_pose_freq, dump_stats_freq;
};
// }}}1

// Application {{{1
int main(int argc, char* argv []) {
	NEW_OPT(kale::in::pdb, "Input PDB file", "");
	NEW_OPT(kale::out::pdb, "Output data file", "");
	NEW_OPT(kale::out::quiet, "Hide progress updates", false);
	NEW_OPT(kale::mc::iterations, "Monte Carlo iterations", 1000);
	NEW_OPT(kale::mc::temperature, "Monte Carlo temperature", -1);
	NEW_OPT(kale::kic::pivots, "Pivot residues", 0);
	NEW_OPT(kale::kic::closure_move, "Closure move", "");
	NEW_OPT(kale::kic::breadth_move, "Breadth move", "");
	NEW_OPT(kale::kic::sampling_move, "Sampling move", "");
	NEW_OPT(kale::kic::move_weights, "Closure/breadth/sidechain move frequencies", 0.5);
	NEW_OPT(kale::kic::score_function, "Score function", "");
	NEW_OPT(kale::kic::dump_pose, "Pose dump frequency", 0);
	NEW_OPT(kale::kic::dump_stats, "Statistics dump frequency", 1);

	devel::init(argc, argv);

	SamplingManager sampler;
	OutputManager output(sampler);

	sampler.setup();
	output.write_header();
	output.dump_pose(0, true);
	output.dump_statistics(0);

	int iterations = sampler.get_iterations();

	for (int i = 1; i <= iterations; i++) {
		sampler.update();
		output.write_progress(i);
		output.dump_pose(i);
		output.dump_statistics(i);
	}

	output.write_footer();
	output.dump_pose(iterations, true);
}
// }}}1

// vim: foldmethod=marker
