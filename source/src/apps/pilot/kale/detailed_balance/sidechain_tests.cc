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
#include <core/conformation/Residue.hh>
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
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/loop_modeling/LoopMover.hh>

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
using core::scoring::get_score_function;

using protocols::moves::Mover;
using protocols::moves::MoverOP;
using protocols::moves::MonteCarlo;
using protocols::moves::MonteCarloOP;
using protocols::simple_moves::sidechain_moves::SidechainMover;
using protocols::simple_moves::sidechain_moves::SidechainMoverOP;
using protocols::loops::Loop;

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.kale.monte_carlo" );

// Options {{{1
OPT_2GRP_KEY(File, kale, in, pdb)
OPT_2GRP_KEY(IntegerVector, kale, in, residues)
OPT_2GRP_KEY(Integer, kale, in, iterations)
OPT_2GRP_KEY(Boolean, kale, out, quiet)
OPT_2GRP_KEY(Integer, kale, out, trajectory)
// }}}1

class SamplingManager { // {{{1

	public:

		SamplingManager() {
			using core::pack::task::TaskFactory;
			using core::pack::task::TaskFactoryOP;
			using core::pack::task::operation::RestrictToRepacking;
			using core::pack::task::operation::PreventRepacking;
			using core::pack::task::operation::PreventRepackingOP;

			pdb_path = option[OptionKeys::kale::in::pdb]();
			first_residue = option[OptionKeys::kale::in::residues]()[1];
			last_residue = option[OptionKeys::kale::in::residues]()[2];
			iterations = option[OptionKeys::kale::in::iterations]();

			core::import_pose::pose_from_file(pose, pdb_path, core::import_pose::PDB_file);

			TaskFactoryOP task_factory = new TaskFactory;
			PreventRepackingOP prevent_repacking = new PreventRepacking();

			for (Size i = 1; i <= pose.total_residue(); i++) {
				if (i < first_residue or i > last_residue) {
					prevent_repacking->include_residue(i);
				}
			}

			task_factory->push_back(new RestrictToRepacking);
			task_factory->push_back(prevent_repacking);

			sidechain_mover = new SidechainMover;
			sidechain_mover->set_preserve_detailed_balance(true);
			sidechain_mover->set_task_factory(task_factory);

			score_function = new ScoreFunction;
			monte_carlo = new MonteCarlo(pose, *score_function, 1);
		}

		void update() {
			sidechain_mover->apply(pose);
			Real proposal_ratio = sidechain_mover->last_proposal_density_ratio();
			monte_carlo->boltzmann(pose, "sidechain", proposal_ratio);
		}

	public:
		Pose pose;
		SidechainMoverOP sidechain_mover;
		ScoreFunctionOP score_function;
		MonteCarloOP monte_carlo;

		string pdb_path;
		Size first_residue, last_residue;
		Size iterations;

};

class OutputManager { // {{{1

	public:

		OutputManager(SamplingManager &master) : sampler(master) { // {{{2
			log_coordinates.open("coordinates.dat");
			log_pivots.open("pivots.dat");

			quiet = option[OptionKeys::kale::out::quiet]();
			dump_pose_freq = option[OptionKeys::kale::out::trajectory]();
		}

		~OutputManager() { // {{{2
			log_coordinates.close();
		}
		// }}}2

		void write_header() { // {{{2
			cout << "Peptide:        " << sampler.pdb_path << endl;
			cout << "Iterations:     " << sampler.iterations << endl;
			cout << "Random Seed:    " << numeric::random::rg().get_seed() << endl;

			log_pivots << sampler.first_residue << " ";
			log_pivots << 0 << " ";
			log_pivots << sampler.last_residue << endl;
		}

		void write_progress(Size iterations_so_far) { // {{{2
			if (quiet == true) return;
			cerr << "\r[" << iterations_so_far << "/" << sampler.iterations << "]";
		}

		void write_footer() { // {{{2
			cerr << endl;
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
			log_coordinates << endl << "iteration  " << iteration;
			log_coordinates << endl << "sequence   " << sampler.pose.sequence();

			for (int i = sampler.first_residue; i <= sampler.last_residue; i++) {
				log_coordinates << endl << "chi        ";
				Size num_chi = sampler.pose.residue(i).nchi();

				for (int j = 1; j <= num_chi; j++) {
					log_coordinates << setw(10) << sampler.pose.chi(j, i) << " ";
				}
			}

			log_coordinates << endl;
		}
		// }}}2

	private:
		SamplingManager & sampler;
		ofstream log_coordinates, log_pivots;
		int dump_pose_freq;
		bool quiet;
};
// }}}1

// Application {{{1
int main(int argc, char* argv []) {
	NEW_OPT(kale::in::pdb, "Input PDB file", "");
	NEW_OPT(kale::in::residues, "Residues to focus on.", 0);
	NEW_OPT(kale::in::iterations, "How many iterations should be used", 1000);
	NEW_OPT(kale::out::quiet, "Suppress progress updates", false);
	NEW_OPT(kale::out::trajectory, "How often poses should be dumped", 0);

	devel::init(argc, argv);

	SamplingManager sampler;
	OutputManager output(sampler);

	output.write_header();
	output.dump_pose(0, true);
	output.dump_statistics(0);

	for (int i = 1; i <= sampler.iterations; i++) {
		sampler.update();
		output.write_progress(i);
		output.dump_pose(i);
		output.dump_statistics(i);
	}

	output.write_footer();
	output.dump_pose(sampler.iterations, true);
}
// }}}1

// vim: foldmethod=marker

