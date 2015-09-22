// Headers {{{1

#include <devel/cycles/SetupMover.hh>
#include <devel/cycles/RandomReindexingMover.hh>
#include <devel/balanced_kic/KinematicMover.hh>
#include <devel/balanced_kic/KinematicPerturber.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>
#include <numeric/random/random.hh>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/backrub.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

using namespace std;
using namespace core;
using namespace basic::options;
using namespace protocols::moves;
using namespace devel::balanced_kic;

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.kale.monte_carlo" );

// Options {{{1
OPT_2GRP_KEY(File, kale, in, pdb)
OPT_2GRP_KEY(File, kale, out, pdb)
OPT_2GRP_KEY(Boolean, kale, out, quiet)
OPT_2GRP_KEY(Integer, kale, mc, iterations)
OPT_2GRP_KEY(Integer, kale, mc, temperature)
OPT_2GRP_KEY(String, kale, kic, topology)
OPT_2GRP_KEY(String, kale, kic, closure_move)
OPT_2GRP_KEY(String, kale, kic, breadth_move)
OPT_2GRP_KEY(IntegerVector, kale, kic, pivots)
OPT_2GRP_KEY(Real, kale, kic, closure_ratio)
OPT_2GRP_KEY(Integer, kale, kic, dump_pose)
OPT_2GRP_KEY(Integer, kale, kic, dump_stats)
// }}}1

class StepwiseTorsionMover : public Mover { // {{{1

	public:

		StepwiseTorsionMover(Size first_index, Size last_index) { // {{{2
			first_index_ = first_index;
			last_index_ = last_index;
		};

		string get_name() const { // {{{2
			return "StepwiseTorsionMover";
		}

		void apply (pose::Pose &pose) { // {{{2
			Size index = numeric::random::random_range(first_index_, last_index_);
			Size which_angle = numeric::random::random_range(1, 3);
			Real angle_value = 360 * numeric::random::uniform();

			switch (which_angle) {
				case 1: pose.set_phi(index, angle_value); break;
				case 2: pose.set_psi(index, angle_value); break;
				case 3: pose.set_omega(index, angle_value); break;
			}
		}
		// }}}2

	private:
		Size first_index_, last_index_;

};

class ConcertedTorsionMover : public Mover { // {{{1

	public:

		ConcertedTorsionMover(Size first_index, Size last_index) { // {{{2
			first_index_ = first_index;
			last_index_ = last_index;
		};

		string get_name() const { // {{{2
			return "StepwiseTorsionMover";
		}

		void apply (pose::Pose &pose) { // {{{2
			for (int i = first_index_; i <= last_index_; ++i) {
				Real phi = 360 * numeric::random::uniform();
				Real psi = 360 * numeric::random::uniform();
				Real omega = 360 * numeric::random::uniform();

				pose.set_phi(i, phi);
				pose.set_psi(i, psi);
				pose.set_omega(i, omega);
			}
		}
		// }}}2

	private:
		Size first_index_, last_index_;

};

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

			topology = option[OptionKeys::kale::kic::topology]();
			closure_move = option[OptionKeys::kale::kic::closure_move]();
			breadth_move = option[OptionKeys::kale::kic::breadth_move]();
			closure_ratio = option[OptionKeys::kale::kic::closure_ratio]();
			pdb_path = option[OptionKeys::kale::in::pdb]();
			first_index = option[OptionKeys::kale::kic::pivots]()[1];
			cut_index = option[OptionKeys::kale::kic::pivots]()[2];
			last_index = option[OptionKeys::kale::kic::pivots]()[3];
			temperature = option[OptionKeys::kale::mc::temperature]();
			iterations = option[OptionKeys::kale::mc::iterations]();

			if (temperature < 0) {
				temperature = numeric_limits<Real>::infinity();
			}
		}

		void setup() { // {{{2
			import_pose::pose_from_pdb(pose, pdb_path);

			Size first = first_index;
			Size cut = cut_index;
			Size last = last_index;

			closure_mover = NULL;
			breadth_mover = NULL;

			string closure_error = "Unknown closure move '" + closure_move + "'.";
			string breadth_error = "Unknown breadth move '" + breadth_move + "'.";

			if (closure_move == "balanced")
				closure_mover = new BalancedKinematicMover(first, cut, last);
			else if (closure_move == "naive")
				closure_mover = new NaiveKinematicMover(first, cut, last);
			else if (closure_move == "off")
				closure_mover = NULL;
			else
				utility_exit_with_message(closure_error);

			if (topology == "linear") {
				if (breadth_move == "on" || breadth_move == "concerted")
					breadth_mover = new ConcertedTorsionMover(first, last);
				else if (breadth_move == "stepwise")
					breadth_mover = new StepwiseTorsionMover(first, last);
				else if (breadth_move == "off") 
					breadth_mover = NULL;
				else
					utility_exit_with_message(breadth_error);
			}

			if (topology == "cyclic") {
				using namespace devel::cycles;

				SetupMover setup_mover;
				setup_mover.apply(pose);

				if (breadth_move == "on")
					breadth_mover = new RandomReindexingMover();
				else if (breadth_move == "off") 
					breadth_mover = NULL;
				else
					utility_exit_with_message(breadth_error);
			}
		}

		void update() { // {{{2

			if (closure_mover && breadth_mover) {
				if (numeric::random::uniform() < closure_ratio) 
					closure_mover->apply(pose);
				else
					breadth_mover->apply(pose);
			}
			else if (closure_mover) {
				closure_mover->apply(pose);
			}
			else if (breadth_mover) {
				breadth_mover->apply(pose);
			}
			else {
				utility_exit_with_message("No movers instantiated.");
			}
		}

		int get_iterations() { // {{{2
			return iterations;
		}
		// }}}2

	private:

		// Pose and movers {{{2
		pose::Pose pose;
		KinematicMoverOP closure_mover;
		MoverOP breadth_mover;

		// Command line options {{{2
		string topology;
		string closure_move;
		string breadth_move;
		Real closure_ratio;
		string pdb_path;
		Size first_index, cut_index, last_index;
		Size iterations;
		Real temperature;
		// }}}2

};

class OutputManager { // {{{1

	public:

		OutputManager(SamplingManager &master) : sampler(master) { // {{{2
			quiet = option[OptionKeys::kale::out::quiet]();
		  dump_pose_freq = option[OptionKeys::kale::kic::dump_pose]();
		  dump_stats_freq = option[OptionKeys::kale::kic::dump_stats]();

			log_coordinates.open("coordinates.dat");
			log_solutions.open("solutions.dat");
			log_pivots.open("pivots.dat");
		}

		~OutputManager() { // {{{2
			log_coordinates.close();
			log_solutions.close();
			log_pivots.close();
		}
		// }}}2

		void write_header() { // {{{2
			cout << "Topology:      " << sampler.topology << endl;
			cout << "Closure Move:  " << sampler.closure_move << endl;
			cout << "Breadth Move:  " << sampler.breadth_move << endl;
			cout << "Closure Ratio: " << sampler.closure_ratio << endl;
			cout << "Peptide:       " << sampler.pdb_path << endl;
			cout << "Pivots:        " << sampler.first_index << "/"
				                        << sampler.cut_index << "/"
																<< sampler.last_index << endl;
			cout << "Random Seed:   " << numeric::random::rg().get_seed() << endl;
			cout << "Iterations:    " << sampler.iterations << endl;
			cout << "Temperature:   " << sampler.temperature << endl;

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

			if (sampler.closure_mover != NULL) {
				Size solution_count = sampler.closure_mover->solutions_to_last_move();
				log_coordinates << endl << "solutions  " << solution_count;
			}

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
	NEW_OPT(kale::kic::topology, "Linear or cyclic", "");
	NEW_OPT(kale::kic::closure_move, "Closure move", "");
	NEW_OPT(kale::kic::breadth_move, "Breadth move", "");
	NEW_OPT(kale::kic::pivots, "Pivot residues", 0);
	NEW_OPT(kale::kic::closure_ratio, "Closure to breadth move ratio", 0.5);
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

		// Optionally:
		// 1. Repacking (every 20 steps)
		// 2. Rotamer trials (what is this?)
		// 3. Gradient minimization

		// Acceptance criterion logic.

		output.write_progress(i);
		output.dump_pose(i);
		output.dump_statistics(i);
	}

	output.write_footer();
	output.dump_pose(iterations, true);
}
// }}}1
