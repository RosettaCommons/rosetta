// Comments {{{1
//
// Input
// 1. PDB structure.
// 2. Loop residues (optional): If not given, all residues will be sampled.  No 
//    design will be done either way.
// 3. Algorithm (optional): String that specifies how the simulation should be 
//    configured.  Don't foresee needing anything more complex at the moment.
// 4. Backbone:sidechain move ratio (optional): Defaults to 1:10
//
// Output
// 1. Trajectory.  This will be large.  Output to database.
// 2. Entropy score.  Maybe later.

// Should I write a new application, or should I make the old one work?
// 1. No need for breadth mover.
// 2. No need for unbalanced moves.
// 3. Want to start using database.
// 4. Will take a while to get going.
// 5. Have to add moving pivots.

// Headers {{{1

// Core headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/kinematic_closure/samplers/BalancedKicSampler.hh>
#include <protocols/kinematic_closure/perturbers/RamaPerturber.hh>
#include <protocols/kinematic_closure/perturbers/WalkingPerturber.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/jd2/JobDistributor.hh>

// Utility headers
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <numeric/random/random.hh>
#include <utility/vector1.hh>
#include <boost/foreach.hpp>

// C++ headers
#include <string>

#define foreach BOOST_FOREACH

// Global Names {{{1

using namespace std;
using namespace basic::options;

using core::Size;
using core::Real;
using core::pose::Pose;
using core::import_pose::pose_from_pdb;
using core::pack::task::TaskFactory;
using core::pack::task::TaskFactoryOP;
using core::pack::task::operation::RestrictToRepacking;
using core::pack::task::operation::PreventRepacking;
using core::pack::task::operation::PreventRepackingOP;
using core::scoring::ScoreFunction;
using core::scoring::ScoreFunctionOP;

using protocols::loops::Loop;
using protocols::moves::Mover;
using protocols::moves::MoverOP;
using protocols::moves::MonteCarlo;
using protocols::moves::MonteCarloOP;
using protocols::simple_moves::sidechain_moves::SidechainMover;
using protocols::simple_moves::sidechain_moves::SidechainMoverOP;
using protocols::kinematic_closure::samplers::BalancedKicSampler;
using protocols::kinematic_closure::samplers::BalancedKicSamplerOP;
using protocols::kinematic_closure::perturbers::RamaPerturber;
using protocols::kinematic_closure::perturbers::WalkingPerturber;

typedef utility::vector1<Size> IndexList;
basic::Tracer tr("apps.native_ensemble");

// Options {{{1

OPT_1GRP_KEY(IntegerVector, native_ensemble, loop)
OPT_1GRP_KEY(IntegerVector, native_ensemble, weights)
OPT_1GRP_KEY(String, native_ensemble, algorithm)
OPT_1GRP_KEY(Integer, native_ensemble, iterations)
OPT_1GRP_KEY(Real, native_ensemble, temperature)
OPT_1GRP_KEY(Boolean, native_ensemble, quiet)

// }}}1

class NativeEnsemble; // {{{1

typedef utility::pointer::owning_ptr<NativeEnsemble> NativeEnsembleOP;
typedef utility::pointer::owning_ptr<NativeEnsemble const> NativeEnsembleCOP;

class NativeEnsemble : public Mover {

public:
	NativeEnsemble()
		: loop_(0, 0),
		  algorithm_("rama"),
			iterations_(1000),
			temperature_(1),
			backbone_ratio_(0.1),
			quiet_(false) {}

	void apply(Pose & pose);
	string get_name() const { return "Native Ensemble"; }
	MoverOP fresh_instance() const { return new NativeEnsemble(*this); }

	void set_loop(Loop loop) { loop_ = loop; }
	void set_loop(IndexList values) { set_loop(Loop(values[1], values[2])); }
	void set_algorithm(string value) { algorithm_ = value; }
	void set_iterations(Size value) { iterations_ = value; }
	void set_temperature(Real value) { temperature_ = value; }
	void set_weights(Real bb, Real sc) { backbone_ratio_ = bb / (bb + sc); }
	void set_weights(IndexList values) { set_weights(values[1], values[2]); }
	void set_quiet(bool value) { quiet_ = value; }

private:
	void show_header() const;
	void show_progress(Size current_iteration) const;

private:
	Loop loop_;
	string algorithm_;
	Size iterations_;
	Real temperature_;
	Real backbone_ratio_;
	bool quiet_;
};

void NativeEnsemble::apply(Pose & pose) { // {{{1

	// Setup the monte carlo run.

	ScoreFunctionOP score_function = core::scoring::getScoreFunction();
	MonteCarloOP monte_carlo = new MonteCarlo(pose, *score_function, 1);
	if (loop_.start() <= 1) loop_ = Loop(2, pose.total_residue() - 1);

	show_header();

	// Setup the backbone mover.

	BalancedKicSamplerOP backbone_mover = new BalancedKicSampler;

	if (algorithm_ == "rama")
		backbone_mover->add_perturber(new RamaPerturber);
	else if (algorithm_ == "walking")
		backbone_mover->add_perturber(new WalkingPerturber);
	else
		utility_exit_with_message("Unknown algorithm: " + algorithm_);

	// Setup the sidechain mover.

	SidechainMoverOP sidechain_mover = new SidechainMover;
	TaskFactoryOP task_factory = new TaskFactory;
	PreventRepackingOP prevent_repacking = new PreventRepacking();

	for (Size i = 1; i <= pose.total_residue(); i++) {
		if (i < loop_.start() or i > loop_.stop()) {
			prevent_repacking->include_residue(i);
		}
	}

	task_factory->push_back(new RestrictToRepacking);
	task_factory->push_back(prevent_repacking);

	sidechain_mover->set_preserve_detailed_balance(true);
	sidechain_mover->set_task_factory(task_factory);

	// Run the monte carlo loop.

	for (int i = 1; i <= iterations_; i++) {
		string move_name;
		Real proposal_ratio;

		if (backbone_ratio_ < numeric::random::uniform()) {
			backbone_mover->apply(pose, loop_);
			move_name = "backbone";
			proposal_ratio = 1;
		}
		else {
			sidechain_mover->apply(pose);
			move_name = "sidechain";
			proposal_ratio = sidechain_mover->last_proposal_density_ratio();
		}

		monte_carlo->boltzmann(pose, move_name, proposal_ratio);
		show_progress(i);
	}

	cerr << endl;
}

void NativeEnsemble::show_header() const { // {{{1
	tr << "Residues:    " << loop_.start() << "/" << loop_.stop() << endl;
	tr << "Algorithm:   " << algorithm_ << endl;
	tr << "Iterations:  " << iterations_ << endl;
	tr << "Temperature: " << temperature_ << endl;
	tr << "% BB Moves:  " << backbone_ratio_ << endl;
}

void NativeEnsemble::show_progress(Size current_iteration) const { // {{{1
	if (quiet_ == true) return;
	cerr << "\r[" << current_iteration << "/" << iterations_ << "]";
}

// }}}1

int main(int argc, char * argv[]) { // {{{1
	IndexList default_loop(2);
	default_loop[1] = 0;
	default_loop[2] = 0;

	IndexList default_weights(2);
	default_weights[1] = 1;
	default_weights[2] = 9;

	NEW_OPT(native_ensemble::loop, "Residues to sample", default_loop);
	NEW_OPT(native_ensemble::weights, "BB-to-SC move ratio", default_weights);
	NEW_OPT(native_ensemble::algorithm, "Backbone sampling algorithm", "rama");
	NEW_OPT(native_ensemble::iterations, "Monte Carlo iterations", 1000);
	NEW_OPT(native_ensemble::temperature, "Monte Carlo temperature", 1);
	NEW_OPT(native_ensemble::quiet, "Hide progress updates", false);

	devel::init(argc, argv);

	NativeEnsembleOP app = new NativeEnsemble;
	Pose pose;

	if (not option[OptionKeys::in::file::s].active()) {
		utility_exit_with_message("Specify an input PDB using the '-s' flag.");
	}

	pose_from_pdb(pose, option[OptionKeys::in::file::s][1].name());

	cout << option[OptionKeys::native_ensemble::loop].active() << endl;
	if (option[OptionKeys::native_ensemble::loop].active()) {
		app->set_loop(option[OptionKeys::native_ensemble::loop]());
	}

	app->set_weights(option[OptionKeys::native_ensemble::weights]());
	app->set_algorithm(option[OptionKeys::native_ensemble::algorithm]());
	app->set_iterations(option[OptionKeys::native_ensemble::iterations]());
	app->set_temperature(option[OptionKeys::native_ensemble::temperature]());
	app->set_quiet(option[OptionKeys::native_ensemble::quiet]());

	//protocols::jd2::JobDistributor::get_instance()->go(app);
	app->apply(pose);
}
// }}}1
