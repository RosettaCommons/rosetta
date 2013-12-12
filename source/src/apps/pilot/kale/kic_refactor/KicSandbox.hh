#pragma once

// Core headers {{{1
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Protocol headers {{{1
#include <protocols/loops/Loop.hh>
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopProtocol.hh>
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/loop_modeling/LoopMoverTask.hh>
#include <protocols/loop_modeling/samplers/KicSampler.hh>
#include <protocols/loop_modeling/samplers/BalancedKicSampler.hh>
#include <protocols/loop_modeling/samplers/LegacyKicSampler.hh>
#include <protocols/loop_modeling/refiners/LocalMinimizationRefiner.hh>
#include <protocols/loop_modeling/refiners/RotamerTrialsRefiner.hh>
#include <protocols/loop_modeling/refiners/RepackingRefiner.hh>
#include <protocols/loop_modeling/utilities/RepeatedTask.hh>
#include <protocols/loop_modeling/utilities/PeriodicTask.hh>
#include <protocols/loop_modeling/loggers/Logger.hh>
#include <protocols/loop_modeling/loggers/ProgressBar.hh>
#include <protocols/loop_modeling/loggers/AcceptanceRates.hh>
#include <protocols/loop_modeling/loggers/PdbLogger.hh>
#include <protocols/loop_modeling/loggers/ScoreVsRmsd.hh>

// Utility headers {{{1
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/exit.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <numeric/random/random.hh>
#include <devel/init.hh>

// C++ headers {{{1
#include <iostream>
// }}}1

// Options {{{1
OPT_1GRP_KEY(String, kale, algorithm)
OPT_1GRP_KEY(File, kale, structure)
OPT_1GRP_KEY(IntegerVector, kale, loop)
OPT_1GRP_KEY(IntegerVector, kale, iterations)
OPT_2GRP_KEY(Boolean, kale, out, mc_poses)
OPT_2GRP_KEY(Boolean, kale, out, all_poses)
OPT_2GRP_KEY(Boolean, kale, out, progress_bar)
// }}}1

namespace apps {
namespace pilot {

// Namespaces {{{1
using namespace std;
using core::pose::Pose;
using protocols::loops::Loop;
using protocols::loop_modeling::LoopProtocolOP;
using protocols::loop_modeling::LoopMoverOP;
using utility::pointer::ReferenceCount;
using protocols::loop_modeling::IndexList;
// }}}1

class KicSandbox : public ReferenceCount {

public: 
	KicSandbox(int argc, char *argv[]);

public:
	Pose pose;
	Loop loop;
	LoopMoverOP mover;
	LoopProtocolOP protocol;
	IndexList iterations;
};

typedef utility::pointer::owning_ptr<KicSandbox> KicSandboxOP;
typedef utility::pointer::owning_ptr<KicSandbox const> KicSandboxCOP;

KicSandbox::KicSandbox(int argc, char *argv[]) {

	// Namespaces (within main) {{{1
	using core::Size;
	using core::import_pose::pose_from_pdb;

	using protocols::loops::Loop;
	using protocols::loop_modeling::LoopProtocol;
	using protocols::loop_modeling::LoopMover;
	using protocols::loop_modeling::samplers::KicSampler;
	using protocols::loop_modeling::samplers::BalancedKicSampler;
	using protocols::loop_modeling::samplers::LegacyKicSampler;
	using protocols::loop_modeling::refiners::LocalMinimizationRefiner;
	using protocols::loop_modeling::refiners::RotamerTrialsRefiner;
	using protocols::loop_modeling::refiners::RepackingRefiner;
	using protocols::loop_modeling::utilities::RepeatedTask;
	using protocols::loop_modeling::utilities::PeriodicTask;
	using protocols::loop_modeling::loggers::LoggerOP;
	using protocols::loop_modeling::loggers::ProgressBar;
	using protocols::loop_modeling::loggers::AcceptanceRates;
	using protocols::loop_modeling::loggers::PdbLogger;
	using protocols::loop_modeling::loggers::ScoreVsRmsd;

	using namespace basic::options;

	// Options (within main) a::a {{{1
	IndexList default_iterations(3, 1);

	NEW_OPT(kale::algorithm, "Sampling algorithm", "balanced-unrefined");
	NEW_OPT(kale::structure, "Input PDB file", "bamf");
	NEW_OPT(kale::loop, "Residues to sample", 0);
	NEW_OPT(kale::iterations, "Iterations in each loop", default_iterations);
	NEW_OPT(kale::out::mc_poses, "Log a pose after every MC check", true);
	NEW_OPT(kale::out::all_poses, "Log a pose after every change", true);
	NEW_OPT(kale::out::progress_bar, "Draw a progress bar", true);

	devel::init(argc, argv);

	if (option[OptionKeys::kale::algorithm].active() == false) {
		utility_exit_with_message("No sampling algorithm specified.");
	}

	if (option[OptionKeys::kale::structure].user() == false) {
		utility_exit_with_message("No input PDB file specified.");
	}

	if (option[OptionKeys::kale::loop].user() == false) {
		utility_exit_with_message("No loop region specified.");
	}

	string algorithm = option[OptionKeys::kale::algorithm]();
	string pdb_path = option[OptionKeys::kale::structure]();
	IndexList loop_indices = option[OptionKeys::kale::loop]();
	iterations = option[OptionKeys::kale::iterations]();
	bool log_mc_poses = option[OptionKeys::kale::out::mc_poses]();
	bool log_all_poses = option[OptionKeys::kale::out::all_poses]();
	bool progress_bar = option[OptionKeys::kale::out::progress_bar]();

	cout << "Algorithm:     " << algorithm << endl;
	cout << "Structure:     " << pdb_path << endl;
	cout << "Loop Indices:  " << loop_indices[1] << " to "
	                          << loop_indices[2] << endl;
	cout << "Random Seed:   " << numeric::random::RG.get_seed() << endl;
	cout << "Iterations:    " << iterations[1] << "/"
	                          << iterations[2] << "/"
	                          << iterations[3] << endl;
	// }}}1

	protocol = new LoopProtocol;

	// Configure the pose {{{1

	pose_from_pdb(pose, pdb_path);

	// Configure the loop region {{{1

	loop = Loop(loop_indices[1], loop_indices[2]);

	// Configure the logging output {{{1

	LoggerOP acceptance_logger = new AcceptanceRates;

	protocol->add_logger(new ScoreVsRmsd(pose, loop));
	protocol->add_logger(acceptance_logger);

	if (progress_bar) protocol->add_logger(new ProgressBar);
	if (log_all_poses) protocol->add_logger(new PdbLogger(true));
	else if (log_mc_poses) protocol->add_logger(new PdbLogger(false));

	// Configure the sampling algorithm {{{1
	
	mover = new LoopMover;

	if (algorithm == "refactor") {
		mover->add_task(new KicSampler(acceptance_logger));
		mover->add_task(new PeriodicTask(new RepackingRefiner, 20));
		mover->add_task(new RotamerTrialsRefiner);
		mover->add_task(new LocalMinimizationRefiner);
	}
	else if (algorithm == "refactor-only") {
		mover->add_task(new KicSampler(acceptance_logger));
	}
	/*
	else if (algorithm == "balanced") {
		mover->add_task(new BalancedKicSampler);
	}
	else if (algorithm == "balanced-sidechain") {
		mover->add_task(new BalancedKicSampler);
		mover->add_task(new BalancedSidechainRefiner);
	}
	else if (algorithm == "balanced-config-bias") {
		mover->add_task(new BalancedKicSampler);
		mover->add_task(new ConfigurationalBiasRefiner);
	}
	*/
	else {
		utility_exit_with_message("Unknown algorithm: " + algorithm);
	}

	// Finish configuring the loop modeling protocol {{{1

	protocol->set_loop(loop);
	protocol->set_mover(mover);
	protocol->set_iterations(iterations);
	// }}}1
	
}

}
}

