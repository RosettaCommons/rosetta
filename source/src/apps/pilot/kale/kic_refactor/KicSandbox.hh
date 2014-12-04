// Core headers {{{1
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Protocol headers {{{1
#include <protocols/loops/Loop.hh>
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopProtocol.hh>
#include <protocols/loop_modeling/samplers/LegacyKicSampler.hh>
#include <protocols/loop_modeling/refiners/MinimizationRefiner.hh>
#include <protocols/loop_modeling/refiners/RotamerTrialsRefiner.hh>
#include <protocols/loop_modeling/refiners/RepackingRefiner.hh>
#include <protocols/kinematic_closure/KicMover.hh>

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

typedef utility::pointer::shared_ptr<KicSandbox> KicSandboxOP;
typedef utility::pointer::shared_ptr<KicSandbox const> KicSandboxCOP;

KicSandbox::KicSandbox(int argc, char *argv[]) {

	// Namespaces (within main) {{{1
	using core::Size;
	using core::import_pose::pose_from_pdb;

	using protocols::loops::Loop;
	using protocols::loop_modeling::LoopProtocol;
	using protocols::loop_modeling::samplers::LegacyKicSampler;
	using protocols::loop_modeling::refiners::MinimizationRefiner;
	using protocols::loop_modeling::refiners::RotamerTrialsRefiner;
	using protocols::loop_modeling::refiners::RepackingRefiner;
	using protocols::kinematic_closure::KicMover;

	using namespace basic::options;

	// Options (within main) {{{1
	IndexList default_iterations(3, 1);

	NEW_OPT(kale::algorithm, "Sampling algorithm", "balanced-unrefined");
	NEW_OPT(kale::structure, "Input PDB file", "bamf");
	NEW_OPT(kale::loop, "Residues to sample", 0);
	NEW_OPT(kale::iterations, "Iterations in each loop", default_iterations);

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

	cout << "Algorithm:     " << algorithm << endl;
	cout << "Structure:     " << pdb_path << endl;
	cout << "Loop Indices:  " << loop_indices[1] << " to "
	                          << loop_indices[2] << endl;
	cout << "Random Seed:   " << numeric::random::rg().get_seed() << endl;
	cout << "Iterations:    " << iterations[1] << "/"
	                          << iterations[2] << "/"
	                          << iterations[3] << endl;
	// }}}1

	protocol = LoopProtocolOP( new LoopProtocol );

	// Configure the pose {{{1

	pose_from_pdb(pose, pdb_path);

	// Configure the loop region {{{1

	loop = Loop(loop_indices[1], loop_indices[2], loop_indices[2]);

	// Configure the sampling algorithm {{{1
	
	if (algorithm == "refactor") {
		protocol->add_mover(LoopMoverOP( new KicMover ));
		protocol->add_mover(LoopMoverOP( new RepackingRefiner(20) ));
		protocol->add_mover(LoopMoverOP( new RotamerTrialsRefiner ));
		protocol->add_mover(LoopMoverOP( new MinimizationRefiner ));
	}
	else if (algorithm == "refactor-only") {
		protocol->add_mover(LoopMoverOP( new KicMover ));
	}
	else {
		utility_exit_with_message("Unknown algorithm: " + algorithm);
	}

	// Finish configuring the loop modeling protocol {{{1

	protocol->set_loop(loop);
	//protocol->set_iterations(iterations);
	
	// }}}1
}

}
}

