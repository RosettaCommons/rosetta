// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Headers {{{1
#include <protocols/loop_modeling/LoopModeler.hh>
#include <protocols/loop_modeling/LoopModelerCreator.hh>
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/loop_modeling/LoopBuilder.hh>
#include <protocols/loop_modeling/LoopProtocol.hh>
#include <protocols/loop_modeling/utilities/rosetta_scripts.hh>
#include <protocols/loop_modeling/utilities/LoopMoverGroup.hh>
#include <protocols/loop_modeling/utilities/AcceptanceCheck.hh>

// Protocols headers
#include <protocols/filters/Filter.hh>
#include <protocols/kinematic_closure/KicMover.hh>
#include <protocols/kinematic_closure/perturbers/OmegaPerturber.hh>
#include <protocols/kinematic_closure/perturbers/Rama2bPerturber.hh>
#include <protocols/kinematic_closure/perturbers/FragmentPerturber.hh>
#include <protocols/loop_modeling/refiners/MinimizationRefiner.hh>
#include <protocols/loop_modeling/refiners/RepackingRefiner.hh>
#include <protocols/loop_modeling/refiners/RotamerTrialsRefiner.hh>
#include <protocols/loop_modeling/utilities/PrepareForCentroid.hh>
#include <protocols/loop_modeling/utilities/PrepareForFullatom.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/loop_mover/LoopMover.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/util.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/simple_filters/BuriedUnsatHbondFilter.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>

// Core headers
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Conformation.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/util/SwitchResidueTypeSet.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <boost/algorithm/string.hpp>

// C++ headers
#include <iostream>
#include <cmath>
#include <ctime>

// Namespaces {{{1
using namespace std;
using namespace basic::options;

using core::Real;
using core::Size;
using core::chemical::CENTROID;
using core::chemical::FA_STANDARD;
using core::pose::Pose;
using core::import_pose::pose_from_file;
using core::scoring::ScoreFunctionCOP;
using core::scoring::ScoreFunctionOP;
using core::pack::task::TaskFactoryOP;
using core::pack::task::TaskFactoryCOP;

using protocols::filters::FilterOP;
using protocols::loop_modeling::utilities::PrepareForCentroid;
using protocols::loop_modeling::utilities::PrepareForCentroidOP;
using protocols::loop_modeling::utilities::PrepareForFullatom;
using protocols::loop_modeling::utilities::PrepareForFullatomOP;
using protocols::loops::Loop;
using protocols::loops::Loops;
using protocols::moves::MonteCarloOP;
using protocols::moves::MoverOP;

using utility::tag::TagCOP;
using basic::datacache::DataMap;
using protocols::filters::Filters_map;
using protocols::moves::Movers_map;
// }}}1

namespace protocols {
namespace loop_modeling {

static basic::Tracer TR("protocols.loop_modeling.LoopModeler");

MoverOP LoopModelerCreator::create_mover() const { // {{{1
	return MoverOP( new LoopModeler );
}

string LoopModelerCreator::keyname() const { // {{{1
	return "LoopModeler";
}
// }}}1

LoopModeler::LoopModeler() { // {{{1
	build_stage_ = add_child( LoopBuilderOP( new LoopBuilder ) );
	centroid_stage_ = add_child( LoopProtocolOP( new LoopProtocol ) );
	fullatom_stage_ = add_child( LoopProtocolOP( new LoopProtocol ) );
	prepare_for_centroid_ = add_child( PrepareForCentroidOP( new PrepareForCentroid ) );
	prepare_for_fullatom_ = add_child( PrepareForFullatomOP( new PrepareForFullatom ) );

	is_build_stage_enabled_ = true;
	is_centroid_stage_enabled_ = true;
	is_fullatom_stage_enabled_ = true;

	set_fa_scorefxn(loops::get_fa_scorefxn());
	set_cen_scorefxn(loops::get_cen_scorefxn());

	setup_kic_config();
}

LoopModeler::~LoopModeler() = default; // {{{1

MoverOP LoopModeler::clone() const { // {{{1
	// This is dangerous.  It only works if something else is holding a shared
	// pointer to this object.  The proper thing to do would be to construct a
	// new loop modeler on the fly, but I think this is impossible because
	// there's no way to make a deep copy of a data map.
	return utility::pointer::const_pointer_cast<Mover>(shared_from_this());
}

void LoopModeler::parse_my_tag( // {{{1
	TagCOP tag,
	DataMap & data,
	Filters_map const & filters,
	Movers_map const & movers,
	Pose const & pose) {

	LoopMover::parse_my_tag(tag, data, filters, movers, pose);
	utilities::set_task_factory_from_tag(*this, tag, data);

	// Parse the 'config' option.

	string config = tag->getOption<string>("config", "kic");

	if ( config == "empty" ) {
		setup_empty_config();
	} else if ( config == "kic" ) {
		setup_kic_config();
	} else if ( config == "kic_with_frags" ) {
		setup_kic_with_fragments_config();
	} else {
		stringstream message;
		message << "Unknown <LoopModeler> config option: '" << config << "'";
		throw utility::excn::EXCN_Msg_Exception(message.str());
	}

	// Parse the 'auto_refine' option.

	if ( ! tag->getOption<bool>("auto_refine", true) ) {
		centroid_stage_->clear_refiners();
		fullatom_stage_->clear_refiners();
	}

	// Pares the 'scorefxn_fa' and 'scorefxn_cen' options.

	using protocols::rosetta_scripts::parse_score_function;

	if ( tag->hasOption("scorefxn_fa") ) {
		set_fa_scorefxn(parse_score_function(tag, "scorefxn_fa", data, ""));
	}
	if ( tag->hasOption("scorefxn_cen") ) {
		set_cen_scorefxn(parse_score_function(tag, "scorefxn_cen", data, ""));
	}

	// Parse the 'fast' option.

	if ( tag->getOption<bool>("fast", false) ) {
		centroid_stage_->mark_as_test_run();
		fullatom_stage_->mark_as_test_run();
	}

	// Parse subtags.

	for ( TagCOP subtag : tag->getTags() ) {

		// Ignore <Loop> subtags (parsed by parent class).

		if ( subtag->getName() == "Loop" ) { continue; }

		// Parse <Build> subtags.

		else if ( subtag->getName() == "Build" ) {
			build_stage_->parse_my_tag(subtag, data, filters, movers, pose);
			is_build_stage_enabled_ = ! subtag->getOption<bool>("skip", ! is_build_stage_enabled_);
		} else if ( subtag->getName() == "Centroid" ) {
			// Parse <Centroid> subtags.
			centroid_stage_->parse_my_tag(subtag, data, filters, movers, pose);
			is_centroid_stage_enabled_ = ! subtag->getOption<bool>("skip", ! is_centroid_stage_enabled_);
		} else if ( subtag->getName() == "Fullatom" ) {
			// Parse <Fullatom> subtags.
			fullatom_stage_->parse_my_tag(subtag, data, filters, movers, pose);
			is_fullatom_stage_enabled_ = ! subtag->getOption<bool>("skip", ! is_fullatom_stage_enabled_);
		} else {
			// Parse LoopMover subtags.
			LoopMoverOP centroid_mover = utilities::loop_mover_from_tag(subtag, data, filters, movers, pose);
			LoopMoverOP fullatom_mover = utilities::loop_mover_from_tag(subtag, data, filters, movers, pose);

			centroid_stage_->add_mover(centroid_mover);
			fullatom_stage_->add_mover(fullatom_mover);
		}
	}
}

bool LoopModeler::do_apply(Pose & pose) { // {{{1
	Pose native_pose(pose);
	if ( option[OptionKeys::in::file::native].user() ) {
		pose_from_file(native_pose, option[OptionKeys::in::file::native](), core::import_pose::PDB_file);
	}
	prepare_for_fullatom_->set_original_pose(pose);

	// Converting the pose to centroid mode can sometimes cause headaches with
	// non-canonical residue type sets, so this conditional makes sure that this
	// conversion only happens if necessary.

	if ( is_build_stage_enabled_ || is_centroid_stage_enabled_ ) {
		prepare_for_centroid_->apply(pose);
	}

	// Build stage

	if ( is_build_stage_enabled_ ) {
		long start_time = time(nullptr);
		TR << "Build Stage" << endl;
		build_stage_->apply(pose);
		long end_time = time(nullptr);
		TR << "Build Time: " << end_time - start_time << " sec" << endl;
	}

	// Output score and RMSD to debug tracer?  That sounds good.  I need two
	// versions of a function that outputs score and rmsd.  One version uses the
	// standard tracer, fills up the pose extra scores, and deals with buried
	// H-bonds.  The other version just outputs to a debug tracer.
	//
	// Actually, I'm starting to think that I should move the score and rmsd code
	// to LoopProtocol.  LoopModeler shouldn't do anything but combine other bits
	// of code, and right now there's no other way to get this output behavior.
	// If I were to move the code into LoopProtocol, both fullatom and centroid
	// would report scores and RMSDs to the log.  The full atom extra pose scores
	// would overwrite the centroid ones, which is good.  But LoopBuilder isn't a
	// LoopProtocol, so I will have to move all this stuff into its own function
	// in order to reuse it.  That's not so bad.  This really would be a
	// generally useful function for evaluating the quality of a loop.  The
	// argument would have to be a pose, a reference pose, a loop, and a tracer.
	//
	// Actually, only the LoopModeler knows about the original pose, doesn't it.
	// So I probably can't move this function to the LoopProtocol.  Also, I
	// should consider the fact that in real application the RMSD is meaningless
	// and should probably not be displayed.  You know, I could have my utility
	// function query -in:file:native.  Then I don't need to pass it a reference
	// structure and I can easily omit RMSDs if there is no starting structure.
	// That's perfect!

	// Centroid stage

	if ( is_centroid_stage_enabled_ ) {
		long start_time = time(nullptr);
		TR << "Centroid Stage" << endl;
		centroid_stage_->apply(pose);
		long end_time = time(nullptr);
		TR << "Centroid Time: " << end_time - start_time << " sec" << endl;
	}

	// Fullatom stage

	if ( is_fullatom_stage_enabled_ ) {
		long start_time = time(nullptr);
		TR << "Fullatom Stage" << endl;
		prepare_for_fullatom_->apply(pose);
		fullatom_stage_->apply(pose);
		long end_time = time(nullptr);
		TR << "Fullatom Time: " << end_time - start_time << " sec" << endl;
	}

	// Report score and rmsd

	LoopsCOP loops = get_loops();
	ObjexxFCL::FArray1D_bool loop_residues (pose.size(), false);

	for ( Size i = 1; i <= pose.size(); i++ ) {
		if ( loops->is_loop_residue(i) ) { loop_residues[i-1] = true; }
	}

	using core::scoring::rmsd_no_super_subset;
	using core::scoring::is_protein_backbone;
	using core::scoring::chainbreak;

	Real total_score = pose.is_fullatom() ?
		fullatom_stage_->get_score_function()->score(pose):
		centroid_stage_->get_score_function()->score(pose);
	Real chainbreak_score = pose.energies().total_energies()[chainbreak];
	Real backbone_rmsd = rmsd_no_super_subset(
		native_pose, pose, loop_residues, is_protein_backbone);

	TR << "Total Score: " << total_score << endl;
	TR << "Chainbreak Score: " << chainbreak_score << endl;
	TR << "Loop Backbone RMSD: " << backbone_rmsd << endl;

	core::pose::setPoseExtraScore(pose, "total_score", total_score);
	core::pose::setPoseExtraScore(pose, "chainbreak_score", chainbreak_score);
	core::pose::setPoseExtraScore(pose, "loop_backbone_rmsd", backbone_rmsd);

	if ( pose.is_fullatom() ) {
		protocols::simple_filters::BuriedUnsatHbondFilter unsats(20, 0);
		Real delta_unsats = unsats.compute(pose) - unsats.compute(native_pose);
		TR << "Delta Buried Unsatisfied H-bonds: " << showpos << delta_unsats << endl;
		core::pose::setPoseExtraScore(pose, "delta_buried_unsats", delta_unsats);
	}

	return true;
}

void LoopModeler::setup_empty_config() { // {{{1
	centroid_stage_->set_sfxn_cycles(5);
	centroid_stage_->set_temp_cycles(20, true);
	centroid_stage_->set_mover_cycles(1);
	centroid_stage_->set_temperature_schedule(2.0, 1.0);
	centroid_stage_->clear_movers_and_refiners();

	fullatom_stage_->set_sfxn_cycles(5);
	fullatom_stage_->set_temp_cycles(10, true);
	fullatom_stage_->set_mover_cycles(2);
	fullatom_stage_->clear_movers_and_refiners();
}

void LoopModeler::setup_kic_config() { // {{{1
	using namespace protocols::kinematic_closure;
	using namespace protocols::loop_modeling::refiners;
	using namespace protocols::loop_modeling::utilities;

	setup_empty_config();

	centroid_stage_->add_mover(LoopMoverOP( new KicMover ));
	centroid_stage_->add_refiner(LoopMoverOP( new MinimizationRefiner ));
	centroid_stage_->mark_as_default();

	fullatom_stage_->set_temperature_ramping(true);
	fullatom_stage_->set_rama_term_ramping(true);
	fullatom_stage_->set_repulsive_term_ramping(true);
	fullatom_stage_->add_mover(LoopMoverOP( new KicMover ));
	fullatom_stage_->add_refiner(LoopMoverOP( new RepackingRefiner(40) ));
	fullatom_stage_->add_refiner(LoopMoverOP( new RotamerTrialsRefiner ));
	fullatom_stage_->add_refiner(LoopMoverOP( new MinimizationRefiner ));
	fullatom_stage_->mark_as_default();
}

void LoopModeler::setup_kic_with_fragments_config() { // {{{1
	using namespace basic::options;
	using namespace protocols::kinematic_closure;
	using namespace protocols::kinematic_closure::perturbers;

	setup_kic_config();

	// Read fragment data from the command line.

	vector1<core::fragment::FragSetOP> frag_libs;

	if ( option[ OptionKeys::loops::frag_files ].user() ) {
		loops::read_loop_fragments(frag_libs);
	} else {
		stringstream message;
		message << "Must specify the -loops:frag_sizes and -loops:frag_files ";
		message << "options in order to use the FragmentPerturber." << endl;
		throw utility::excn::EXCN_Msg_Exception(message.str());
	}

	// Enable fragments during the build stage.

	build_stage_->use_fragments(frag_libs);

	// Create a centroid "KIC with fragments" mover (see note).

	KicMoverOP centroid_kic_mover( new KicMover );
	centroid_kic_mover->add_perturber(perturbers::PerturberOP( new FragmentPerturber(frag_libs) ));
	centroid_stage_->add_mover(centroid_kic_mover);
	centroid_stage_->mark_as_default();

	// Create a fullatom "KIC with fragments" mover (see note).

	KicMoverOP fullatom_kic_mover( new KicMover );
	fullatom_kic_mover->add_perturber(perturbers::PerturberOP( new FragmentPerturber(frag_libs) ));
	fullatom_stage_->add_mover(fullatom_kic_mover);
	fullatom_stage_->mark_as_default();

	// Note: Because loop movers can query their parents for attributes like
	// score functions and task operations, no loop mover can have more than one
	// parent.  This is why two separate centroid and fullatom KicMover objects
	// must be created.
}

LoopBuilderOP LoopModeler::build_stage() { // {{{1
	return build_stage_;
}

LoopProtocolOP LoopModeler::centroid_stage() { // {{{1
	return centroid_stage_;
}

LoopProtocolOP LoopModeler::fullatom_stage() { // {{{1
	return fullatom_stage_;
}

void LoopModeler::enable_build_stage() { // {{{1
	is_build_stage_enabled_ = true;
}

void LoopModeler::disable_build_stage() { // {{{1
	is_build_stage_enabled_ = false;
}

void LoopModeler::enable_centroid_stage() { // {{{1
	is_centroid_stage_enabled_ = true;
}

void LoopModeler::disable_centroid_stage() { // {{{1
	is_centroid_stage_enabled_ = false;
}

void LoopModeler::enable_fullatom_stage() { // {{{1
	is_fullatom_stage_enabled_ = true;
}

void LoopModeler::disable_fullatom_stage() { // {{{1
	is_fullatom_stage_enabled_ = false;
}

TaskFactoryOP LoopModeler::get_task_factory() { // {{{1
	return get_tool<TaskFactoryOP>(ToolboxKeys::TASK_FACTORY);
}

TaskFactoryOP LoopModeler::get_task_factory(TaskFactoryOP fallback) { // {{{1
	return get_tool<TaskFactoryOP>(ToolboxKeys::TASK_FACTORY, fallback);
}

void LoopModeler::set_task_factory(TaskFactoryOP task_factory) { // {{{1
	set_tool(ToolboxKeys::TASK_FACTORY, task_factory);
}

void LoopModeler::set_fa_scorefxn(ScoreFunctionOP function) { // {{{1
	prepare_for_fullatom_->set_score_function(function);
	fullatom_stage_->set_score_function(function);
}

void LoopModeler::set_cen_scorefxn(ScoreFunctionOP function) { // {{{1
	build_stage_->set_score_function(function);
	centroid_stage_->set_score_function(function);
}

// }}}1

}
}

