// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#define private public
#define protected public

// Headers {{{1
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/protocols/loop_modeling/utilities.hh>

// Core headers
#include <core/pose/Pose.hh>

// Protocol headers
#include <protocols/filters/Filter.hh>
#include <protocols/filters/BasicFilters.hh>
#include <protocols/filters/FilterFactory.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>
#include <protocols/rosetta_scripts/RosettaScriptsParser.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/loop_modeling/LoopBuilder.hh>
#include <protocols/loop_modeling/LoopModeler.hh>
#include <protocols/loop_modeling/LoopProtocol.hh>
#include <protocols/loop_modeling/utilities/LoopMoverGroup.hh>
#include <protocols/loop_modeling/utilities/PeriodicMover.hh>

// Utility headers
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/exit.hh>

// C++ headers
#include <string>
#include <iostream>

// Namespaces {{{1
using namespace std;
using namespace protocols::loop_modeling;
using namespace protocols::rosetta_scripts;

using utility::excn::EXCN_Msg_Exception;

// Brainstorming {{{1
// Look at the TaskOperationRegistrator class to see how to register KIC 
// helpers.  (TaskOperation doesn't derive from Mover, but the Registrator/ 
// Factory machinery is actually agnostic to all that.)  The question is still 
// how to distinguish perturbers from pivot_pickers from solution_pickers.
//
// There are three classes: Registrator, Factory, and Creator
//
// Registrator is just a C++ trick to associate creators with factories before 
// main.
//
// Factory is a singleton that manages a group of creators.  For example, there 
// is a MoverFactory and a separate Creator for each mover.  Each creator has a 
// method create an object and a method to provide a key name.  For each 
// creator it receives, the factory looks at the key name and uses it to 
// construct a key-name -> creator map, which it can use later on to 
// instantiate things.
//
// Right now, the factories just return default constructed objects and 
// parse_my_tag() is called to set them up.  It's nice that the factory doesn't 
// depend on rosetta scripts.
//
// For KIC, I don't know what type of object I'm going to get if I pass the 
// factory a name.  I could create a type for all KIC helpers, that would 
// basically have just parse_my_tag(tag, kic_mover) that would setup the mover 
// and add it to the kic mover.  Or I could split that into two functions: 
// parse_my_tag(tag) and add_to_kic_mover().  The former could be written once 
// for each thing.  There's not even one KicMover, though.  There are 3, and 
// they don't have all the same methods.  Could write one add_to_kic_mover() 
// for each, although this sounds like the kind of can't-redefine-overloaded- 
// methods-in-subclass issues that I might have run across before.
//
// If I have a different factory for perturbers, pivot_pickers, and 
// solution_pickers, I'll have to manually deal with misspelled names.  But 
// appealing because KicMover.parse_my_tag() can do things like: No 
// solution_pickers for balanced kic mover.
//
// Only one solution picker and pivot picker.
//
// Could do:
// <KicMover>
//   <Perturber name=rama options=...>
//   <Perturber name=fragments options=...>
//   <SolutionPicker name=standard options=...>
// </KicMover>
//
// or:
// <KicMover>
//   <Perturbers>
//     <RamaPerturber options=...>
//     <FragmentPerturber options=...>
//   </Perturbers>
//   <SolutionPicker>		// Misleading because more than one picker is illegal.
//     <StandardPivots options=...>
//   </SolutionPicker>
// </KicMover>
//
// As opposed to:
// <KicMover>
//   <RamaPerturber options=...>
//   <FragmentPerturber options=...>
//   <StandardPivots options=...>
// </KicMover>
//
// The thing is that specifying a pivot picker or a solution picker won't 
// change the perturbers, and likewise specifying a perturber won't affect the 
// pickers.  So maybe there's value in clearly distinguishing these things.
//
// A heuristic approach: Perturbers must end with "Perturber", pivot pickers 
// must end with "Pivots", and solution pickers must end with "Solutions".  
// Could be enforced by a single factory which manages three types of creators.  
// Then KicMover.parse_my_tag() can know what type to ask for.  I'll have to 
// rename IdealizeNonPhiPsi and PerturberSet.  This approach seems inflexible, 
// but from a YAGNI perspective it might be the way to go.  The KicMovers will 
// have to do the string manip logic, though.
// 
// I could have the factory have a method that takes a name and returns the 
// type as an enum.  Then the KicMover can drop into a switch and request the 
// type it wants.  Or the factory could return a ReferenceCount object (i.e.  
// base class) and an enum type, and KicMover could drop into a switch to cast.
// 
// I could have the factory just build up lists of perturbers etc.  But then I 
// won't know which tag goes with which instance.
//
// I think the cleanest system is to have the subclasses manage everything 
// themselves.  Best not to even know what class I'm dealing with.  Now the 
// question is: what do I call it?  ClosureAlgorithm?  KicPlugin.  That's not 
// bad.  
//
// Anyways, the poit is that the whole registrator system is flexible and not 
// limited to movers.
// }}}1

class LoopModelerTests : public CxxTest::TestSuite {

public:

	void setUp() { protocols_init(); } // {{{1
	// }}}1

	static LoopMovers get_centroid_movers(LoopModelerOP modeler) { // {{{1
		return modeler->centroid_stage_->movers_->nested_movers_;
	}

	static LoopMovers get_centroid_refiners(LoopModelerOP modeler) { // {{{1
		return modeler->centroid_stage_->refiners_->nested_movers_;
	}

	static LoopMovers get_fullatom_movers(LoopModelerOP modeler) { // {{{1
		return modeler->fullatom_stage_->movers_->nested_movers_;
	}

	static LoopMovers get_fullatom_refiners(LoopModelerOP modeler) { // {{{1
		return modeler->fullatom_stage_->refiners_->nested_movers_;
	}
	// }}}1

	static void assert_default_build_params(LoopModelerOP modeler) { // {{{1
		TS_ASSERT_EQUALS(modeler->is_build_stage_enabled_, true);
		TS_ASSERT_EQUALS(modeler->build_stage_->max_attempts_, 10000);
	}

	static void assert_default_build_stage(LoopModelerOP modeler) { // {{{1
		assert_default_build_params(modeler);
	}

	static void assert_default_centroid_params(LoopModelerOP modeler) { // {{{1
		TS_ASSERT_EQUALS(modeler->is_centroid_stage_enabled_, true);
		TS_ASSERT_EQUALS(modeler->centroid_stage_->sfxn_cycles_, 3);
		TS_ASSERT_EQUALS(modeler->centroid_stage_->temp_cycles_, 20);
		TS_ASSERT_EQUALS(modeler->centroid_stage_->scale_temp_cycles_, true);
		TS_ASSERT_EQUALS(modeler->centroid_stage_->mover_cycles_, 1);
		TS_ASSERT_EQUALS(modeler->centroid_stage_->ramp_sfxn_rep_, false);
		TS_ASSERT_EQUALS(modeler->centroid_stage_->ramp_sfxn_rama_, false);
		TS_ASSERT_EQUALS(modeler->centroid_stage_->ramp_temp_, true);
		TS_ASSERT_EQUALS(modeler->centroid_stage_->initial_temp_, 1.5);
		TS_ASSERT_EQUALS(modeler->centroid_stage_->final_temp_, 0.5);
	}

	static void assert_default_centroid_movers(LoopModelerOP modeler) { // {{{1
		LoopMovers const & movers = get_centroid_movers(modeler);

		TS_ASSERT_EQUALS(movers.size(), 1);
		TS_ASSERT_EQUALS(movers[1]->get_name(), "KicMover");
	}

	static void assert_default_centroid_refiners(LoopModelerOP modeler) { // {{{1
		LoopMovers const & refiners = get_centroid_refiners(modeler);

		TS_ASSERT_EQUALS(refiners.size(), 1);
		TS_ASSERT_EQUALS(refiners[1]->get_name(), "MinimizationRefiner");
	}

	static void assert_default_centroid_stage(LoopModelerOP modeler) { // {{{1
		assert_default_centroid_params(modeler);
		assert_default_centroid_movers(modeler);
		assert_default_centroid_refiners(modeler);
	}

	static void assert_default_fullatom_params(LoopModelerOP modeler) { // {{{1
		TS_ASSERT_EQUALS(modeler->is_fullatom_stage_enabled_, true);
		TS_ASSERT_EQUALS(modeler->fullatom_stage_->sfxn_cycles_, 3);
		TS_ASSERT_EQUALS(modeler->fullatom_stage_->temp_cycles_, 10);
		TS_ASSERT_EQUALS(modeler->fullatom_stage_->scale_temp_cycles_, true);
		TS_ASSERT_EQUALS(modeler->fullatom_stage_->mover_cycles_, 2);
		TS_ASSERT_EQUALS(modeler->fullatom_stage_->ramp_sfxn_rep_, false);
		TS_ASSERT_EQUALS(modeler->fullatom_stage_->ramp_sfxn_rama_, false);
		TS_ASSERT_EQUALS(modeler->fullatom_stage_->ramp_temp_, true);
		TS_ASSERT_EQUALS(modeler->fullatom_stage_->initial_temp_, 1.5);
		TS_ASSERT_EQUALS(modeler->fullatom_stage_->final_temp_, 0.5);
	}

	static void assert_default_fullatom_movers(LoopModelerOP modeler) { // {{{1
		LoopMovers const & movers = get_fullatom_movers(modeler);

		TS_ASSERT_EQUALS(movers.size(), 1);
		TS_ASSERT_EQUALS(movers[1]->get_name(), "KicMover");
	}

	static void assert_default_fullatom_refiners(LoopModelerOP modeler) { // {{{1
		LoopMovers const & refiners = get_fullatom_refiners(modeler);

		TS_ASSERT_EQUALS(refiners.size(), 3);
		TS_ASSERT_EQUALS(refiners[1]->get_name(), "PeriodicMover");
		TS_ASSERT_EQUALS(refiners[2]->get_name(), "RotamerTrialsRefiner");
		TS_ASSERT_EQUALS(refiners[3]->get_name(), "MinimizationRefiner");

		utilities::PeriodicMoverOP periodic_refiner =
			utility::pointer::dynamic_pointer_cast< utilities::PeriodicMover > ( refiners[1] );

		TS_ASSERT_EQUALS(periodic_refiner->period_, 20);
		TS_ASSERT_EQUALS(periodic_refiner->mover_->get_name(), "RepackingRefiner");

	}

	static void assert_default_fullatom_stage(LoopModelerOP modeler) { // {{{1
		assert_default_fullatom_params(modeler);
		assert_default_fullatom_movers(modeler);
		assert_default_fullatom_refiners(modeler);
	}

	static void assert_default_params(LoopModelerOP modeler) { // {{{1
		assert_default_build_params(modeler);
		assert_default_centroid_params(modeler);
		assert_default_fullatom_params(modeler);
	}
	// }}}1

	void test_default_config() { // {{{1
		string tag = "<LoopModeler/>";
		LoopModelerOP modeler = parse_tag<LoopModeler>(tag);

		assert_default_build_stage(modeler);
		assert_default_centroid_stage(modeler);
		assert_default_fullatom_stage(modeler);
	}

	void test_empty_config() { // {{{1
		string tag = "<LoopModeler config=empty/>";
		LoopModelerOP modeler = parse_tag<LoopModeler>(tag);

		assert_default_params(modeler);

		TS_ASSERT(get_centroid_movers(modeler).empty());
		TS_ASSERT(get_centroid_refiners(modeler).empty());
		TS_ASSERT(get_fullatom_movers(modeler).empty());
		TS_ASSERT(get_fullatom_refiners(modeler).empty());
	}
	
	void test_modeler_options() { // {{{1
		string tag =
			"<LoopModeler"
			"  loops_file=\"protocols/loop_modeling/inputs/loop.txt\""
			"  auto_refine=no/>";
		LoopModelerOP modeler = parse_tag<LoopModeler>(tag);

		assert_default_build_params(modeler);
		assert_default_centroid_params(modeler);
		assert_default_fullatom_params(modeler);
		assert_default_centroid_movers(modeler);
		assert_default_fullatom_movers(modeler);

		TS_ASSERT_EQUALS(modeler->get_loops().size(), 1);
		TS_ASSERT_EQUALS(modeler->get_loops()[1].start(), 34);
		TS_ASSERT_EQUALS(modeler->get_loops()[1].cut(), 41);
		TS_ASSERT_EQUALS(modeler->get_loops()[1].stop(), 46);
		TS_ASSERT_EQUALS(modeler->get_loops()[1].skip_rate(), 0);
		TS_ASSERT_EQUALS(modeler->get_loops()[1].is_extended(), true);

		TS_ASSERT(get_centroid_refiners(modeler).empty());
		TS_ASSERT(get_fullatom_refiners(modeler).empty());
	}

	void test_build_options() { // {{{1
		string tag = 
			"<LoopModeler>"
			"  <Build skip=yes max_attempts=5000/>"
			"</LoopModeler>";
		LoopModelerOP modeler = parse_tag<LoopModeler>(tag);

		assert_default_centroid_stage(modeler);
		assert_default_fullatom_stage(modeler);

		TS_ASSERT_EQUALS(modeler->is_build_stage_enabled_, false);
		TS_ASSERT_EQUALS(modeler->build_stage_->get_max_attempts(), 5000);
	}

	void test_centroid_options() { // {{{1
		string tag = 
			"<LoopModeler>"
			"  <Centroid"
      "    sfxn_cycles=21 temp_cycles=53 mover_cycles=11"
      "    ramp_rama=yes ramp_rep=yes ramp_temp=no"
			"    initial_temp=2.6 final_temp=0.6"
      "    skip=yes auto_refine=no/>"
			"</LoopModeler>";
		LoopModelerOP modeler = parse_tag<LoopModeler>(tag);

		assert_default_build_stage(modeler);
		assert_default_centroid_movers(modeler);
		assert_default_fullatom_stage(modeler);

		TS_ASSERT_EQUALS(modeler->is_centroid_stage_enabled_, false);
		TS_ASSERT_EQUALS(modeler->centroid_stage_->sfxn_cycles_, 21);
		TS_ASSERT_EQUALS(modeler->centroid_stage_->temp_cycles_, 53);
		TS_ASSERT_EQUALS(modeler->centroid_stage_->scale_temp_cycles_, false);
		TS_ASSERT_EQUALS(modeler->centroid_stage_->mover_cycles_, 11);
		TS_ASSERT_EQUALS(modeler->centroid_stage_->ramp_sfxn_rep_, true);
		TS_ASSERT_EQUALS(modeler->centroid_stage_->ramp_sfxn_rama_, true);
		TS_ASSERT_EQUALS(modeler->centroid_stage_->ramp_temp_, false);
		TS_ASSERT_EQUALS(modeler->centroid_stage_->initial_temp_, 2.6);
		TS_ASSERT_EQUALS(modeler->centroid_stage_->final_temp_, 0.6);

		TS_ASSERT(get_centroid_refiners(modeler).empty())
	}

	void test_fullatom_options() { // {{{1
		string tag = 
			"<LoopModeler>"
			"  <Fullatom"
      "    sfxn_cycles=28 temp_cycles=55 mover_cycles=18"
      "    ramp_rama=yes ramp_rep=yes ramp_temp=no"
			"    initial_temp=2.4 final_temp=0.4"
      "    skip=yes auto_refine=no/>"
			"</LoopModeler>";
		LoopModelerOP modeler = parse_tag<LoopModeler>(tag);

		assert_default_build_stage(modeler);
		assert_default_centroid_stage(modeler);
		assert_default_fullatom_movers(modeler);

		TS_ASSERT_EQUALS(modeler->is_fullatom_stage_enabled_, false);
		TS_ASSERT_EQUALS(modeler->fullatom_stage_->sfxn_cycles_, 28);
		TS_ASSERT_EQUALS(modeler->fullatom_stage_->temp_cycles_, 55);
		TS_ASSERT_EQUALS(modeler->fullatom_stage_->scale_temp_cycles_, false);
		TS_ASSERT_EQUALS(modeler->fullatom_stage_->mover_cycles_, 18);
		TS_ASSERT_EQUALS(modeler->fullatom_stage_->ramp_sfxn_rep_, true);
		TS_ASSERT_EQUALS(modeler->fullatom_stage_->ramp_sfxn_rama_, true);
		TS_ASSERT_EQUALS(modeler->fullatom_stage_->ramp_temp_, false);
		TS_ASSERT_EQUALS(modeler->fullatom_stage_->initial_temp_, 2.4);
		TS_ASSERT_EQUALS(modeler->fullatom_stage_->final_temp_, 0.4);

		TS_ASSERT(get_fullatom_refiners(modeler).empty())
	}

	void test_temp_cycles_regex() { // {{{1
		string tag; LoopModelerOP modeler;

		// Case 1: Number of cycles explicitly given.

		tag = "<LoopModeler> <Centroid temp_cycles=53/> </LoopModeler>";
		modeler = parse_tag<LoopModeler>(tag);

		TS_ASSERT_EQUALS(modeler->centroid_stage_->temp_cycles_, 53);
		TS_ASSERT_EQUALS(modeler->centroid_stage_->scale_temp_cycles_, false);

		// Case 2: Number of cycles multiplied by loop length.

		tag = "<LoopModeler> <Centroid temp_cycles=11x/> </LoopModeler>";
		modeler = parse_tag<LoopModeler>(tag);

		TS_ASSERT_EQUALS(modeler->centroid_stage_->temp_cycles_, 11);
		TS_ASSERT_EQUALS(modeler->centroid_stage_->scale_temp_cycles_, true);

		// Case 3: Illegal value.

		tag = "<LoopModeler> <Centroid temp_cycles=x/> </LoopModeler>";
		TS_ASSERT_THROWS(parse_tag<LoopModeler>(tag), EXCN_Msg_Exception);
	}

	void test_custom_movers() { // {{{1
		string tag = 
			"<LoopModeler>"
			"  <LegacyKicSampler/>"
			"  <KicMover/>"
			"  <RotamerTrialsRefiner/>"
			"  <RepackingRefiner/>"
			"  <MinimizationRefiner/>"
			"  <LoopBuilder/>"
			"  <LoopProtocol/>"
			"</LoopModeler>";
		LoopModelerOP modeler = parse_tag<LoopModeler>(tag);

		assert_default_build_stage(modeler);
		assert_default_centroid_params(modeler);
		assert_default_centroid_refiners(modeler);
		assert_default_fullatom_params(modeler);
		assert_default_fullatom_refiners(modeler);

		LoopMovers const & movers = get_centroid_movers(modeler);

		TS_ASSERT_EQUALS(movers.size(), 7);
		TS_ASSERT_EQUALS(movers[1]->get_name(), "LegacyKicSampler");
		TS_ASSERT_EQUALS(movers[2]->get_name(), "KicMover");
		TS_ASSERT_EQUALS(movers[3]->get_name(), "RotamerTrialsRefiner");
		TS_ASSERT_EQUALS(movers[4]->get_name(), "RepackingRefiner");
		TS_ASSERT_EQUALS(movers[5]->get_name(), "MinimizationRefiner");
		TS_ASSERT_EQUALS(movers[6]->get_name(), "LoopBuilder");
		TS_ASSERT_EQUALS(movers[7]->get_name(), "LoopProtocol");
	}

	void test_illegal_movers() { // {{{1
		string tag = 
			"<LoopModeler>"
			"  <MutateResidue/>"
			"</LoopModeler>";

		// Can only use LoopMover subclasses within LoopModeler, not regular 
		// Movers.

		TS_ASSERT_THROWS(parse_tag<LoopModeler>(tag), EXCN_Msg_Exception);
	}

	void test_custom_centroid_movers() { // {{{1
		// Only test one nested mover here.  The more comprehensive test is 
		// test_custom_movers().
		
		string tag = 
			"<LoopModeler>"
			"  <Centroid>"
			"    <LegacyKicSampler/>"
			"  </Centroid>"
			"</LoopModeler>";
		LoopModelerOP modeler = parse_tag<LoopModeler>(tag);

		assert_default_build_stage(modeler);
		assert_default_centroid_params(modeler);
		assert_default_centroid_refiners(modeler);
		assert_default_fullatom_stage(modeler);

		LoopMovers const & movers = get_centroid_movers(modeler);

		TS_ASSERT_EQUALS(movers.size(), 1);
		TS_ASSERT_EQUALS(movers[1]->get_name(), "LegacyKicSampler");
	}

	void test_custom_fullatom_movers() { // {{{1
		// Only test one nested mover here.  The more comprehensive test is 
		// test_custom_movers().
		
		string tag = 
			"<LoopModeler>"
			"  <Fullatom>"
			"    <LegacyKicSampler/>"
			"  </Fullatom>"
			"</LoopModeler>";
		LoopModelerOP modeler = parse_tag<LoopModeler>(tag);

		assert_default_build_stage(modeler);
		assert_default_centroid_stage(modeler);
		assert_default_fullatom_params(modeler);
		assert_default_fullatom_refiners(modeler);

		LoopMovers const & movers = get_fullatom_movers(modeler);

		TS_ASSERT_EQUALS(movers.size(), 1);
		TS_ASSERT_EQUALS(movers[1]->get_name(), "LegacyKicSampler");
	}
	// }}}1

};

