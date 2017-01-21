// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   struct_frag.cc
/// @brief  Creating fragments from supplied structures.
/// @author andreas scheck (graziano.bud@gmail.com), Correia's LPDI/EPFL

#include <protocols/jd2/JobDistributor.hh>

#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <utility/io/ozstream.hh>

#include <protocols/struct_fragment/StructFragmentMover.hh>
// #include <protocols/fold_from_loops/LoopsSelectorManager.hh>
// #include <protocols/fold_from_loops/NubInitio.hh>

// C++ headers
#include <string>

// Rosetta Options
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/ConstDataMap.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>

// #include <basic/options/keys/fold_from_loops.OptionKeys.gen.hh>

// Rosetta Objects
// #include <protocols/moves/MoverContainer.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
// #include <protocols/constraint_generator/AtomPairConstraintGenerator.hh>

static THREAD_LOCAL basic::Tracer TR( "apps.struct_fragment_app" );

// Declaring namespaces
// Alias namespaces are used to keep track of the source of the functions/objects
namespace bo  = basic::options;
namespace oi  = bo::OptionKeys::in::file;
// namespace of  = bo::OptionKeys::fold_from_loops;
namespace ue  = utility::excn;
// namespace ffl = protocols::fold_from_loops;


////////////////////////////////
/// ========= MAIN ========= ///
////////////////////////////////
int
main( int argc, char* argv [] )
{
	try {
		// Initialize the option system
		// register_options();
		devel::init( argc, argv );


		// core::pose::Pose pose;
		// core::import_pose::pose_from_file( pose , "4u0q.pdb" );
		// ffl::LoopsSelectorManager lrmanager = ffl::LoopsSelectorManager();
		// lrmanager.input_loops_file("4u0q.loops");
		// lrmanager.apply(pose);
		// TR << (*lrmanager.get_loops()) << std::endl;

		// Load a core::pose::Pose from the Template's PDB file
		// The Template PDB is provided with -in:file:s option
		// This option is MANDATORY and will rise an exeption if not provided
		// if ( !bo::option[ oi::s ].user() ) {
		// 	std::string msg = "A Template PDB file must be provided through the -in:file:s option";
		// 	throw ue::EXCN_Msg_Exception( msg );
		// }
    core::pose::Pose pose;
		core::import_pose::pose_from_file( pose , "obj01.pdb" );
		// // Prepare Sequence of Events
		// if ( bo::option[ of::native_ca_cst ].user() ) {}
		//
		//
		// // 1. Load the template Pose

		//
		// // 2. Add constraints


		// Initialize a new FoldFromLoopsMover
		protocols::struct_fragment::StructFragmentMoverOP ffl_mover( new protocols::struct_fragment::StructFragmentMover );
		// ffl_mover->parse_command_line();
		ffl_mover->set_small_frag_file("finalTest.10.2mers");
		ffl_mover->set_large_frag_file("finalTest.10.16mers");
		ffl_mover->apply( pose );
		// Execute the FFLMover
		// protocols::jd2::JobDistributor::get_instance()->go( ffl_mover );

		basic::datacache::ConstDataMap testMap = pose.const_data_cache();
		core::fragment::ConstantLengthFragSet small = testMap.get<core::fragment::ConstantLengthFragSet >( "StructFragmentMover", "small");
		core::fragment::ConstantLengthFragSet large = testMap.get<core::fragment::ConstantLengthFragSet >( "StructFragmentMover", "large");


		TR << small.size() << std::endl;
		TR << small.nr_frames() << std::endl;
		TR << large.size() << std::endl;
		TR << large.nr_frames() << std::endl;


		// ffl_mover->set_small_frag_file("customWeightTest.10.2mers");
		// ffl_mover->set_large_frag_file("customWeightTest.10.16mers");
		// ffl_mover->apply( pose );
		// // Execute the FFLMover
		// // protocols::jd2::JobDistributor::get_instance()->go( ffl_mover );
		//
		// testMap = pose.const_data_cache();
		//  small = testMap.get<core::fragment::ConstantLengthFragSet >( "StructFragmentMover", "small");
		//  large = testMap.get<core::fragment::ConstantLengthFragSet >( "StructFragmentMover", "large");
		//
		//
		// TR << small.size() << std::endl;
		// TR << small.nr_frames() << std::endl;
		// TR << large.size() << std::endl;
		// TR << large.nr_frames() << std::endl;

		// if ( testMap.has( "StructFragmentMover" ) ) {
		// 	TR << "category found" << std::endl;
		// }
		//
		// if ( testMap.has( "StructFragmentMover", "large" ) ) {
		// 	TR << "category and name found" << std::endl;
		// 	// core::fragment::ConstantLengthFragSet output = testMap.get<core::fragment::ConstantLengthFragSet>( "StructFragmentMover", "Fragments");
		// 	core::fragment::ConstantLengthFragSet small = testMap.get<core::fragment::ConstantLengthFragSet >( "StructFragmentMover", "small");
		// 	core::fragment::ConstantLengthFragSet large = testMap.get<core::fragment::ConstantLengthFragSet >( "StructFragmentMover", "large");
		//
		//
		// 	TR << small.size() << std::endl;
		// 	TR << small.nr_frames() << std::endl;
		// 	TR << large.size() << std::endl;
		// 	TR << large.nr_frames() << std::endl;
		//
		// }


	} catch ( ue::EXCN_Base const & e ) {
		TR.Fatal << TR.Red   << "[    FATAL ERROR    ]" << TR.Reset << std::endl;
		TR.Fatal << TR.Red   << "[ EXECUTION STOPPED ]" << TR.Reset << std::endl;
		TR.Fatal << TR.Bold  << e.msg() << TR.Reset << std::endl;
		return -1;
	}

	return 0;
}
