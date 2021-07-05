// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/vmullig/test_trRosettaProtocolMover.cc
/// @brief A unit test (albeit one in the integration test suite) for the
/// trRosettaProtocolMover.  This mainly tests that the native pose is set
/// properly and retained when the mover is cloned.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// devel headers
#include <devel/init.hh>

// core headers

// protocol headers

// utility headers
#include <utility/excn/Exceptions.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/trRosetta.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>

#ifdef USE_TENSORFLOW
#include <protocols/trRosetta_protocols/movers/trRosettaProtocolMover.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/rms_util.hh>

#include <basic/citation_manager/CitationManager.hh>
#else
#include <basic/tensorflow_manager/util.hh>
#endif


static basic::Tracer TR("apps.pilot.vmullig.test_trRosettaProtocolMover");

#define FLOAT_COMPARISON_THRESHOLD 0.002

/// @brief Indicate which commandline flags are relevant to this application.
void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	option.add_relevant( basic::options::OptionKeys::in::file::fasta );
	option.add_relevant( basic::options::OptionKeys::in::file::native );
	option.add_relevant( basic::options::OptionKeys::trRosetta::msa_file );

}

#ifdef USE_TENSORFLOW

bool
do_test(
	std::string const & native_file,
	std::string const & fasta_file,
	std::string const & msa_file
) {
	using namespace protocols::trRosetta_protocols::movers;

	bool success(true);

	trRosettaProtocolMoverOP tr_mover( utility::pointer::make_shared< trRosettaProtocolMover >() );
	tr_mover->set_fasta_file(fasta_file);
	tr_mover->set_msa_file(msa_file);
	core::pose::PoseOP natpose( core::import_pose::pose_from_file( native_file ) );
	tr_mover->set_native_pose(natpose);

	{
		basic::citation_manager::CitationCollectionList citations;
		tr_mover->provide_citation_info(citations);
		basic::citation_manager::CitationManager::get_instance()->add_citations(citations);
	}

	if( tr_mover->get_native_pose().get() != natpose.get() ) {
		TR.Error << "Initial setting of native pose failed!" << std::endl;
		success = false;
	}

	{
		core::pose::PoseOP pose( utility::pointer::make_shared< core::pose::Pose >() ); //Empty pose.
		TR << "Predicting the structure of ubiqutin..." << std::endl;
		tr_mover->apply(*pose); //Predict structure of ubiqutin.
		core::Real const gdt( core::scoring::native_CA_gdtmm( *natpose, *pose ) );
		core::Real const rmsd( core::scoring::native_CA_rmsd( *natpose, *pose ) );
		TR << "GDT = " << gdt << std::endl;
		if( gdt < 0.5 ) {
			TR.Error << "GDT is under 0.5!" << std::endl;
			success = false;
		}
		TR << "RMSD = " << rmsd << " A" << std::endl;
		if( rmsd > 3.5 ) {
			TR.Error << "RMSD is over 3.5 A!" << std::endl;
			success = false;
		}
		pose->dump_pdb( "OUTPUT.pdb" );
	}

	if( tr_mover->get_native_pose().get() != natpose.get() ) {
		TR.Error << "Native pose not preserved!" << std::endl;
		success = false;
	}

	TR << "Testing clone..." << std::endl;
	protocols::moves::MoverOP tr_mover_clone( tr_mover->clone() );
	if( tr_mover_clone->get_native_pose().get() != natpose.get() ) {
		TR.Error << "Native pose not preserved on clone!" << std::endl;
		success = false;
	}
	if( utility::pointer::dynamic_pointer_cast< protocols::trRosetta_protocols::movers::trRosettaProtocolMover >(tr_mover_clone) == nullptr ) {
		TR.Error << "Clone did not preserve type!" << std::endl;
		success = false;
	}

	return success;
}

#endif //USE_TENSORFLOW

/// @brief Program entry point.
int
main(
	int argc,
	char * argv []
) {
	try {
#ifdef USE_TENSORFLOW
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		register_options();
		devel::init( argc, argv );

		TR << "Starting test_trRosettaProtocolMover." << std::endl;
		TR << "This is a unit test (albeit one in the integration test suite) for the "
			"trRosettaProtocolMover.  This mainly tests that the native pose is set "
			"properly and retained when the mover is cloned." << std::endl;
		TR << "Test created 21 March 2021 by Vikram K. Mulligan, Flatiron "
			"Institute (vmulligan@flatironinstitute.org)." << std::endl;

		runtime_assert_string_msg(
			option[ basic::options::OptionKeys::trRosetta::msa_file ].user() && !option[ basic::options::OptionKeys::trRosetta::msa_file ]().empty(),
			"A multiple sequence alignment must be provided with the -trRosetta::msa_file flag."
		);
		runtime_assert_string_msg(
			option[ in::file::fasta ].user() && (option[ in::file::fasta ]().size() == 1),
				"A FASTA file must be provided with the -in:file:fasta flag."
		);
		runtime_assert_string_msg(
			option[ in::file::native ].user() && !(option[ in::file::native ]().empty()),
				"A native file must be provided with the -in:file:native flag."
		);

		bool const success( do_test( option[ in::file::native ](), option[in::file::fasta]()[1], option[trRosetta::msa_file]() ) );
		TR << "Test " << (success ? "PASSED" : "FAILED" ) << "!" << std::endl;

		basic::citation_manager::CitationManager::get_instance()->write_all_citations_and_unpublished_author_info();

		if( !success ) {
			utility_exit(); //Exit with failure status if we don't succeed.
		}

#else
		devel::init( argc, argv );
		utility_exit_with_message(
			"The test_trRosettaProtocolMover must be compiled with -extras=tensorflow or -extras=tensorflow_gpu!\n\n"
			+ basic::tensorflow_manager::get_tensorflow_compilation_instructions( "test_trRosettaProtocolMover application" )
		);
#endif

	} catch ( utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;
}
