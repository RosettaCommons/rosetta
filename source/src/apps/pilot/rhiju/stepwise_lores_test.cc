// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file stepwise_lores_test.cc
/// @author Rhiju Das (rhiju@stanford.edu)

// libRosetta headers
#include <core/types.hh>
#include <core/chemical/util.hh>
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/scoring/CachedVdwScreenInfo.hh>
#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/rna/denovo/RNA_FragmentMonteCarlo.hh>
#include <protocols/rna/denovo/options/RNA_FragmentMonteCarloOptions.hh>
#include <protocols/rna/denovo/fragments/FullAtomRNA_Fragments.hh>
#include <protocols/rna/denovo/setup/RNA_DeNovoPoseInitializer.hh>
#include <protocols/rna/denovo/libraries/RNA_JumpLibrary.hh>
#include <protocols/rna/denovo/libraries/RNA_ChunkLibrary.hh>
#include <protocols/rna/denovo/libraries/BasePairStepLibrary.hh>
#include <protocols/toolbox/AtomLevelDomainMap.hh>
#include <protocols/rna/denovo/util.hh>
#include <basic/database/open.hh>

#include <protocols/viewer/viewers.hh>

#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>

//////////////////////////////////////////////////////////
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <list>

using namespace protocols;
using utility::vector1;

static basic::Tracer TR( "apps.pilot.rhiju.stepwise_lores_test" );

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
stepwise_lores_test()
{
	using namespace core::pose;
	using namespace core::scoring;
	using namespace core::chemical;
	using namespace core::pose::full_model_info;
	using namespace protocols::stepwise::modeler;
	using namespace protocols::stepwise::setup;
	using namespace protocols::rna::denovo;
	using namespace protocols::rna::denovo::options;
	using namespace utility::file;

	// Following is 'standard' setup from stepwise.cc
	ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	PoseOP native_pose, align_pose;
	PoseOP pose_op = initialize_pose_and_other_poses_from_command_line( rsd_set );
	protocols::scoring::fill_vdw_cached_rep_screen_info_from_command_line( *pose_op );
	initialize_native_and_align_pose( native_pose, align_pose, rsd_set, pose_op );
	Pose & pose = *pose_op;

	Vector center_vector = ( align_pose != 0 ) ? get_center_of_mass( *align_pose ) : Vector( 0.0 );
	protocols::viewer::add_conformation_viewer ( pose.conformation(), "current", 500, 500, false, ( align_pose != 0 ), center_vector );

	// Following checks if FARFAR can now be run on that pose, obeying its FullModelInfo.
	RNA_FragmentMonteCarloOptionsOP options( new RNA_FragmentMonteCarloOptions );
	options->initialize_from_command_line();
	RNA_FragmentMonteCarlo rna_fragment_monte_carlo( options );
	rna_fragment_monte_carlo.set_native_pose( native_pose );
	rna_fragment_monte_carlo.set_denovo_scorefxn( ScoreFunctionFactory::create_score_function( "rna/denovo/rna_lores.wts" ) );
	rna_fragment_monte_carlo.set_hires_scorefxn( get_rna_hires_scorefxn() );
	rna_fragment_monte_carlo.set_refine_pose( true ); // no heating, etc.
	rna_fragment_monte_carlo.apply( pose );

	pose.dump_pdb( "output.pdb" ); // could/should switch to silent file.

}

///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	clock_t const my_main_time_start( clock() );
	stepwise_lores_test();
	protocols::viewer::clear_conformation_viewers();
	std::cout << "Total time to run " << static_cast<core::Real>( clock() - my_main_time_start ) / CLOCKS_PER_SEC << " seconds." << std::endl;
	exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		devel::init(argc, argv);

		protocols::viewer::viewer_main( my_main );

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}


