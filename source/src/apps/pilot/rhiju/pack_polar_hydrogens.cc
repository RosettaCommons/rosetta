// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file swa_monte_carlo.cc
/// @author Rhiju Das (rhiju@stanford.edu)

// libRosetta headers
#include <core/types.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/stepwise/modeler/polar_hydrogens/util.hh>
#include <protocols/viewer/viewers.hh>

//////////////////////////////////////////////////
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh> // for option[ out::file::silent  ] and etc.
#include <basic/options/keys/in.OptionKeys.gen.hh> // for option[ in::file::tags ] and etc.
#include <basic/options/option_macros.hh>

#include <basic/Tracer.hh>
#include <utility/file/FileName.hh>


// C++ headers
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <list>

using namespace protocols;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using utility::vector1;

static basic::Tracer TR( "apps.pilot.rhiju.pack_polar_hydrogens" );

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
pack_polar_hydrogens()
{
	using namespace core::pose;
	using namespace core::scoring;
	using namespace core::chemical;
	using namespace core::import_pose;
	using namespace utility::file;

	ResidueTypeSetCOP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	// Following could go to a FullModelSetup class.
	// read starting pose(s) from disk
	utility::vector1< std::string > const & input_files = option[ in::file::s ]();
	PoseOP pose_op = get_pdb_and_cleanup( input_files[1], rsd_set );
	Pose & pose = *pose_op;

	protocols::stepwise::modeler::polar_hydrogens::pack_polar_hydrogens( pose );

	std::string const outfile = option[ out::file::o ]();
	if ( outfile.size() > 0 ) pose.dump_pdb( outfile );

}

///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	clock_t const my_main_time_start( clock() );
	pack_polar_hydrogens();
	protocols::viewer::clear_conformation_viewers();
	std::cout << "Total time to run " << static_cast<core::Real>( clock() - my_main_time_start ) / CLOCKS_PER_SEC << " seconds." << std::endl;
	exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

		utility::vector1< core::Size > blank_size_vector;
		devel::init(argc, argv);
		protocols::viewer::viewer_main( my_main );

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}


