// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/public/rosetta_scripts/rosetta_scripts.cc
/// @brief The application file for rosetta_scripts, aka jd2_scripting or the parser
/// @author Sarel Fleishman (sarelf@u.washington.edu)

// Project Headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/InnerJob.hh>
#include <protocols/jd2/MPIWorkPoolJobDistributor.hh>
#include <protocols/jd2/MPIFileBufJobDistributor.hh>
#include <protocols/jd2/BOINCJobDistributor.hh>
#include <protocols/viewer/viewers.hh>
#include <core/types.hh>
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <core/pose/Pose.hh>
#include <protocols/rosetta_scripts/RosettaScriptsParser.hh>

#include <devel/init.hh>
#include <basic/options/option.hh>

// Utility Headers

// Unit Headers
#include <protocols/moves/Mover.fwd.hh>

#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>

#include <utility/vector1.hh>
#include <utility/excn/EXCN_Base.hh>
#include <basic/Tracer.hh>

// Tracer
static THREAD_LOCAL basic::Tracer TR( "apps.public.rosetta_scripts.rosetta_scripts" );

// FUNCTION PROTOTYPES
void* my_main( void *);

// FUNCTION DECLARATIONS
void*
my_main( void *)
{
	protocols::moves::MoverOP mover;//note that this is not instantiated and will crash if the job distributor actually tries to use it. That means that this can only be used with parser=true
	protocols::jd2::JobDistributor::get_instance()->go(mover);
	return 0 ;
}

/// @details dock_design_scripting provides an xml-based scripting capability
/// to run rosetta movers and filters defined in a text file provided by the
/// user. A full documentation of dock_design_scripting is available at:
/// manual_doxygen/applications/app_dock_design.dox
int
main( int argc, char * argv [] )
{
	try{
		protocols::abinitio::ClassicAbinitio::register_options();
		// setup random numbers and options
		devel::init(argc, argv);
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		using namespace protocols::jd2;
		using namespace core::pose;

		protocols::moves::MoverOP mover;//note that this is not instantiated and will crash if the job distributor actually tries to use it.

		protocols::rosetta_scripts::RosettaScriptsParser rs;
		rs.generate_mover_from_pose(
			JobCOP( new Job( InnerJobOP( new InnerJob( "foo", 1 ) ), 1 ) ),
			*PoseOP( new Pose),
			mover,
			true,
			option[ parser::protocol ].value(),
			true
		);

	} catch( utility::excn::EXCN_Base& excn ) {
		basic::Error()
			<< "ERROR: Exception caught by rosetta_scripts application:"
			<< excn << std::endl;
		std::exit( 1 );
	}
}

