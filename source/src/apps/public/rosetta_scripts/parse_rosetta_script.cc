// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/public/rosetta_scripts/parse_rosetta_script.cc
/// @brief  Using the same command line arguments as the rosetta_scripts application,
///         create all of the objects defined within a RosettaScript using the parse-my-tag
///         and parse-tag functions of those objects -- but do not call apply on the resulting
///         ParsedProtocolMover that is produced.
/// @author Andy Watkins (andy.watkins2@gmail.com)

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
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>

// Tracer
static basic::Tracer TR( "apps.public.rosetta_scripts.parse_rosetta_script" );

/// @details Validate an input XML script against the internally-generated XSD and then construct
/// all of the Movers/Filters/etc. that are defined within that script. Exits with a 0 exit status
/// if it succeeds and an exit status of 1 if not.
int
main( int argc, char * argv [] )
{
	try{
		protocols::abinitio::ClassicAbinitio::register_options();
		// setup random numbers and options
		devel::init(argc, argv);
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		using namespace core::pose;

		protocols::rosetta_scripts::RosettaScriptsParser rs;
		rs.generate_mover_and_apply_to_pose( *PoseOP( new Pose), option[ parser::protocol ].value());
		TR << "Successfully constructed and initialized all objects specified in the " << option[ parser::protocol ]() << " script" << std::endl;

	} catch (utility::excn::Exception& excn ) {
		basic::Error()
			<< "ERROR: Exception caught by parse_rosetta_script application:"
			<< excn << std::endl;
		std::exit( 1 );
	}
}

