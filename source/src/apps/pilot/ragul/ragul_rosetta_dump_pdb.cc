// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @author Ragul Gowthaman

#include <string>

#include <basic/Tracer.hh>
#include <basic/options/util.hh>

#include <core/chemical/AA.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/pose/Pose.hh>




#include <basic/options/keys/OptionKeys.hh> // AUTO IWYU For OptionKeys
#include <utility/excn/Exceptions.hh> // AUTO IWYU For Exception


using namespace core;
using namespace basic::options;
using namespace core::scoring;
using namespace basic::options::OptionKeys;

static basic::Tracer TR( "apps.pilot.dump_pdb.main" );


int main( int argc, char * argv [] )
{
	try{
		devel::init(argc, argv);

		// create native pose from pdb
		pose::Pose pose_init;
		std::string const input_pdb_name( basic::options::start_file() );
		core::import_pose::pose_from_file( pose_init, input_pdb_name , core::import_pose::PDB_file);
		std::string out_pdb_name = "rosetta_" + input_pdb_name;
		pose_init.dump_pdb(out_pdb_name);

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;

}

