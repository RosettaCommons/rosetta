// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


// Project Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>






// Mover headers

#include <core/chemical/ChemicalManager.hh>

// Filter headers
//#include <core/pose/metrics/PoseMetricContainer.fwd.hh>

#include <protocols/pose_metric_calculators/PackstatCalculator.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/util.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <string>
#include <sstream>

#include <core/conformation/Residue.hh> // AUTO IWYU For Pose::Residue Residue::ResidueType
#include <basic/options/keys/OptionKeys.hh> // AUTO IWYU For OptionKeys
#include <utility/file/FileName.fwd.hh> // AUTO IWYU For FileName

//The original author used a lot of using declarations here.  This is a stylistic choice.
// Namespaces
using namespace core;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace pose;
using namespace protocols;
using namespace protocols::pose_metric_calculators;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::id;
using basic::Error;
using basic::Warning;
using utility::file::FileName;

// tracer - used to replace cout
static basic::Tracer TR("test_d_l_readin" );

int
main( int argc, char* argv[] )
{
	try {

		//utility::vector1< core::Size > empty_vector(0);

		// init command line options
		devel::init(argc, argv);

		ResidueTypeSetCOP rts( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD ) );

		Pose pose;
		import_pose::pose_from_file( pose, "inputs/a98.pdb" , core::import_pose::PDB_file);
		TR << "Name of residue 2 is A98: " << pose.residue( 2 ).type().name() << std::endl;

		pose.clear();
		import_pose::pose_from_file( pose, "inputs/d98.pdb" , core::import_pose::PDB_file);
		TR << "Name of residue 2 is A98 (it's not chiral!): " << pose.residue( 2 ).type().name() << std::endl;

		pose.clear();
		import_pose::pose_from_file( pose, "inputs/nlu.pdb" , core::import_pose::PDB_file);
		TR << "Name of residue 2 is NLU: " << pose.residue( 2 ).type().name() << std::endl;

		pose.clear();
		import_pose::pose_from_file( pose, "inputs/dnlu.pdb" , core::import_pose::PDB_file);
		TR << "Name of residue 2 is DNLU: " << pose.residue( 2 ).type().name() << std::endl;


	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}//main

