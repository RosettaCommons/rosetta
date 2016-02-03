// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Grant


// Unit headers
#include <devel/init.hh>
#include <devel/denovo_protein_design/util.hh>

//project Headers
#include <core/io/pdb/pdb_writer.hh>
#include <basic/options/util.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <basic/Tracer.hh>
#include <core/scoring/EnergyMap.hh>

#include <protocols/jobdist/standard_mains.hh>


#include <basic/options/option.hh>
// Utility Headers

// Numeric Headers

// ObjexxFCL Headers

// C++ headers
#include <vector>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <core/chemical/AA.hh>
#include <devel/denovo_protein_design/DesignRelaxMover.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <utility/excn/Exceptions.hh>


// Auto-header: duplicate removed #include <protocols/jobdist/standard_mains.hh>

int
main( int argc, char * argv [] )
{

	try {

  using namespace basic::options;
  using utility::file::FileName;

	devel::init(argc, argv);

	core::pose::Pose pose;

	core::import_pose::pose_from_file( pose, basic::options::start_file() , core::import_pose::PDB_file);

	core::scoring::ScoreFunctionOP fullfxn(core::scoring::get_score_function());

	core::pack::task::TaskFactoryOP designtaskfactory( new core::pack::task::TaskFactory );

	devel::denovo_protein_design::design_setup( pose, designtaskfactory );

	devel::denovo_protein_design::DesignRelaxMover designrelaxmover( designtaskfactory );

	protocols::jobdist::main_plain_pdb_mover( designrelaxmover, fullfxn);

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

  return 0;
}
