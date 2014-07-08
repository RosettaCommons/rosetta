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
#include <core/io/pdb/pose_io.hh>
#include <basic/options/util.hh>
// AUTO-REMOVED #include <core/pose/util.hh>
// AUTO-REMOVED #include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
// AUTO-REMOVED #include <core/pack/pack_rotamers.hh>
#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <basic/Tracer.hh>
// AUTO-REMOVED #include <basic/MetricValue.hh>
// AUTO-REMOVED #include <protocols/loops/loops_main.hh> //for getting ss from dssp
#include <core/scoring/EnergyMap.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionInfo.hh>
// AUTO-REMOVED #include <protocols/evaluation/RmsdEvaluator.hh>

#include <protocols/jobdist/standard_mains.hh>
// AUTO-REMOVED #include <protocols/simple_moves/PackRotamersMover.hh>
// AUTO-REMOVED #include <protocols/relax_protocols.hh>


// AUTO-REMOVED #include <core/scoring/TenANeighborGraph.hh>
#include <basic/options/option.hh>
// Utility Headers
// AUTO-REMOVED #include <utility/file/file_sys_util.hh>

// Numeric Headers

// ObjexxFCL Headers

// C++ headers
#include <vector>
#include <iostream>
#include <iomanip>
// AUTO-REMOVED #include <fstream>
#include <string>
#include <sstream>
// AUTO-REMOVED #include <utility/assert.hh> //ASSERT_ONLY makes release build happy
#include <core/chemical/AA.hh>
// AUTO-REMOVED #include <devel/DenovoProteinDesign/CreateStartingStructureMover.hh>
#include <devel/denovo_protein_design/DesignRelaxMover.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <utility/excn/Exceptions.hh>


// Auto-header: duplicate removed #include <protocols/jobdist/standard_mains.hh>
// AUTO-REMOVED #include <protocols/moves/MoverContainer.hh>

int
main( int argc, char * argv [] )
{

	try {

  using namespace basic::options;
  using utility::file::FileName;

	devel::init(argc, argv);

	core::pose::Pose pose;

	core::import_pose::pose_from_pdb( pose, basic::options::start_file() );

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
