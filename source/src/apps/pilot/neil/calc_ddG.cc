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


#include <basic/database/open.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <basic/options/keys/holes.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/id/AtomID_Map.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>
#include <fstream>
#include <iostream>
#include <math.h>
#include <numeric/random/random.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <protocols/docking/util.hh>
#include <protocols/jobdist/not_universal_main.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/viewer/viewers.hh>
#include <sstream>
#include <string>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/mpistream.hh>
#include <utility/string_util.hh>
#include <protocols/protein_interface_design/movers/ddG.hh>
// Added 100922
// Added 101102
// Added 101103
// AUTO-REMOVED #include <basic/prof.hh>
// AUTO-REMOVED #include <core/id/AtomID_Map.Pose.hh>
// AUTO-REMOVED #include <numeric/model_quality/rms.hh>
// AUTO-REMOVED #include <core/scoring/dssp/Dssp.hh>
// AUTO-REMOVED #include <protocols/simple_moves/MinMover.hh>
// AUTO-REMOVED #include <pstream.h>
// AUTO-REMOVED #include <time.h>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
//Auto Headers

static basic::Tracer TR("design_symm");

using std::string;
using ObjexxFCL::string_of;
using namespace core;
using utility::vector1;
using core::id::AtomID;
using namespace ObjexxFCL::fmt;
using core::scoring::ScoreFunctionOP;
using core::pose::Pose;
typedef vector1<Size> Sizes;

void
*dostuff(void*) {
  using namespace core;
  using namespace basic;
  using namespace options;
  using namespace OptionKeys;
  using namespace pose;
  using namespace core::conformation::symmetry;
  using namespace scoring;
  using namespace utility;
  using basic::options::option;
	
	chemical::ResidueTypeSetCAP resi_set = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");

	utility::vector1<std::string> files = option[in::file::s]();
	for(Size ifile = 1; ifile <= files.size(); ++ifile) {
		std::string file = files[ifile];
		
		// Read in pose
		Pose pose;
		import_pose::pose_from_pdb(pose, file, resi_set);
	
		// Handle all of the symmetry stuff
		core::pose::symmetry::make_symmetric_pose(pose);

		// Get scorefxn
		ScoreFunctionOP scorefxn = getScoreFunction();	
		scorefxn->score(pose);

		// Calculate binding energy
		protocols::protein_interface_design::movers::ddG ddG_mover = protocols::protein_interface_design::movers::ddG::ddG(scorefxn, 1, true);
		ddG_mover.calculate(pose);
		Real ddG = ddG_mover.sum_ddG();
		TR << files[ifile] << " ddG = " << ddG << std::endl;
		ddG_mover.report_ddG(TR);

		// Write the pdb file of the design
		/*
		utility::io::ozstream out( files[ifile] + ".ddG" );
		pose.dump_pdb(out);
		out.close();
		*/
			
	} // ifile

	return NULL;

}

int
main (int argc, char *argv[])
{
	try{
	devel::init(argc,argv);

	void* (*func)(void*) = &dostuff;

	if (basic::options::option[ basic::options::OptionKeys::parser::view ]()) {
		protocols::viewer::viewer_main( func );
	} else {
		func(NULL);
	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}
}



