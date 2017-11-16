// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <core/io/pdb/pdb_writer.hh>
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
#include <core/scoring/packing/compute_holes_score.hh>
#include <core/scoring/packing/HolesParams.hh>
#include <core/scoring/packstat/compute_sasa.hh>
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
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>
#include <protocols/relax/FastRelax.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>

static basic::Tracer TR( "relax_building_block" );

using std::string;
using ObjexxFCL::string_of;
using namespace core;
using utility::vector1;
using core::id::AtomID;
using namespace core::scoring::packing;
using namespace ObjexxFCL::format;
using core::scoring::ScoreFunctionOP;
using core::pose::Pose;
typedef numeric::xyzVector<Real> Vec;
typedef numeric::xyzMatrix<Real> Mat;
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
	core::io::silent::SilentFileData sfd;

	// Get a random number tag for the tmp dir for sc calculations
	std::string sctag = string_of(numeric::random::uniform()).substr(2,4);

	// Create a score function object, get the fa_rep weight from it (default = 0.44)
	ScoreFunctionOP sf = get_score_function();

	utility::vector1<std::string> files = option[in::file::s]();
	for ( Size ifile = 1; ifile <= files.size(); ++ifile ) {
		std::string file = files[ifile];

		// Read in pose
		Pose pose;
		import_pose::pose_from_file(pose, file, resi_set, core::import_pose::PDB_file);
		Pose mono = pose;

		// Parse the input filename so that the output filenames can be constructed
		/*
		std::vector<std::string> path_fn_vector = string_split(string_of(file), '/');
		std::vector<std::string> fn_vector = string_split(path_fn_vector[(path_fn_vector.size()-1)], '_');
		Real input_trans = (Real) atoi(fn_vector[3].c_str()); // THIS REALLY NEEDS TO BE FIXED TO HANDLE NON-INT VALUES.
		Real input_angle = (Real) atoi(fn_vector[4].c_str()); // THIS REALLY NEEDS TO BE FIXED TO HANDLE NON-INT VALUES.
		//TR << "Symm: " << fn_vector[0] << " bb: " << fn_vector[1] << " pdb: " << fn_vector[2] << " radial_disp: " << fn_vector[3].c_str() << " angle:  " << fn_vector[4]  << std::endl;
		//TR << input_trans << " " << input_angle << std::endl;
		*/

		// Handle all of the symmetry stuff
		core::pose::symmetry::make_symmetric_pose(pose);
		SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
		std::map<Size,SymDof> dofs = sym_info->get_dofs();
		int sym_jump = 0;
		for ( std::map<Size,SymDof>::iterator i = dofs.begin(); i != dofs.end(); i++ ) {
			Size jump_num = i->first;
			if ( sym_jump == 0 ) {
				sym_jump = jump_num;
			} else {
				utility_exit_with_message("Can only handle one subunit!");
			}
		}
		if ( sym_jump == 0 ) {
			utility_exit_with_message("No jump defined!");
		}
		//Vec start_trans = pose.jump(sym_jump).get_translation();
		//Mat start_rot   = pose.jump(sym_jump).get_rotation();

		Pose const pose_for_relax = pose;

		// Fast relax that sucker
		protocols::relax::FastRelax frelax(sf);
		TR << "Relaxing now" << std::endl;
		frelax.apply(pose_for_relax);
		TR << "Finished relaxing" << std::endl;

		/*
		std::string tag = string_of(numeric::random::uniform()).substr(2,4);
		std::string fn = string_of(fn_vector[0])+"_"+string_of(fn_vector[1])+"_"+string_of(fn_vector[2])+"_"+string_of(input_trans+trans.x())+"_"+string_of(input_angle+delta_ang)+"_"+tag+"_final.pdb.gz";
		*/

		// Rescore with scorefxn
		//ScoreFunctionOP scorefxn = get_score_function();
		//scorefxn->score(pose_for_design);

		// Write the pdb file of the design
		utility::io::ozstream out( option[out::file::o]() + "/" + fn );
		pose_for_relax.dump_pdb(out);
		core::io::pdb::extract_scores(pose_for_relax,out);
		out.close();

	} // ifile

	return NULL;

}

int
main (int argc, char *argv[])
{
	try{
		devel::init(argc,argv);

		void* (*func)(void*) = &dostuff;

		if ( basic::options::option[ basic::options::OptionKeys::parser::view ]() ) {
			protocols::viewer::viewer_main( func );
		} else {
			func(NULL);
		}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}


