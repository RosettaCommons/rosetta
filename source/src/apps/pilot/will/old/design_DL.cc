// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /src/apps/pilat/will/genmatch.cc
/// @brief ???

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Stub.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/scoring/ImplicitFastClashCheck.hh>
#include <sstream>
#include <utility/io/izstream.hh>
// #include <devel/init.hh>

// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <apps/pilot/will/mynamespaces.ihh>
#include <apps/pilot/will/will_util.ihh>

using core::kinematics::Stub;
using protocols::scoring::ImplicitFastClashCheck;

static THREAD_LOCAL basic::Tracer TR( "gentetra" );
static core::io::silent::SilentFileData sfd;


int main (int argc, char *argv[]) {

	try {


	devel::init(argc,argv);

	core::chemical::ResidueTypeSetCAP rs = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	for(Size ifile = 1; ifile <= basic::options::option[basic::options::OptionKeys::in::file::s]().size(); ++ifile) {
		std::string infile = basic::options::option[basic::options::OptionKeys::in::file::s]()[ifile];
		Pose pose;
		core::import_pose::pose_from_file(pose,infile, core::import_pose::PDB_file);
		infile = utility::file_basename(infile);

		ScoreFunctionOP sf = core::scoring::get_score_function();

		using namespace core::pack::task;
		PackerTaskOP task = TaskFactory::create_packer_task(pose);
		vector1< bool > aas(20,false);
		aas[core::chemical::aa_ala] = true;
		// aas[core::chemical::aa_gly] = true;
		// aas[core::chemical::aa_pro] = true;

		for(Size i = 1; i <= pose.n_residue(); ++i) {
			if(pose.residue(i).name3()=="PRO" || pose.residue(i).name3()=="DPR") {
				task->nonconst_residue_task(i).prevent_repacking();
			} else {
				task->nonconst_residue_task(i).restrict_absent_canonical_aas(aas);
				task->nonconst_residue_task(i).allow_noncanonical_aa("DAL",*rs);
				// task->nonconst_residue_task(i).allow_noncanonical_aa("DPR",*rs);
			}
		}

		protocols::simple_moves::PackRotamersMover repack( sf, task );
		repack.apply(pose);

		pose.dump_pdb(infile+"_DL.pdb");
	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


