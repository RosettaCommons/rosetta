// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody_design/antibody_designer.cc
/// @brief Just a quickly written app for basic antibody grafting/modeling for now.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/jd2/JobDistributor.hh>
#include <basic/options/util.hh>
#include <devel/init.hh>
#include <protocols/viewer/viewers.hh>

#include <protocols/antibody/design/AntibodyDesignMoverGenerator.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/grafting/CCDEndsGraftMover.hh>
#include <protocols/grafting/simple_movers/KeepRegionMover.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>

#include <core/types.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/electron_density/util.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>

#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <numeric/constants.hh>

#include <ObjexxFCL/string.functions.hh>
#include <core/kinematics/Jump.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

// C++ headers
#include <iostream>
#include <fstream>
#include <string>

using namespace basic;
using namespace utility;
using namespace protocols;

OPT_1GRP_KEY(String, ab, L1_pose)
OPT_1GRP_KEY(String, ab, L2_pose)
OPT_1GRP_KEY(String, ab, L3_pose)
OPT_1GRP_KEY(String, ab, H1_pose)
OPT_1GRP_KEY(String, ab, H2_pose)
OPT_1GRP_KEY(String, ab, H3_pose)

OPT_1GRP_KEY(String, ab, L1_seq)
OPT_1GRP_KEY(String, ab, L2_seq)
OPT_1GRP_KEY(String, ab, L3_seq)
OPT_1GRP_KEY(String, ab, H1_seq)
OPT_1GRP_KEY(String, ab, H2_seq)
OPT_1GRP_KEY(String, ab, H3_seq)

OPT_1GRP_KEY(String, ab, DE_pose)
OPT_1GRP_KEY(String, ab, DE_seq)

OPT_1GRP_KEY(Boolean, ab, relax_cdrs)
OPT_1GRP_KEY(Boolean, ab, snugdock)

namespace protocols {
namespace antibody {

void
model_cdrs(core::pose::Pose & pose, AntibodyInfoCOP ab_info){

	utility::vector1<std::string> pdb_files(6, "");
	utility::vector1<std::string> sequences(6, "");

	pdb_files[ l1 ] = option [ ab::L1_pose ]();
	pdb_files[ l2 ] = option [ ab::L2_pose ]();
	pdb_files[ l3 ] = option [ ab::L3_pose ]();
	pdb_files[ h1 ] = option [ ab::H1_pose ]();
	pdb_files[ h2 ] = option [ ab::H2_pose ]();
	pdb_files[ h3 ] = option [ ab::H3_pose ]();

	sequences[ l1 ] = option [ ab::L1_seq ]();
	sequences[ l2 ] = option [ ab::L2_seq ]();
	sequences[ l3 ] = option [ ab::L3_seq ]();
	sequences[ h1 ] = option [ ab::H1_seq ]();
	sequences[ h2 ] = option [ ab::H2_seq ]();
	sequences[ h3 ] = option [ ab::H3_seq ]();

	for (core::Size i = 1; i <= 6; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);

		//Graft CDR?
		if (pdbfiles[ i ] != "" ){

			core::pose::PoseOP cdr_pose = pose_from_file(pdbfiles[ i ], core::import_pose::PDB_file);
			keep_cdr(*cdr_pose, *ab_info, cdr)
			graft_in_cdr(pose, *ab_info, cdr, cdr_pose);
		}

		//Copy a new sequence?
		if (sequences[ i ] != "" ){
			copy_in_seq(pose, *ab_info, cdr, sequences[ i ]);

		}


	}

	design::AntibodyDesignMoverGenerator generator = design::AntibodyDesignMoverGenerator(ab_info);
	generator.set_cdr_range(h1, l3, true);
	generator.set_include_neighbor_sc(true);
	generator.set_min_sc(true);
	generator.neighbor_detection_dis(6.0);
	generator.stem_size(2);


	//Repack all CDRs, neighbors, and stem.
	generator.generate_repack_cdrs(pose);
	generator.apply(pose);


	//Add dihedral constraints so we don't screw up the CDRs too much.  Have a bit more give sense they could be wrong for the sequence.
	//Add 10 degrees to phi/psi sd.
	for (core::Size i = 1; i <= 6; ++i ){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		add_harmonic_dihedral_cst_general(ab_info, pose, cdr, 30, 40);
	}



	//Optionally Run Relax on CDRs with dihedral csts
	if (option [ ab::relax_cdrs ] ){
		//Output other pose for before relax - CMD line can control dualspace and nonideal here.
		generator.generate_relax(pose);
		generator.apply(pose);

	}

	//Optionally Relax the whole antibody - Should respect cmd-line settings
	//if (option [ ab::relax_ab ] ) {

	//}

	//Optionally Run SnugDock
	//if (option [ ab::snugdock] ) {
		//output other pose for snugdock
	//}




}


void
keep_cdr(core::pose::Pose & pose, AntibodyInfo & ab_info, CDRNameEnum cdr) {

	using namespace protocols::grafting::simple_movers;
	KeepRegionMover keeper = KeepRegionMover(ab_info.get_CDR_start(cdr, pose), ab_info.get_CDR_end(cdr, pose));
	keeper.apply(pose);

}

void
graft_in_cdr(core::pose::Pose & pose, AntibodyInfo & ab_info, CDRNameEnum cdr, core::pose::PoseOP template_pose){

}

void
copy_in_seq(core::pose::Pose & pose, AntibodyInfo & ab_info, CDRNameEnum cdr, std::string seq){

}





	} //antibody
} //protocols











int main(int argc, char* argv[])
{
	try{


		//NEW_OPT(denstools::lowres, "low res limit", 1000.0);

		NEW_OPT(L1_pose, "L1 PDB file", "");
		NEW_OPT(L1_seq, "L1 sequence", "");

		NEW_OPT(L2_pose, "L2 PDB file", "");
		NEW_OPT(L2_seq, "L2 sequence", "");

		NEW_OPT(L3_pose, "L3 PDB file", "");
		NEW_OPT(L3_seq, "L3 sequence", "");

		NEW_OPT(H1_pose, "H1 PDB file", "");
		NEW_OPT(H1_seq, "H1 sequence", "");

		NEW_OPT(H2_pose, "H2 PDB file", "");
		NEW_OPT(H2_seq, "H2 sequence", "");

		NEW_OPT(H3_pose, "H3 PDB file", "");
		NEW_OPT(H3_seq, "H3 sequence", "");

		OPT(DE_pose);
		OPT(DE_seq);

		NEW_OPT(relax_cdrs, "Relax CDRs and 2 residues into framework with dihedral constraints", true);
		NEW_OPT(snugdock, "Run SnugDock after relax", false);

		devel::init(argc, argv);

		utility::vector1< protocols::jobdist::BasicJobOP > input_jobs = protocols::jobdist::load_s_and_l();
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, input_jobs[ 1 ]->input_tag() , core::import_pose::PDB_file);

		using namespace protocols::antibody;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;


		AntibodyInfoOP ab_info( new AntibodyInfoOP(pose));

		//protocols::jd2::JobDistributor::get_instance()->go(protocols::moves::MoverOP( new AntibodyDesignProtocol ));
	}catch (utility::excn::Exception & excn){
		std::cout << "Exception: "<<std::endl;
		excn.show(std::cerr);
		return -1;
	}

	return(0);

