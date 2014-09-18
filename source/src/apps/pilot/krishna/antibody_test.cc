// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/apps/pilot/krishna/antibody_test.cc
/// @brief
/// @author krishna (kkpraneeth@jhu.edu)
/// 02/16/2013



#include <protocols/antibody2/AntibodyInfo.hh>
#include <protocols/antibody2/AntibodyUtil.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <protocols/rigid/RigidBodyMover.hh>

#include <devel/init.hh>
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>
//#include <utility/tools/make_vector1.hh>

// option key includes
#include <basic/options/option.hh>
#include <string>
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/file/FileName.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/io/pdb/pose_io.hh>


static thread_local basic::Tracer TR( "protocols.antibody2" );

int
main( int argc, char * argv [] )
{
	using namespace basic::options;
	using namespace protocols::antibody2;
	//using namespace protocols::jd2;

	utility::file::FileName pdb_name;
	utility::file::FileName native_pdb_name;
	vector1<core::Size> align_residues, rmsd_residues;
	//vector1<core::Size> temp_residues;

	//protocols::jd2::register_options();
	// initialize core
	devel::init(argc, argv);

	if (option[ OptionKeys::in::file::s ].user()){
		pdb_name = option[ OptionKeys::in::file::s ]()[1];
	}

	if (option[ OptionKeys::in::file::native ].user()){
		native_pdb_name = option[ OptionKeys::in::file::native ]();
	}

	core::pose::Pose pose;
	core::pose::Pose native_pose;
	core::import_pose::pose_from_pdb(pose, pdb_name);
	core::import_pose::pose_from_pdb(native_pose, native_pdb_name);

/*	Size align_residue_list [] = {10,20,30,40};
	Size rmsd_residue_list [] = {50,60,70,80};
	for (Size i=0; i<=1; i++){
		temp_residues.assign(align_residue_list+2*i,align_residue_list+2*(i+1));
		align_residues.push_back(temp_residues);
		temp_residues.clear();
	}
	for (Size j=0; j<=1; j++){
		temp_residues.assign(rmsd_residue_list+2*j,rmsd_residue_list+2*(j+1));
		rmsd_residues.push_back(temp_residues);
		temp_residues.clear();
	}*/

	align_residues.push_back(10);
	align_residues.push_back(20);
	rmsd_residues.push_back(50);
	rmsd_residues.push_back(60);

	core::Real rmsd = align_segment_with_native_and_calc_rmsd_over_req_res(pose,native_pose,align_residues,rmsd_residues);
	TR << "Alignment RMSD is " << rmsd << std::endl;

//	AntibodyInfoOP ab_info = new AntibodyInfo();
//	core::kinematics::FoldTreeOP my_foldtree ( ab_info->LA_H_foldtree(pose) );
//	my_foldtree->show(std::cout);
//	pose.fold_tree(*my_foldtree);
//
//	protocols::rigid::RigidBodyTransMoverOP translate_away ( new protocols::rigid::RigidBodyTransMover( pose, 1 ) );
//	translate_away->step_size( 100 );
//	translate_away->apply( pose );
//	core::io::pdb::dump_pdb( pose, "separated_pose.pdb" );
//	exit(-1);
}



