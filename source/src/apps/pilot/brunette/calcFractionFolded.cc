// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.AlignmentCluster

/// @file   apps/pilot/brunette/dumpPhiPsiStats
///
/// @brief  dumps pos abego dssp phi psi omega for each protein called on

/// @usage: -in:file:s <protein name> auto dumps output to <protein name>.stats

/// @author TJ Brunette


// Utility Headers
#include <basic/Tracer.hh>

// Core Headers
#include <core/chemical/util.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Residue.hh>

#include <core/id/NamedAtomID.hh>
#include <core/id/SequenceMapping.hh>

#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/scoring/rms_util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/util.hh>

#include <core/types.hh>

#include <core/util/ABEGOManager.hh>

#include <devel/init.hh>

#include <utility/vector1.hh>

 #include <protocols/jumping/util.hh>
//basic & utility
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <utility/io/ozstream.hh>
#include <iostream>
#include <string>
#include <ObjexxFCL/format.hh>

using namespace ObjexxFCL::fmt;
using utility::vector1;
using core::Size;
using core::Real;
using core::pose::Pose;
using namespace core::sequence;

basic::Tracer tr( "calcFractionFolded" );

Real two_region_rmsd(Pose & mod_pose, Pose const & ref_pose, SequenceAlignment aln, Size start_res, Size end_res, Size region2_start_res, Size region2_end_res){
	using namespace core::sequence;
	using namespace core::scoring;
	using namespace core::id;
	vector1<Size> mod_pose_positions;
	vector1<Size> ref_pose_positions;
	SequenceMapping map = aln.sequence_mapping(2,1);
	for (Size ii=start_res; ii<= end_res; ++ii){
		mod_pose_positions.push_back(ii);
		ref_pose_positions.push_back(map[ii]);
	}
	if(region2_start_res != 0){
		for(Size ii=region2_start_res; ii<=region2_end_res; ++ii){
			mod_pose_positions.push_back(ii);
			ref_pose_positions.push_back(map[ii]);
		}
	}
	core::kinematics::FoldTree f_mod(mod_pose_positions.size());
	core::kinematics::FoldTree f_ref(mod_pose_positions.size());
	Pose sub_mod_pose;
	Pose sub_ref_pose;
	core::pose::create_subpose(mod_pose,mod_pose_positions,f_mod,sub_mod_pose);
	core::pose::create_subpose(ref_pose,ref_pose_positions,f_ref,sub_ref_pose);
	return(CA_rmsd(sub_mod_pose,sub_ref_pose));
}

Real region_rmsd(Pose & mod_pose, Pose const & ref_pose, SequenceAlignment aln, Size start_res, Size end_res){
    return(two_region_rmsd(mod_pose,ref_pose,aln,start_res,end_res,0,0));
}
int main( int argc, char * argv [] ) {
    try {
    using namespace core::sequence;
	using namespace core::chemical;
	using namespace core::import_pose::pose_stream;
	using core::import_pose::pose_from_pdb;
	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
    devel::init(argc, argv);
    Real threshold = 0.5; 
	ResidueTypeSetCAP rsd_set = rsd_set_from_cmd_line();
    Pose native_pose;
    pose_from_pdb(
			native_pose,
			*rsd_set,
			option[ in::file::native ]()
			);
    SequenceOP native_sequence = new Sequence(native_pose);
	//create vector of input poses.
	MetaPoseInputStream input = streams_from_cmd_line();
    vector1 < vector1 < Real> > rmsd_by_position(native_pose.total_residue()-8); 
 	while(input.has_another_pose()){
		core::pose::PoseOP input_poseOP;
		input_poseOP = new core::pose::Pose();
		input.fill_pose(*input_poseOP,*rsd_set);
        SequenceOP inputPose_sequence = new Sequence(*input_poseOP);
        SequenceAlignment native_alignment =  align_naive(native_sequence,inputPose_sequence);
        for(Size ii=1; ii<=(Size)input_poseOP->total_residue()-8; ++ii){
            Real rmsd =  region_rmsd(*input_poseOP,native_pose,native_alignment,ii,ii+8);
            rmsd_by_position[ii].push_back(rmsd);
        }
    }
    std::string out_name = option[ out::file::o ];
	utility::io::ozstream output(out_name+".foldingStats");
    output << "pos avgRmsd pctUnder0.5" << std::endl;
    for(Size ii = 1; ii<= (Size)rmsd_by_position.size(); ++ii){
        Real ctBelowThreshold = 0;
        Real totalRmsd = 0;
        for(Size kk =1; kk <= (Size)rmsd_by_position[ii].size(); ++ kk){
            //std::cout << ii << " " << rmsd_by_position[ii][kk] << std::endl;
            totalRmsd += rmsd_by_position[ii][kk];
            if(rmsd_by_position[ii][kk] <= threshold)
                ctBelowThreshold += 1.0;
        }
        Real avgRmsd = totalRmsd/(Real)rmsd_by_position[ii].size();
        Real pctBelowThresh = ctBelowThreshold/(Real)rmsd_by_position[ii].size();
        output << I(4,ii) << F(8,2,avgRmsd) << F(8,2,pctBelowThresh) << std::endl;
    }
    } catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
    }
	return 0;
}
