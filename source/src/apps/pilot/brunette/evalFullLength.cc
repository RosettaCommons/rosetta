// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.AlignmentCluster

/// @file   apps/pilot/brunette/evalFullLength
///
/// @brief  analyzes AA composition of Full length protein
//
/// @usage:

/// @author TJ Brunette


// Utility Headers
#include <basic/Tracer.hh>
#include <numeric/conversions.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.hh>

// Core Headers
#include <core/chemical/util.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Residue.hh>

#include <core/id/NamedAtomID.hh>

#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <core/scoring/packing/HolesParams.hh>
#include <core/scoring/packstat/compute_sasa.hh>

#include <core/types.hh>

#include <core/util/ABEGOManager.hh>
#include <core/util/SwitchResidueTypeSet.hh>

#include <devel/init.hh>

#include <utility/vector1.hh>

//protocols
#include <protocols/jumping/util.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

//basic & utility
#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/holes.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/remodel.OptionKeys.gen.hh>
#include <utility/io/ozstream.hh>
#include <iostream>
#include <ObjexxFCL/format.hh>

using namespace ObjexxFCL::format;
using utility::vector1;
using core::Size;
using core::Real;

static thread_local basic::Tracer tr( "evalFullLength" );


Size ala_ct(const core::pose::Pose& pose){
  Size score = 0;
	for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		if(pose.residue(ii).name3() == "ALA")
            score++;
	}
	std::cout << "score: " << score << std::endl;
	return score;
}

Size glu_ct(const core::pose::Pose& pose){
  Size score = 0;
	for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		if(pose.residue(ii).name3() == "GLU")
            score++;
	}
	std::cout << "score: " << score << std::endl;
	return score;
}

Size tyr_ct(const core::pose::Pose& pose){
  Size score = 0;
	for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		if(pose.residue(ii).name3() == "TYR")
            score++;
	}
	std::cout << "score: " << score << std::endl;
	return score;
}


Real get_holes_score(const core::pose::Pose& pose){
	core::scoring::packing::HolesParams hp_resl,hp_dec,hp_dec15;
	hp_dec15.read_data_file(basic::database::full_name("scoring/rosettaholes/decoy15.params"));
	Real  holes_result = core::scoring::packing::compute_dec15_score(pose);
	return holes_result;
}


int main( int argc, char * argv [] ) {
    try {
	using namespace core::chemical;
	using namespace core::import_pose::pose_stream;
	using core::import_pose::pose_from_pdb;
	using protocols::loops::Loops;
	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	devel::init(argc, argv);
	std::string out_nametag = option[ out::file::o ];
	std::string out_file_name_str( out_nametag + ".scores");
	utility::io::ozstream output(out_file_name_str);
	ResidueTypeSetCAP rsd_set = rsd_set_from_cmd_line();
	std::cout << "rsd_name: " << rsd_set->name() << std::endl;
	//create vector of input poses.
	MetaPoseInputStream input = streams_from_cmd_line();
	vector1<core::pose::PoseOP> poses;
	output << "score  holes  alaCt  gluCt tyrCt tag" << std::endl;
 	while(input.has_another_pose()){
		core::pose::PoseOP input_poseOP;
		input_poseOP = new core::pose::Pose();
		input.fill_pose(*input_poseOP,*rsd_set);
		std::string tag = core::pose::tag_from_pose(*input_poseOP);
	    Real holesScore = 0;
        Real score = 0;
        if("fa_standard" == rsd_set->name()){
            holesScore = get_holes_score(*input_poseOP);
            core::scoring::ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function(TALARIS_2013 ));
			score = scorefxn->score(*input_poseOP);
        }else{
            core::scoring::ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function("score3"));
            score = scorefxn->score(*input_poseOP);
        }
        Size alaCt = ala_ct(*input_poseOP);
        Size gluCt = glu_ct(*input_poseOP);
        Size tyrCt = tyr_ct(*input_poseOP);
        output << F(8,3,score) << " " << F(8,3,holesScore) <<" "<< I(4,alaCt) << " "<< I(4,gluCt) << " " << I(4,tyrCt) << " " << tag << " " << std::endl;
        }
    } catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

