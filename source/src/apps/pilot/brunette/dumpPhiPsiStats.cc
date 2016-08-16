// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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

#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/types.hh>

#include <core/sequence/ABEGOManager.hh>

#include <devel/init.hh>

#include <utility/vector1.hh>

 #include <protocols/jumping/util.hh>
//basic & utility
#include <basic/options/option.hh>
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

static THREAD_LOCAL basic::Tracer tr( "dumpPhiPsiStats" );


int main( int argc, char * argv [] ) {
    try {
	using namespace core::chemical;
	using namespace core::import_pose::pose_stream;
	using core::import_pose::pose_from_file;
	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
    devel::init(argc, argv);
	ResidueTypeSetCAP rsd_set = rsd_set_from_cmd_line();
    std::cout << "rsd_name: " << rsd_set->name() << std::endl;
	//create vector of input poses.
	MetaPoseInputStream input = streams_from_cmd_line();
 	while(input.has_another_pose()){
		core::pose::PoseOP input_poseOP;
		input_poseOP = new core::pose::Pose();
		input.fill_pose(*input_poseOP,*rsd_set);
		std::string tag = core::pose::tag_from_pose(*input_poseOP);
        std::string tag_no_pdb = tag.substr(0,tag.size()-4);
		utility::io::ozstream output(tag_no_pdb+".stats");
        output << "pos abego  dssp    phi     psi     omega" << std::endl;
        utility::vector1< std::string >  abego_vector = core::sequence::get_abego(*input_poseOP,1);
         protocols::jumping::assign_ss_dssp( *input_poseOP );
        for(int ii=1; ii<=(int)input_poseOP->total_residue(); ++ii){
            output << I(4,ii)<<"  " <<abego_vector[ii] << "    " << input_poseOP->secstruct(ii) << "   " << F(8,1,input_poseOP->phi(ii)) << F(8,1,input_poseOP->psi(ii)) << F(8,1,input_poseOP->omega(ii)) << std::endl;
        }
    }
    } catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
    }
	return 0;
}
