// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   apps/pilot/brunette/outputLayerDesignBurial
///
/// @brief  outputs in fasta like format >pdb name burial in format C:core
///       B:boundary  and S:surface, note uses the default SASA values.
//
/// @usage:

/// @author TJ Brunette


// Utility Headers
#include <basic/Tracer.hh>

// Core Headers
#include <core/chemical/util.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Residue.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/id/AtomID_Map.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <core/sequence/ABEGOManager.hh>
#include <core/types.hh>
#include <devel/init.hh>
#include <utility/vector1.hh>

//protocols

//basic & utility
#include <utility/io/ozstream.hh>
#include <iostream>
#include <ObjexxFCL/format.hh>

using namespace ObjexxFCL::format;
using utility::vector1;
using core::Size;
using core::Real;
using core::pose::Pose;

static THREAD_LOCAL basic::Tracer tr( "ouputAbego" );

int main( int argc, char * argv [] ) {
	try{
		using namespace core::chemical;
		using namespace core::import_pose::pose_stream;
		using core::import_pose::pose_from_file;
		devel::init(argc, argv);
		MetaPoseInputStream input = streams_from_cmd_line();
    		ResidueTypeSetCOP rsd_set( rsd_set_from_cmd_line() );
		while(input.has_another_pose()){
			core::pose::PoseOP input_poseOP;
			input_poseOP = core::pose::PoseOP( new core::pose::Pose() );
			input.fill_pose(*input_poseOP,*rsd_set);
			std::string tag = core::pose::tag_from_pose(*input_poseOP);
        		std::string outFile = (tag + ".abego");
        		utility::io::ozstream output(outFile);
       		utility::vector1< std::string >  abego_vector = core::sequence::get_abego(*input_poseOP,1);
       		 for (Size ii=1; ii<=abego_vector.size(); ++ii){
				output << abego_vector[ii];
			}
       		output.close();
   		 }
    } catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
    }
    return 0;
}
