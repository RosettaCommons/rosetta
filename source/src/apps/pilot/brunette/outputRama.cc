// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/types.hh>

#include <devel/init.hh>
#include <utility/vector1.hh>

//protocols

//basic & utility
#include <utility/io/ozstream.hh>
#include <iostream>
#include <ObjexxFCL/format.hh>

using namespace ObjexxFCL::fmt;
using namespace core;
using namespace core::scoring;
using utility::vector1;
using core::Size;
using core::Real;
using core::pose::Pose;

static thread_local basic::Tracer tr( "ouputrama" );

vector1< Real> calc_rama(Pose & pose){
    core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
    scorefxn->set_weight( core::scoring::rama, 1.0 );
    (*scorefxn)( pose );
    vector1<Real> rama_v;
    for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
	    Real rama( pose.energies().residue_total_energies( ii )[ core::scoring::rama ] );
        rama_v.push_back(rama);
    }
    return(rama_v);
}

vector1< Real> calc_rama2b(Pose & pose){
    core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
    scorefxn->set_weight( core::scoring::rama, 0.0 );
    scorefxn->set_weight( core::scoring::rama2b, 1.0 );
    (*scorefxn)( pose );
    vector1<Real> rama2b_v;
    for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
        Real rama2b( pose.energies().residue_total_energies( ii )[ core::scoring::rama2b ] );
        rama2b_v.push_back(rama2b);
    }
    return(rama2b_v);
}


int main( int argc, char * argv [] ) {
	try{
    using namespace core::chemical;
	using namespace core::import_pose::pose_stream;
	using core::import_pose::pose_from_pdb;
	devel::init(argc, argv);
	MetaPoseInputStream input = streams_from_cmd_line();
    ResidueTypeSetCAP rsd_set = rsd_set_from_cmd_line();
      while(input.has_another_pose()){
		core::pose::PoseOP input_poseOP;
		input_poseOP = new core::pose::Pose();
		input.fill_pose(*input_poseOP,*rsd_set);
		std::string tag = core::pose::tag_from_pose(*input_poseOP);
        std::string outFile = (tag + ".rama");
        utility::io::ozstream output(outFile);
        output << ">" << tag << std::endl;
        vector1<Real> rama_v =  calc_rama(*input_poseOP);
        for(Size ii=1; ii<=rama_v.size(); ++ii)
            output << ii << " " << rama_v[ii] << std::endl;
        output.close();
        std::string outFile2 = (tag + ".rama2b");
        utility::io::ozstream output2(outFile2);
        output2 << ">" << tag << std::endl;
        vector1<Real> rama2b_v =  calc_rama2b(*input_poseOP);
        for(Size ii=1; ii<=rama2b_v.size(); ++ii)
            output2 << ii << " " << rama2b_v[ii] << std::endl;
        output2.close();
        std::cout <<"finished " << tag << std::endl;
    }
    } catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
    }
    return 0;
}
