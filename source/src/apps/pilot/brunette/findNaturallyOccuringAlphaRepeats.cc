// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/brunette/findNaturallyOccuringRepeats
///
/// @brief  looks for repeating helices

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

#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <core/scoring/packing/HolesParams.hh>
#include <core/scoring/packstat/compute_sasa.hh>

#include <core/types.hh>

#include <core/sequence/ABEGOManager.hh>

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
#include <utility/io/ozstream.hh>
#include <iostream>
#include <ObjexxFCL/format.hh>

using namespace ObjexxFCL::format;
using utility::vector1;
using core::Size;
using core::Real;

static THREAD_LOCAL basic::Tracer tr( "alphaRepeatIdenfier" );

void avg_ca_position(
	const core::pose::Pose& pose,
	const protocols::loops::Loop& region,
	numeric::xyzVector<double>* point
) {
	assert(point);

	point->zero();
	for (unsigned i = region.start(); i <= region.stop(); ++i) {
		(*point) += pose.xyz(core::id::NamedAtomID("CA", i));
	}

	(*point) /= region.length();
}

void get_helices(core::pose::Pose& pose, protocols::loops::Loops* helices){
	using protocols::loops::Loop;
	using numeric::xyzVector;
	protocols::jumping::assign_ss_dssp( pose );
	char lastSecStruct = pose.secstruct(1);
	Size startHelix = 0;
	Size endHelix = 0;
	if(pose.secstruct(1) == 'H')
		startHelix = 1;
	for ( core::Size ii = 2; ii <= pose.total_residue(); ++ii ) {
		if(pose.secstruct(ii) == 'H' && lastSecStruct != 'H')
			startHelix = ii;
		if(pose.secstruct(ii) != 'H' && lastSecStruct == 'H'){
			endHelix = ii-1;
			if(endHelix-startHelix >= 2)
				helices->add_loop(Loop(startHelix,endHelix));
		}
		lastSecStruct = pose.secstruct(ii);
	}
}
void get_sheets(core::pose::Pose& pose, protocols::loops::Loops* sheets){
	using protocols::loops::Loop;
	using numeric::xyzVector;
	protocols::jumping::assign_ss_dssp( pose );
	char lastSecStruct = pose.secstruct(1);
	Size startSheet = 0;
	Size endSheet = 0;
	if(pose.secstruct(1) == 'E')
		startSheet = 1;
	for ( core::Size ii = 2; ii <= pose.total_residue(); ++ii ) {
		if(pose.secstruct(ii) == 'E' && lastSecStruct != 'E')
			startSheet = ii;
		if(pose.secstruct(ii) != 'E' && lastSecStruct == 'E'){
			endSheet = ii-1;
			sheets->add_loop(Loop(startSheet,endSheet));
		}
		lastSecStruct = pose.secstruct(ii);
	}
}


Real get_distance(const core::pose::Pose& pose, const protocols::loops::Loop helix1, const protocols::loops::Loop helix2){
	using protocols::loops::Loop;
	using numeric::xyzVector;
	//std::cout <<"helix length" << helix1.length() << "," << helix2.length() << std::endl;
	if((helix1.length() < 3) || (helix2.length() < 3)){
		tr.Warning << "distance invoked with helix shorter than 3 residues" << std::endl;
		return(-1);
	}
	xyzVector<double> a, b;
	avg_ca_position(pose, Loop(helix1.stop()-2, helix1.stop()),&a);
	avg_ca_position(pose, Loop(helix2.start(), helix2.start()+2),&b);
	return(a.distance(b));
}
Real get_distance_endpoint(const core::pose::Pose& pose, const protocols::loops::Loop helix1, const protocols::loops::Loop helix2){
	using protocols::loops::Loop;
	using numeric::xyzVector;
	//std::cout <<"helix length" << helix1.length() << "," << helix2.length() << std::endl;
	if((helix1.length() < 3) || (helix2.length() < 3)){
		tr.Warning << "distance invoked with helix shorter than 3 residues" << std::endl;
		return(-1);
	}
	xyzVector<double> a, b;
	avg_ca_position(pose, Loop(helix2.stop()-2, helix2.stop()),&a);
	avg_ca_position(pose, Loop(helix1.start(), helix1.start()+2),&b);
	return(a.distance(b));
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
	//create vector of input poses.
	MetaPoseInputStream input = streams_from_cmd_line();
	vector1<core::pose::PoseOP> poses;
	output << "tag" << " " <<   "distance" << " " << "distance_endpoint" <<" " << "numb_helices numb_sheets" << std::endl;
	while(input.has_another_pose()){
		core::pose::PoseOP input_poseOP;
		input_poseOP = new core::pose::Pose();
		input.fill_pose(*input_poseOP,*rsd_set);
		std::string tag = core::pose::tag_from_pose(*input_poseOP);
		std::cout << "working on" << tag << std::endl;
		Loops helices;
		Loops sheets;
		get_helices(*input_poseOP,&helices);
		get_sheets(*input_poseOP,&sheets);
		//for(int ii=1; ii<=helices.size(); ++ii)
		//	std::cout << helices[ii].start() <<"," << helices[ii].stop() << std::endl;
		if(helices.num_loop() > 4)
			if((helices[1].length() >= 3)&&(helices[2].length() >= 3)){//assumes repeat proteins with > 1 loops
				Real distance = get_distance(*input_poseOP,helices[1],helices[2]);
				Real distance_endpoint = get_distance_endpoint(*input_poseOP,helices[1],helices[2]);
				output << I(6,tag) <<" "  << F(8,3,distance) << " " << F(8,3,distance_endpoint) << " " << I(4,helices.size()) << " " << I(4,sheets.size())  << std::endl;
			}
	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}

