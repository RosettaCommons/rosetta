// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.AlignmentCluster

/// @file   apps/pilot/brunette/evalRepeats
///
/// @brief  analyzes helix oritentation in repeat proteins

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
#include <utility/io/ozstream.hh>
#include <iostream>
#include <ObjexxFCL/format.hh>

using namespace ObjexxFCL::format;
using utility::vector1;
using core::Size;
using core::Real;

static thread_local basic::Tracer tr( "evalRepeats" );

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
			helices->add_loop(Loop(startHelix,endHelix));
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
	avg_ca_position(pose, Loop(helix1.start(), helix1.stop()),&a);
	avg_ca_position(pose, Loop(helix2.start(), helix2.stop()),&b);
	return(a.distance(b));
}

Real get_distance_midCA(const core::pose::Pose& pose, const protocols::loops::Loop helix1, const protocols::loops::Loop helix2){
	using protocols::loops::Loop;
	using numeric::xyzVector;
	//std::cout <<"helix length" << helix1.length() << "," << helix2.length() << std::endl;
	if((helix1.length() < 3) || (helix2.length() < 3)){
		tr.Warning << "distance invoked with helix shorter than 3 residues" << std::endl;
		return(-1);
	}
	xyzVector<double> a, b;
	Size midHelix1 = (Size)floor((helix1.stop()-helix1.start())/2.0)+helix1.start();
	Size midHelix2 = (Size)floor((helix2.stop()-helix2.start())/2.0)+helix2.start();
	a = pose.xyz(core::id::NamedAtomID("CA", midHelix1));
	b = pose.xyz(core::id::NamedAtomID("CA", midHelix2));
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

Real get_holes_score(const core::pose::Pose& pose){
	core::scoring::packing::HolesParams hp_resl,hp_dec,hp_dec15;
	hp_dec15.read_data_file(basic::database::full_name("scoring/rosettaholes/decoy15.params"));
	Real  holes_result = core::scoring::packing::compute_dec15_score(pose);
	return holes_result;
}

void get_angles(const core::pose::Pose& pose, const protocols::loops::Loop helix1, const protocols::loops::Loop helix2, Real& theta, Real& sigma, Real& phi){
	using protocols::loops::Loop;
	using numeric::xyzVector;
	using numeric::conversions::degrees;
	if((helix1.length() < 3) || (helix2.length() < 3)){
		tr.Warning << "distance invoked with helix shorter than 3 residues" << std::endl;
	}
	xyzVector<double> h1_start, h1_stop, h2_start, h2_stop;
	avg_ca_position(pose, Loop(helix1.start(), helix1.start()+2),&h1_start); //you want the two helices to start at the same turn.
	avg_ca_position(pose, Loop(helix1.stop()-2, helix1.stop()),&h1_stop);
	avg_ca_position(pose, Loop(helix2.start(), helix2.start()+2),&h2_start);
	avg_ca_position(pose, Loop(helix2.stop()-2, helix2.stop()),&h2_stop);
	numeric::xyzVector<Real> v1(h1_stop - h1_start);
	numeric::xyzVector<Real> v2(h2_stop - h2_start);
	numeric::xyzVector<Real> u1 = v1.normalize();
	numeric::xyzVector<Real> u2 = v2.normalize();
	numeric::xyzVector<Real> r = h1_start-h2_start;
	numeric::xyzVector<Real> ur = r.normalize();
	theta = degrees(numeric::arccos(u1.dot_product(u2)));
	sigma = degrees(numeric::arccos(u1.dot_product(ur)));
	Real phi_top = u2.dot_product((u1.cross_product(ur)));
	Real phi_bottom = u2.dot_product(u1.cross_product(u1.cross_product(ur)));
	phi = degrees(std::atan(phi_top/phi_bottom));
}

Real res_type_score(const core::pose::Pose& pose){
  Real score = 0;
	for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		if((pose.residue(ii).name3() == "VAL") || (pose.residue(ii).name3() == "ILE") || (pose.residue(ii).name3() == "LEU"))
			score++;
		else
			if((pose.residue(ii).name3() == "ALA") || (pose.residue(ii).name3() == "TRP") || (pose.residue(ii).name3() == "MET"))
				score=score-1;
	}
	std::cout << "score: " << score << std::endl;
	return score;
}

Real get_hb_srbb_score(const core::pose::Pose& pose){
	//note hb_srbb flag increases the hbond_sr_bb.
	using namespace core::scoring;
	core::pose::Pose centroid_pose=pose;
	core::util::switch_to_residue_type_set( centroid_pose, core::chemical::CENTROID);
	core::Real score = pose.energies().total_energies()[ScoreType(hbond_sr_bb)];
	std::cout << "testing" << std::endl;
	return(score);
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
	output << "score" <<" " << "holes" << " " << "distance1_4" << " " << "distance2_5" << " "<< "distance3_6" <<" " <<  "distance1_2" << " " << "distance2_3" <<" " << "tag" <<" " <<"resTypeScore" << " " << "hbsrbb" << "  " << "abego" << std::endl;
	while(input.has_another_pose()){
		core::pose::PoseOP input_poseOP;
		input_poseOP = new core::pose::Pose();
		input.fill_pose(*input_poseOP,*rsd_set);
		std::string tag = core::pose::tag_from_pose(*input_poseOP);
		Loops helices;
	//get_helices(*input_poseOP,&helices);
		helices.add_loop(2,11);
    helices.add_loop(13,22);
    helices.add_loop(26,41);
    helices.add_loop(44,53);
    helices.add_loop(55,64);
    helices.add_loop(68,83);
	  Real distance1_4 = get_distance(*input_poseOP,helices[1],helices[4]);
    Real distance2_5 = get_distance(*input_poseOP,helices[2],helices[5]);
    Real distance3_6 = get_distance(*input_poseOP,helices[3],helices[6]);
    Real distance1_2 = get_distance(*input_poseOP,helices[1],helices[2]);
    Real distance2_3 = get_distance(*input_poseOP,helices[2],helices[3]);
		Real holesScore = get_holes_score(*input_poseOP);
		core::scoring::ScoreFunctionOP scorefxn( get_score_function() );
		Real fa_score = scorefxn->score(*input_poseOP);
		Real resTypeScore = res_type_score(*input_poseOP);
		Real hb_srbb_score = get_hb_srbb_score(*input_poseOP);
		utility::vector1< std::string >  abego_vector = core::util::get_abego(*input_poseOP,1);
		output << F(8,3,fa_score) << " " <<F(8,3,holesScore) <<" "<< F(8,3,distance1_4) << " "<< F(8,3,distance2_5) << " "<< F(8,3,distance3_6) << " "<< F(8,3,distance1_2)  <<" " << F(8,3,distance2_3) << " " << I(6,tag) <<" " << I(4,resTypeScore) << " " <<F(8,3,hb_srbb_score) << "  " << std::endl;
  }
  } catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
  }
  std::cout <<"complete" << std::endl;
	return 0;
}

