// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/moves/AnalyzeLoopModeling.fwd.hh
/// @brief connects chains using a very fast RMSD lookback. only works for chains <5 residues. Designed to make loops look within .4 RMSD to naturally occuring loops
/// @author TJ Brunette tjbrunette@gmail.com
///
// Unit headers
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/pose/symmetry/util.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

#include <protocols/pose_length_moves/AnalyzeLoopModeling.hh>
#include <protocols/pose_length_moves/AnalyzeLoopModelingCreator.hh>
#include <protocols/moves/Mover.hh>

#include <core/scoring/dssp/Dssp.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/Loop.hh>

#include <core/indexed_structure_store/SSHashedFragmentStore.hh>
#include <core/indexed_structure_store/FragmentStore.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <protocols/pose_length_moves/NearNativeLoopCloser.hh>

// Core Headers
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <numeric/alignment/QCP_Kernel.hh>
#include <numeric/xyzVector.hh>

#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <ctime>


static THREAD_LOCAL basic::Tracer TR( "protocols.pose_length_moves.AnalyzeLoopModeling" );

namespace protocols {
namespace pose_length_moves {
using namespace core;
using namespace std;
using utility::vector1;



AnalyzeLoopModeling::AnalyzeLoopModeling():moves::Mover("AnalyzeLoopModeling"){}

std::string AnalyzeLoopModelingCreator::keyname() const
{
	return AnalyzeLoopModelingCreator::mover_name();
}

std::string AnalyzeLoopModelingCreator::mover_name(){
	return "AnalyzeLoopModeling";
}


protocols::moves::MoverOP
AnalyzeLoopModelingCreator::create_mover() const {
	return protocols::moves::MoverOP( new AnalyzeLoopModeling );
}

protocols::loops::Loops AnalyzeLoopModeling::get_loops(core::pose::Pose const & pose){
	protocols::loops::Loops pose_loops;
	protocols::loops::Loops pose_loops_pass2;
	core::scoring::dssp::Dssp dssp( pose );
	dssp.dssp_reduced();
	std::string dssp_string = dssp.get_dssp_secstruct();
	string lastSecStruct = dssp_string.substr(0,1);
	Size startLoop = 0;
	Size endLoop = 0;
	Size min_SS_length = 4;
	if ( dssp_string.substr(0,1) == "L" ) {
		startLoop = 1;
	}
	for ( core::Size ii = 2; ii <= pose.total_residue(); ++ii ) {
		if ( dssp_string.substr(ii-1,1) == "L" && lastSecStruct != "L" ) {
			startLoop = ii;
		}
		if ( dssp_string.substr(ii-1,1) != "L" && lastSecStruct == "L" ) {
			endLoop = ii-1;
			if ( (startLoop != 1) && (endLoop!=pose.total_residue()) ) {
				pose_loops.add_loop(startLoop,endLoop);
			}
		}
		lastSecStruct = dssp_string.substr(ii-1,1);
	}
	for ( Size ii=2; ii<pose_loops.num_loop(); ++ii ) {
		Size before_loop_start = pose_loops[ii-1].stop()+1;
		Size before_loop_end = pose_loops[ii].start()-1;
		bool same_res_type = true;
		for(Size jj=before_loop_start+1; jj<=before_loop_end; ++jj) //gets rid of screwing helix-sheet transitions
			if(dssp_string[jj-1]!= dssp_string[before_loop_start-1])
				same_res_type = false;
		Size length_before_loop = before_loop_end-before_loop_start+1;
		Size after_loop_start = pose_loops[ii].stop()+1;
		Size after_loop_end = pose_loops[ii+1].start()-1;
		Size length_after_loop = after_loop_end-after_loop_start+1;
		for(Size jj=after_loop_start+1; jj<=after_loop_end; ++jj)
			if(dssp_string[jj-1]!= dssp_string[after_loop_start-1])
				same_res_type = false;
		if(length_before_loop>=min_SS_length && length_after_loop>=min_SS_length && same_res_type)
			pose_loops_pass2.add_loop(pose_loops[ii].start(),pose_loops[ii].stop());
	}
	//std::cout << dssp_string << std::endl;
	// for (int ii=pose_loops_pass2.num_loop(); ii>=0; --ii ) {
	// 	Size loopStart=pose_loops_pass2[ii].start();
	// 	Size loopEnd=pose_loops_pass2[ii].stop();
	// 	Size loopLength = loopEnd-loopStart+1;
	// 	if(loopLength>=loopLengthRangeLow_ && loopLength<= loopLengthRangeHigh_){
	// 		std::cout << loopStart << "-" << loopEnd << std::endl;
	// 	}
	// }
 	return(pose_loops_pass2);
 }

Real AnalyzeLoopModeling::rmsd_between_coordinates(std::vector< numeric::xyzVector<numeric::Real> > fragCoordinates,std::vector< numeric::xyzVector<numeric::Real> > coordinates){
	numeric::alignment::QCP_Kernel<core::Real>::remove_center_of_mass( &fragCoordinates.front().x() , fragCoordinates.size());
	numeric::alignment::QCP_Kernel<core::Real>::remove_center_of_mass( &coordinates.front().x() , coordinates.size());
	numeric::alignment::QCP_Kernel<core::Real> qcp;
	vector1<Real> rot_vector;
	Real rmsd = qcp.calc_centered_coordinate_rmsd( &coordinates.front().x(), &fragCoordinates.front().x(), fragCoordinates.size(), &rot_vector[1]);
	return(rmsd);
}

Real AnalyzeLoopModeling::get_loop_rmsd(core::pose::Pose native_pose,core::pose::Pose test_pose, Size loopStart, Size loopEnd){
	using namespace chemical;
	std::vector< numeric::xyzVector<numeric::Real> > native_coordinates;
	std::vector< numeric::xyzVector<numeric::Real> > test_coordinates;
	for ( Size ii = loopStart;  ii <= loopEnd; ++ii ) {
			native_coordinates.push_back(native_pose.residue(ii).xyz("CA"));
			test_coordinates.push_back(test_pose.residue(ii).xyz("CA"));
		}
	return(rmsd_between_coordinates(native_coordinates,test_coordinates));
}

Size AnalyzeLoopModeling::get_valid_resid(int resid,core::pose::Pose const pose){
	if ( resid<1 ) {
		TR << "invalid resid encountered n-term" << resid << std::endl;
		resid = 1;
	}
	if ( resid+9-1>(int)pose.total_residue() ) {
		TR << "invalid resid encountered c-term" << resid << std::endl;
		resid = (int)pose.total_residue()-9+1;
	}
	return(Size(resid));
}


Real AnalyzeLoopModeling::generate_lookback_rmsd(core::pose::Pose pose, Size position){
	vector1<Size> resids;
	for(int ii=position-3; ii<=(int)position+2; ii=ii+3){
		Size tmp_resid = get_valid_resid(ii,pose);
			resids.push_back(tmp_resid);
		}
	Real rmsd = SSHashedFragmentStore_->max_rmsd_in_region(pose,resids);
	return(rmsd);
}


void AnalyzeLoopModeling::apply(core::pose::Pose & pose) {
	protocols::loops::Loops pose_loops = get_loops(pose);
	int resAdjustmentStartLow = 0;
	int resAdjustmentStartHigh = 0;
	int resAdjustmentStopLow = 0;
	int resAdjustmentStopHigh = 0;
	int resAdjustmentStartLow_sheet = 0;
	int resAdjustmentStartHigh_sheet = 0;
	int resAdjustmentStopLow_sheet = 0;
	int resAdjustmentStopHigh_sheet = 0;
	Real rmsThreshold=2.5;
	Real max_vdw_change=0;
	bool ideal=false;
	utility::io::ozstream acc_out("loop_accuracy.txt", std::ios_base::app);
	for ( Size ii=pose_loops.num_loop(); ii>=1; --ii ) {
		resAdjustmentStartLow = 0;
		resAdjustmentStartHigh = 0;
		resAdjustmentStopLow = 0;
		resAdjustmentStopHigh = 0;
		resAdjustmentStartLow_sheet = 0;
		resAdjustmentStartHigh_sheet = 0;
		resAdjustmentStopLow_sheet = 0;
		resAdjustmentStopHigh_sheet = 0;
		Size loopStart=pose_loops[ii].start();
		Size loopEnd=pose_loops[ii].stop();
		Size loopLength = loopEnd-loopStart+1;
		std::string tag = jd2::JobDistributor::get_instance()->current_job()->input_tag();
		if(loopLength>=loopLengthRangeLow_ && loopLength<= loopLengthRangeHigh_){
			core::pose::PoseOP kicPoseOP = pose.clone();
			core::pose::PoseOP lookbackPoseOP = pose.clone();
			core::pose::PoseOP lookbackPlusPoseOP = pose.clone();
			NearNativeLoopCloserOP loopCloserKicOP(new NearNativeLoopCloser(resAdjustmentStartLow,resAdjustmentStartHigh,resAdjustmentStopLow,resAdjustmentStopHigh,resAdjustmentStartLow_sheet,resAdjustmentStartHigh_sheet,resAdjustmentStopLow_sheet,resAdjustmentStopHigh_sheet,loopLength,loopLength,pose_loops[ii].start()-1,pose_loops[ii].stop()+1,'A','A',rmsThreshold,max_vdw_change,true,ideal,true,"kic"));
			loopCloserKicOP->apply(*kicPoseOP);
			NearNativeLoopCloserOP lookbackCloserOP(new NearNativeLoopCloser(resAdjustmentStartLow,resAdjustmentStartHigh,resAdjustmentStopLow,resAdjustmentStopHigh,resAdjustmentStartLow_sheet,resAdjustmentStartHigh_sheet,resAdjustmentStopLow_sheet,resAdjustmentStopHigh_sheet,loopLength,loopLength,pose_loops[ii].start()-1,pose_loops[ii].stop()+1,'A','A',rmsThreshold,max_vdw_change,true,ideal,true,"lookback"));
			lookbackCloserOP->apply(*lookbackPoseOP);
			Real kicRmsd = get_loop_rmsd(pose,*kicPoseOP,loopStart-2,loopEnd+2);
			Real lookbackRmsd = get_loop_rmsd(pose,*lookbackPoseOP,loopStart-2,loopEnd+2);
			resAdjustmentStartLow = -3;
			resAdjustmentStartHigh = 3;
			resAdjustmentStopLow = -3;
			resAdjustmentStopHigh = 3;
			NearNativeLoopCloserOP lookbackCloserPlusOP(new NearNativeLoopCloser(resAdjustmentStartLow,resAdjustmentStartHigh,resAdjustmentStopLow,resAdjustmentStopHigh,resAdjustmentStartLow_sheet,resAdjustmentStartHigh_sheet,resAdjustmentStopLow_sheet,resAdjustmentStopHigh_sheet,loopLength,loopLength,pose_loops[ii].start()-1,pose_loops[ii].stop()+1,'A','A',rmsThreshold,max_vdw_change,true,ideal,true,"lookback"));
			lookbackCloserPlusOP->apply(*lookbackPlusPoseOP);
			Real kicVallRmsd = generate_lookback_rmsd(*kicPoseOP,pose_loops[ii].start()-1);
			Real lookbackVallRmsd = generate_lookback_rmsd(*lookbackPoseOP,pose_loops[ii].start()-1);
			Real lookbackPlusVallRmsd = generate_lookback_rmsd(*lookbackPlusPoseOP,pose_loops[ii].start()-1);  //Note: The span is -3 to +3 so some other regions may be considered but this covers the entire region if insertions have occured
			Real poseVallRmsd = generate_lookback_rmsd(pose,pose_loops[ii].start()-1);
			acc_out << tag << "," << pose_loops[ii].start()-1 << "," << loopLength << "," << kicRmsd << "," << lookbackRmsd << "," <<kicVallRmsd << "," << lookbackVallRmsd << "," << lookbackPlusVallRmsd << "," << poseVallRmsd << std::endl;
		}
	}
 }


std::string AnalyzeLoopModeling::get_name() const {
	return "AnalyzeLoopModeling";
}

void
AnalyzeLoopModeling::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & ){
	//start_time_ = time(NULL);
	std::string loopLengthRange( tag->getOption< std::string >( "loopLengthRange", "2,5") );
	utility::vector1< std::string > loopLengthRange_split( utility::string_split( loopLengthRange , ',' ) );
	if ( loopLengthRange_split.size()==2 ) {
		loopLengthRangeLow_ = atoi(loopLengthRange_split[1].c_str());
		loopLengthRangeHigh_ = atoi(loopLengthRange_split[2].c_str());
	}
	SSHashedFragmentStore_ = core::indexed_structure_store::SSHashedFragmentStore::get_instance();
	SSHashedFragmentStore_->set_threshold_distance(.1);
	SSHashedFragmentStore_->init_SS_stub_HashedFragmentStore();
	TR << "database loaded!!" << std::endl;
}

}//pose_length_moves
}//protocols
