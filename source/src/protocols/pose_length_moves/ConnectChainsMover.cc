// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/moves/ConnectChainsMover.fwd.hh
/// @details connects chains and returns the connection with lowest rmsd. Relies on NearestNativeCloser
/// @author TJ Brunette tjbrunette@gmail.com
///
// Unit headers
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/pose/symmetry/util.hh>

#include <protocols/pose_length_moves/ConnectChainsMover.hh>
#include <protocols/pose_length_moves/ConnectChainsMoverCreator.hh>
#include <protocols/moves/Mover.hh>

#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <protocols/pose_length_moves/NearNativeLoopCloser.hh>

// Core Headers
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/ReferenceCount.hh>



#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <ctime>
#include <boost/algorithm/string.hpp>
#include <boost/range/algorithm/count.hpp>


static THREAD_LOCAL basic::Tracer TR( "protocols.pose_length_moves.ConnectChainsMover" );

namespace protocols {
namespace pose_length_moves {
using namespace core;
using namespace std;
using utility::vector1;



ConnectChainsMover::ConnectChainsMover():moves::Mover("ConnectChainsMover"){}

std::string ConnectChainsMoverCreator::keyname() const
{
	return ConnectChainsMoverCreator::mover_name();
}

std::string ConnectChainsMoverCreator::mover_name(){
	return "ConnectChainsMover";
}


protocols::moves::MoverOP
ConnectChainsMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new ConnectChainsMover );
}


void ConnectChainsMover::parse_input(vector1<std::string> & individual_chains,vector1< vector1 <std::string> > & chains_in_poses){
	set<std::string> tmp_chains_set;
	boost::replace_all(output_chains_,"],[","|");
	boost::replace_all(output_chains_,"[","");
	boost::replace_all(output_chains_,"]","");
	utility::vector1< std::string > output_strings( utility::string_split( output_chains_ , '|' ) );
	for ( Size ii=1; ii<=output_strings.size(); ii++ ) {
		utility::vector1< std::string > chains_in_pose(utility::string_split( output_strings[ii] , ',' ) );
		chains_in_poses.push_back(chains_in_pose);
		for ( Size jj=1; jj<=chains_in_pose.size(); ++jj ) {
			tmp_chains_set.insert(chains_in_pose[jj]);
			//get individual chains
			utility::vector1< std::string > single_residue_chains(utility::string_split( chains_in_pose[jj] , '+' ) );
			for ( Size kk=1; kk<=single_residue_chains.size(); ++kk ) {
				tmp_chains_set.insert(single_residue_chains[kk]);
			}
		}
	}
	for (const auto & iter : tmp_chains_set) {
		individual_chains.push_back(iter);
		TR << "chains" << iter << std::endl;
	}

	std::sort(individual_chains.begin(),individual_chains.end(),chain_lt());
}


map<std::string, Chain> ConnectChainsMover::generate_connected_chains(core::pose::Pose const pose,vector1<std::string> individual_chains){
	map<std::string, Chain> connected_chains;
	for ( Size ii=1; ii<=individual_chains.size(); ++ii ) {
		Size chain_length = boost::count(individual_chains[ii],'+');
		if ( chain_length ==0 ) {
			Size chain_id =  get_chain_id_from_chain(individual_chains[ii],pose);
			core::pose::PoseOP chain = pose.split_by_chain(chain_id);
			struct Chain chain_tmp(chain,0);
			connected_chains.insert(std::pair<std::string, Chain>(individual_chains[ii],chain_tmp));
		} else {
			assemble_missing_chain(connected_chains,"",individual_chains[ii]);
		}
	}
	return(connected_chains);
}

void ConnectChainsMover::assemble_missing_chain(map<std::string, Chain> & connected_chains,std::string chain_assembled,std::string chain_remainder){
	Size chain_remainder_length = boost::count(chain_remainder,'+');
	Size chain_assembled_length = boost::count(chain_assembled,'+');
	if ( chain_remainder.size()!=0 ) {
		chain_remainder_length++;
	}
	if ( chain_assembled.size()!=0 ) {
		chain_assembled_length++;
	}
	vector1< std::string > chain_remainder_split(utility::string_split(chain_remainder , '+' ) );
	if ( chain_remainder_length == 0 ) {
		return;
	} else {
		std::string to_assemble_string = "";
		if ( chain_assembled_length==0 ) {
			to_assemble_string = chain_remainder_split[1];
		} else {
			to_assemble_string = chain_assembled + "+" + chain_remainder_split[1];
		}
		if ( connected_chains.find(to_assemble_string) == connected_chains.end() ) {
			TR << "connecting" << chain_assembled << " to " << chain_remainder_split[1] << std::endl;
			//not found. connect chains & put in connected chains
			core::pose::PoseOP chainA = connected_chains.at(chain_assembled).poseOP;
			Real chainA_rmsd = connected_chains.at(chain_assembled).rmsd;
			core::pose::PoseOP chainA_plus = connected_chains.at(chain_assembled).poseOP->clone();
			core::pose::PoseOP chainB = connected_chains.at(chain_remainder_split[1]).poseOP;
			Real chainB_rmsd = connected_chains.at(chain_assembled).rmsd;
			NearNativeLoopCloserOP loopCloserOP(new NearNativeLoopCloser(resAdjustmentStartLow_,resAdjustmentStartHigh_,resAdjustmentStopLow_,resAdjustmentStopHigh_,resAdjustmentStartLow_sheet_,resAdjustmentStartHigh_sheet_,resAdjustmentStopLow_sheet_,resAdjustmentStopHigh_sheet_,loopLengthRangeLow_,loopLengthRangeHigh_,1,1,'A','B',rmsThreshold_,0,true,false,true));
			append_pose_to_pose(*chainA_plus,*chainB,true);
			renumber_pdbinfo_based_on_conf_chains(*chainA_plus,true,false,false,false);
			utility::vector1< char > pdb_chains;
			for ( Size ii=1; ii<=chainA->total_residue(); ++ii ) {
				pdb_chains.push_back('A');
			}
			for ( Size ii=1; ii<=chainB->total_residue(); ++ii ) {
				pdb_chains.push_back('B');
			}
			chainA_plus->pdb_info()->set_chains(pdb_chains);
			if ( chainA_rmsd > rmsThreshold_ || chainB_rmsd > rmsThreshold_ ) { //no closure needed if rmsd is > threshold. But fill in the DB so this is not re-tried
				Real tmp_rmsd = chainA_rmsd;
				if ( tmp_rmsd < chainB_rmsd ) {
					tmp_rmsd = chainB_rmsd;
				}
				struct Chain chain_tmp(chainA_plus,tmp_rmsd);
				connected_chains.insert(std::pair<std::string, Chain>(to_assemble_string,chain_tmp));
			} else {
				Real return_rmsd = loopCloserOP->close_loop(*chainA_plus);
				if ( return_rmsd < chainA_rmsd ) {
					return_rmsd = chainA_rmsd;
				}
				if ( return_rmsd < chainB_rmsd ) {
					return_rmsd = chainB_rmsd;
				}
				struct Chain chain_tmp(chainA_plus,return_rmsd);
				connected_chains.insert(std::pair<std::string, Chain>(to_assemble_string,chain_tmp));
			}
		}
		std::string next_remainder = "";
		if ( chain_remainder_length>1 ) {
			for ( Size ii=2; ii<chain_remainder_length; ii++ ) {
				next_remainder+=chain_remainder_split[ii] + "+";
			}
			next_remainder+=chain_remainder_split[chain_remainder_length];
		}
		assemble_missing_chain(connected_chains,to_assemble_string,next_remainder);
	}
}

void ConnectChainsMover::generate_best_final_pose(core::pose::Pose & pose,vector1< vector1 <std::string> > chains_in_poses,map<std::string, Chain> connected_chains){
	Real low_rmsd = 9999999;
	Size low_rmsd_pose_index = 0;
	std::string low_rmsd_chain_string = "";
	for ( Size ii=1; ii<=chains_in_poses.size(); ++ii ) {
		Real tmp_pose_rmsd = 0;
		std::string tmp_chain_string = "";
		for ( Size jj=1; jj<=chains_in_poses[ii].size(); ++jj ) {
			Real tmp_chain_rmsd = connected_chains.at(chains_in_poses[ii][jj]).rmsd;
			tmp_chain_string += chains_in_poses[ii][jj] + ",";
			if ( tmp_pose_rmsd<tmp_chain_rmsd ) {
				tmp_pose_rmsd = tmp_chain_rmsd;
			}
		}
		if ( tmp_pose_rmsd < low_rmsd ) {
			low_rmsd = tmp_pose_rmsd;
			low_rmsd_pose_index = ii;
			low_rmsd_chain_string = tmp_chain_string;
		}
	}
	if ( low_rmsd > rmsThreshold_ ) {
		utility_exit_with_message("No loop closures found below threshold, exiting");
	} else {
		TR << "low rmsd to nine residue region in VALL: " << low_rmsd << " description of output chain:" << low_rmsd_chain_string << std::endl;
		core::pose::PoseOP return_pose = connected_chains.at(chains_in_poses[low_rmsd_pose_index][1]).poseOP->clone();
		for ( Size jj=2; jj<=chains_in_poses[low_rmsd_pose_index].size(); ++jj ) {
			core::pose::PoseOP append_pose = connected_chains.at(chains_in_poses[low_rmsd_pose_index][jj]).poseOP;
			append_pose_to_pose(*return_pose,*append_pose,true);
		}
		pose = *return_pose;
	}
}

void ConnectChainsMover::apply(core::pose::Pose & pose) {
	vector1<std::string> individual_chains;
	vector1< vector1 <std::string> > chains_in_poses;
	parse_input(individual_chains, chains_in_poses);
	map<std::string, Chain> connected_chains = generate_connected_chains(pose,individual_chains);
	generate_best_final_pose(pose,chains_in_poses,connected_chains);
}


std::string ConnectChainsMover::get_name() const {
	return "ConnectChainsMover";
}

void
ConnectChainsMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & ){
	//start_time_ = time(NULL);
	std::string loopLengthRange( tag->getOption< std::string >( "loopLengthRange", "1,4") );
	rmsThreshold_ = tag->getOption< core::Real >( "RMSthreshold", 0.4 );
	std::string resAdjustmentRange1( tag->getOption< std::string >( "resAdjustmentRangeSide1", "-3,3") );
	std::string resAdjustmentRange2( tag->getOption< std::string >( "resAdjustmentRangeSide2","-3,3") );
	std::string resAdjustmentRange1_sheet( tag->getOption< std::string >( "resAdjustmentRangeSide1_sheet", "-1,1") );
	std::string resAdjustmentRange2_sheet( tag->getOption< std::string >( "resAdjustmentRangeSide2_sheet","-1,1") );
	output_chains_ = tag->getOption<std::string>("chain_connections");
	utility::vector1< std::string > resAdjustmentRange1_split( utility::string_split( resAdjustmentRange1 , ',' ) );
	utility::vector1< std::string > resAdjustmentRange2_split( utility::string_split( resAdjustmentRange2 , ',' ) );
	utility::vector1< std::string > resAdjustmentRange1_sheet_split( utility::string_split( resAdjustmentRange1_sheet , ',' ) );
	utility::vector1< std::string > resAdjustmentRange2_sheet_split( utility::string_split( resAdjustmentRange2_sheet , ',' ) );
	utility::vector1< std::string > loopLengthRange_split( utility::string_split( loopLengthRange , ',' ) );
	if ( resAdjustmentRange1_split.size()==2 ) {
		resAdjustmentStartLow_ = atoi(resAdjustmentRange1_split[1].c_str());
		resAdjustmentStartHigh_ = atoi(resAdjustmentRange1_split[2].c_str());
	}
	if ( resAdjustmentRange1_split.size()==1 ) {
		resAdjustmentStartLow_= atoi(resAdjustmentRange1_split[1].c_str());
		resAdjustmentStartHigh_= atoi(resAdjustmentRange1_split[1].c_str());
	}
	if ( resAdjustmentRange2_split.size()==2 ) {
		resAdjustmentStopLow_ = atoi(resAdjustmentRange2_split[1].c_str());
		resAdjustmentStopHigh_ = atoi(resAdjustmentRange2_split[2].c_str());
	}
	if ( resAdjustmentRange2_split.size()==1 ) {
		resAdjustmentStopLow_ = atoi(resAdjustmentRange2_split[1].c_str());
		resAdjustmentStartHigh_ = atoi(resAdjustmentRange2_split[1].c_str());
	}
	if ( resAdjustmentRange1_sheet_split.size()==2 ) {
		resAdjustmentStartLow_sheet_ = atoi(resAdjustmentRange1_sheet_split[1].c_str());
		resAdjustmentStartHigh_sheet_ = atoi(resAdjustmentRange1_sheet_split[2].c_str());
	}
	if ( resAdjustmentRange1_sheet_split.size()==1 ) {
		resAdjustmentStartLow_sheet_= atoi(resAdjustmentRange1_sheet_split[1].c_str());
		resAdjustmentStartHigh_sheet_= atoi(resAdjustmentRange1_sheet_split[1].c_str());
	}
	if ( resAdjustmentRange2_sheet_split.size()==2 ) {
		resAdjustmentStopLow_sheet_ = atoi(resAdjustmentRange2_sheet_split[1].c_str());
		resAdjustmentStopHigh_sheet_ = atoi(resAdjustmentRange2_sheet_split[2].c_str());
	}
	if ( resAdjustmentRange2_sheet_split.size()==1 ) {
		resAdjustmentStopLow_sheet_ = atoi(resAdjustmentRange2_sheet_split[1].c_str());
		resAdjustmentStartHigh_sheet_ = atoi(resAdjustmentRange2_sheet_split[1].c_str());
	}
	if ( loopLengthRange_split.size()==2 ) {
		loopLengthRangeLow_ = atoi(loopLengthRange_split[1].c_str());
		loopLengthRangeHigh_ = atoi(loopLengthRange_split[2].c_str());
	}
	if ( loopLengthRange_split.size()==1 ) {
		loopLengthRangeLow_ = atoi(loopLengthRange_split[1].c_str());
		loopLengthRangeHigh_ = atoi(loopLengthRange_split[1].c_str());
	}
}

}//pose_length_moves
}//protocols
