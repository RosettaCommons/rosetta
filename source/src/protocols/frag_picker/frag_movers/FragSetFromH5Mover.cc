// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file       protocols/frag_picker/frag_movers/FragSetFromH5Mover.hh
///
/// @brief      Pick fragments from the H5 database.
/// @details    Allows based on ABEGO def,SSdef,PredictedSS,9mers,3mers pushes to stack
///
/// @author     TJ Brunette (tjbrunette@gmail.com)
/// @note

#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>


// Core Headers
#include <protocols/frag_picker/frag_movers/FragSetFromH5Mover.hh>
#include <protocols/frag_picker/frag_movers/FragSetFromH5MoverCreator.hh>

#include <protocols/ss_prediction/SS_predictor.hh>

#include <core/indexed_structure_store/ABEGOHashedFragmentStore.hh>
#include <core/indexed_structure_store/FragmentLookup.hh>
#include <core/indexed_structure_store/FragmentStore.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/io/external/PsiPredInterface.hh>
#include <core/sequence/ABEGOManager.hh>

#include <core/types.hh>

#include <basic/options/keys/remodel.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <string>

static THREAD_LOCAL basic::Tracer TR( "protocols.frag_picker.frag_movers.FragSetFromH5Mover" );

namespace protocols {
namespace frag_picker {
namespace frag_movers {

using utility::vector1;
using namespace core::indexed_structure_store;

std::string FragSetFromH5MoverCreator::keyname() const
{
	return FragSetFromH5MoverCreator::mover_name();
}

std::string FragSetFromH5MoverCreator::mover_name(){
	return "FragSetFromH5Mover";
}

protocols::moves::MoverOP
FragSetFromH5MoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new FragSetFromH5Mover );
}

FragSetFromH5Mover::FragSetFromH5Mover():moves::Mover("FragSetFromH5Mover"){fragSetUnchanged_ = false;}

void FragSetFromH5Mover::convert_ss_to_abego_helper(char tmpChar, vector1<std::string> & abego_strings){
	vector1<std::string> add_abegos;
	if ( tmpChar == 'H' ) {
		add_abegos.push_back("A");
	}
	if ( tmpChar == 'E' ) {
		add_abegos.push_back("B");
	}
	if ( tmpChar == 'L' ) {
		add_abegos.push_back("A");
		add_abegos.push_back("B");
		add_abegos.push_back("E");
		add_abegos.push_back("G");
		add_abegos.push_back("O");
	}
	if ( abego_strings.size()==0 ) {
		for ( Size jj=1; jj<=add_abegos.size(); ++jj ) {
			abego_strings.push_back(add_abegos[jj]);
		}
	} else {
		Size abego_string_init_size = abego_strings.size();
		for ( Size ii=1; ii<=abego_string_init_size; ++ii ) {
			std::string tmpAbegoStr = abego_strings[ii];
			abego_strings[ii]=tmpAbegoStr+add_abegos[1];
			for ( Size jj=2; jj<=add_abegos.size(); ++jj ) {
				abego_strings.push_back(tmpAbegoStr+add_abegos[jj]);
			}
		}
	}
}

vector1 <std::string> FragSetFromH5Mover::convert_ss_to_abegos(std::string ss){
	vector1<std::string> abego_strings;
	std::cout <<"XXHERE:" << ss << std::endl;
	if ( ss.size()>9 ) {
		utility_exit_with_message("convert_ss_to_abegos only operates on short strings < 9 residues because of possible computational explosion.");
	} else {
		//9 loop residues would be equal to 1953125 - 2(all helix and all abego) so lots of structures!
		for ( Size ii=0; ii<ss.size(); ii++ ) {
			char tmpChar = ss.at(ii);
			convert_ss_to_abego_helper(tmpChar,abego_strings);
		}
	}
	//delete all A and B strings if loop present. The all A and B are huge so this will help a lot.
	if ( ss.find("L")!=std::string::npos ) {
		vector1<std::string>::iterator location = std::find(abego_strings.begin(),abego_strings.end(),"AAAAAAAAA");
		if ( location!=abego_strings.end() ) {
			abego_strings.erase(location);
		}
		location = std::find(abego_strings.begin(),abego_strings.end(),"BBBBBBBBB");
		if ( location!=abego_strings.end() ) {
			abego_strings.erase(location);
		}
	}
	return(abego_strings);
}

struct FragLocation{
	std::string abegoID;
	core::indexed_structure_store::FragmentLookupResult lookupResult;
	FragLocation(core::indexed_structure_store::FragmentLookupResult lookupResult_,std::string abegoID_){
		abegoID = abegoID_;
		lookupResult = lookupResult_;
	}
};

vector1< vector1< core::Real > >  FragSetFromH5Mover::get_ss_prediction(const core::pose::Pose & pose){
	using namespace core::io::external;
	utility::vector1< utility::vector1< core::Real > >  prediction_by_pct;
	//The psipred tool does not output %helical etc.  So I'm going to use the svm
	if ( use_svm_ ) {
		std::string sequence;
		for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
			if ( pose.residue( i ).is_protein() ) sequence += pose.residue( i ).name1();
		}
		prediction_by_pct = ss_predictor_->predict_ss( sequence );
	}
	/* else {
		runtime_assert( psipred_interface_ != 0 );
		std::string currentSS = pose.secstruct();
		PsiPredResult const psipred_result = psipred_interface_->run_psipred( pose, currentSS );
		for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
			prediction_by_pct.push_back(psipred_result.psipred_prob[i]);
		}
	}
	*/
	return(prediction_by_pct);
}

bool FragSetFromH5Mover::fragSet_needs_update(){
	if ( datamap_.has("FragSets",fragSetName_) && (!use_pose_) ) {
		return(false);
	} else {
		return(true);
	}
}
core::pose::Pose FragSetFromH5Mover::get_selected_pose(const core::pose::Pose pose){
	if(chain_.size()!=0){
		if(!has_chain(chain_,pose)){
			utility_exit_with_message("chain does not exist");
		}
		else{
			Size chain_id =  get_chain_id_from_chain(chain_,pose);
			return(*pose.split_by_chain(chain_id));
		}
	}
	return(pose);
}
vector1 <std::string> FragSetFromH5Mover::get_ss_strings_for_residue_from_ssPred(vector1<vector1<Real> > ss_prediction,Size position){
	vector1<std::string> allowed_ss;
	if(ss_prediction[position][1]>=ssPred_cutoff_)
		allowed_ss.push_back("H");
	if(ss_prediction[position][2]>=ssPred_cutoff_)
		allowed_ss.push_back("L");
	if(ss_prediction[position][3]>=ssPred_cutoff_)
		allowed_ss.push_back("E");
	return(allowed_ss);
}


vector1 <std::string> FragSetFromH5Mover::get_abego_strings_for_residue_from_ssPred(vector1<vector1<Real> > ss_prediction,Size res,Size fragLength){
	vector1 <std::string> abegos_to_check;
	for(Size ii=res; ii<=res+fragLength; ++ii){
		vector1<std::string> allowed_ss = get_ss_strings_for_residue_from_ssPred(ss_prediction,ii);
		if(abegos_to_check.size() == 0){
			for(Size jj=1; jj<=allowed_ss.size(); ++jj)
				abegos_to_check.push_back(allowed_ss[jj]);
		}
		else{
			Size abegos_to_check_init_size = abegos_to_check.size();
			for(Size kk=1; kk<=abegos_to_check_init_size; ++kk){
				std::string tmpString = abegos_to_check[kk];
				abegos_to_check[kk]=abegos_to_check[kk]+allowed_ss[1];
				for(Size jj=2; jj<=allowed_ss.size(); ++jj){
					abegos_to_check.push_back(tmpString+allowed_ss[jj]);
				}
			}
		}
	}
	return(abegos_to_check);
}


vector1 <std::string> FragSetFromH5Mover::get_abego_strings_for_residue(const core::pose::Pose pose,Size res,Size fragLength){
	vector1 <std::string> abegos_to_check;
	if(use_pose_){
		//use ABEGO manager
		core::sequence::ABEGOManager AM;
		utility::vector1< std::string > abegoSeq = AM.get_symbols( pose,1 );//1 stands for class of ABEGO strings
		std::string tmp;
		for(Size ii=res; ii<=res+fragLength; ++ii)
			tmp += abegoSeq[ii];
		abegos_to_check.push_back(tmp);
	}
	else{
		if(ss_.size()>0){
			abegos_to_check = convert_ss_to_abegos(ss_.substr(res-1,fragLength));
		}
		if(abego_.size()>0){
			abegos_to_check.push_back(abego_.substr(res-1,fragLength));
		}
	}
	return(abegos_to_check);
}


void FragSetFromH5Mover::apply(core::pose::Pose & pose){
	//step0 decide if fragset needs updating
	std::cout <<"*******here" << fragSet_needs_update() << std::endl;
	if(fragSet_needs_update()){
		//step1 get pose & abegos
		vector1 <std::string>  abegos_to_check;
		Size numbRes = 1;
		core::pose::Pose chosenPose = pose;
		if(use_pose_){
			chosenPose = get_selected_pose(pose);
			numbRes = chosenPose.total_residue();
		}
		else{
			if(ss_.size()>0)
				numbRes = ss_.size();
			if(abego_.size()>0)
				numbRes = abego_.size();
		}
		std::cout << "ss_.size()" <<  ss_.size() << "abego_" << abego_.size() << "numbRes" << numbRes << std::endl;
		Size fragLength = 9;
		vector1<vector1< Real> > ss_prediction;
		if(use_pose_ && use_ssPred_)
			ss_prediction = get_ss_prediction(chosenPose);
		for(Size ii=1; ii<=numbRes-fragLength+1; ++ii){
			vector1<string> abegos_to_check;
			if(use_pose_ && use_ssPred_)
				abegos_to_check = get_abego_strings_for_residue_from_ssPred(ss_prediction,ii,fragLength);
			else
				abegos_to_check = get_abego_strings_for_residue(chosenPose,ii,fragLength);
			for(Size jj=1; jj<=abegos_to_check.size(); ++jj)
				std::cout << "pos:" << ii << " abego:" <<  abegos_to_check[jj] << std::endl;
		}

	}
}


std::string FragSetFromH5Mover::get_name() const {
	return "FragSetFromH5Mover";
}

void
FragSetFromH5Mover::parse_my_tag(

	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & ){
	using namespace core::indexed_structure_store;
	/*
	ABEGOHashedFragmentStore_ = ABEGOHashedFragmentStoreOP(new ABEGOHashedFragmentStore(rmsThreshold_));
	*/
	if ( tag->hasOption("abego") ) {
		abego_ = tag->getOption<std::string>("abego");
	} else {
		abego_="";
	}
	if ( tag->hasOption("ss") ) {
		ss_ = tag->getOption<std::string>("ss");
	} else {
		ss_ = "";
	}
	use_ssPred_ = false; //True only if usePose && ssPred are set
	if ( tag->hasOption("use_pose") ) {
		use_pose_ = true;
		if (tag->hasOption("use_ssPred")){
			ssPred_cutoff_ = tag->getOption<Real>( "ssPred_cutoff", 0.3 );
			use_ssPred_=true;
			use_svm_ = true; //psipred interface not yet ready
			/*
			use_svm_ = tag->getOption< bool >( "use_svm_ssPred", use_svm_);
			use_psipred_ = tag->getOption< bool >( "use_psipred_ssPred", use_psipred_);
			std::string psi_cmd_ = tag->getOption< std::string >( "psipred_cmd", "" );
			psipred_interface_ = PsiPredInterfaceOP( new PsiPredInterface( psi_cmd_ ) );
			*/
		}
		else{
			use_svm_ = false;
		}
	} else {
		use_pose_ = false;
		use_ssPred_ = false;
	}
	if ( use_svm_ ) {
		ss_predictor_ = protocols::ss_prediction::SS_predictorOP( new protocols::ss_prediction::SS_predictor( "HLE" ) );
	}
	chain_ = tag->getOption<std::string>("chain","");
	nFrag_ = tag->getOption<Size>( "nFrag", 5 );
	fragSetName_ = tag->getOption<std::string>("name");
}

}//frag_movers
}//frag_picker
}//protocols
