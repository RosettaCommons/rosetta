/// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ss_predictor/HLE_predicot
/// @brief  A simple svm based secondary structure predictor
/// @author TJ Brunette

// Utility Headers
#include <basic/database/open.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>

#include <core/types.hh>
#include <utility/vector1.hh>
// C++ headers
#include <iostream>
#include <cstring>
#include <cmath>
#include <string>
#include <utility/libsvm/Svm_rosetta.hh>
#include <utility/io/izstream.hh>

#include <protocols/ss_prediction/SS_predictor.hh>

namespace protocols{
namespace ss_prediction{

using std::string;
using utility::vector1;
using core::Real;
using core::Size;
using core::SSize;
using namespace utility::libsvm;

basic::Tracer tr( "ss_prediction.ss_predictor" );

/////////////////////////////////////////////////////////////////////////////////
//@brief: This class does single sequence prediction for either HLE or ABEGO
/////////////////////////////////////////////////////////////////////////////////
SS_predictor::SS_predictor(string type){//type = HLE or ABEGO
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
	ss_type = type;
  string ss_predictor_db_dir = basic::database::full_name( "external/svm_models/ss_prediction" );

	string rd1_model_fl = ss_predictor_db_dir + "/" + type + "_rd1_model";
	string rd2_model_fl = ss_predictor_db_dir + "/" + type + "_rd2_model";
	tr << "loading:" << rd1_model_fl<< std::endl;
	tr << "loading:" << rd2_model_fl<< std::endl;

	utility::io::izstream is( rd1_model_fl );
	if ( !is.good() ) {
		utility_exit_with_message("Error: Make sure you download the svm_models from https://svn.rosettacommons.org/source/trunk/svm_models and put it into rosetta_database/external/svm_models");
	}
	load_models(rd1_model_fl,rd2_model_fl);
}

SS_predictor::~SS_predictor(){
}
/////////////////////////////////////////////////////////////////////////////////
//@brief: Loads the models
///////////////////////////////////////////////////////////////////////////////// =
void SS_predictor::load_models(string rd1_model_fl, string rd2_model_fl){
	const char* rd1_model_fl_char= rd1_model_fl.c_str();
	const char* rd2_model_fl_char= rd2_model_fl.c_str();
	rd1_model = Svm_rosettaOP(new Svm_rosetta(rd1_model_fl_char));
	rd2_model = Svm_rosettaOP(new Svm_rosetta(rd2_model_fl_char));
}
/////////////////////////////////////////////////////////////////////////////////
//@brief: Gets 15 residue window from the fasta. ALso replaces non AA with X. This
//happens on termini
/////////////////////////////////////////////////////////////////////////////////
string SS_predictor::get_window_aa(string fasta,SSize position){
	SSize half_window_size = (SSize)floor((Real)WINDOW_SIZE/2.0);
	string window_aa = "";
	for(SSize ii=position-half_window_size; ii<=position+half_window_size; ++ii){
		if(ii<0 || ii>=SSize(fasta.size()))
			window_aa += "X";
		else
			window_aa += fasta.at(ii);
	}
	return(window_aa);
}
/////////////////////////////////////////////////////////////////////////////////
//@brief: Prediction for a single location in rd1
/////////////////////////////////////////////////////////////////////////////////
vector1 <Real> SS_predictor::predict_pos_rd1(string window_aa){
	string aa_order =  "ARNDCQEGHILKMFPSTWYVX";
	vector1< Svm_node_rosettaOP > features;
	for(Size ii=0; ii<window_aa.size(); ++ii){
		Size aa_pos_index = (aa_order.find(window_aa.at(ii))+1);
		Size final_index = ii*21+aa_pos_index;
		Svm_node_rosettaOP tmpNode = Svm_node_rosettaOP(new Svm_node_rosetta(final_index,1));
		features.push_back(tmpNode);
	}
	vector1< Real> probs_to_return = rd1_model->predict_probability(features);
	return(probs_to_return);
}
/////////////////////////////////////////////////////////////////////////////////
//@brief: Prediction for a single location in rd2
/////////////////////////////////////////////////////////////////////////////////
vector1 <Real> SS_predictor::predict_pos_rd2( vector1<vector1<Real> > rd1_preds,SSize position,string fasta){
	string ABEGO_RD1_OUT_ORDER="EBAGO"; // order the probabilities are stored in rd1
	string HLE_RD1_OUT_ORDER="LEH"; //rd2 expects the probabilities to be ABEGO,HLE
	string ABEGO_RD2_IN_ORDER ="ABEGO";
	string HLE_RD2_IN_ORDER = "HLE";
	string ABEGO_RD2_OUT_ORDER ="EBAGO";
	string HLE_RD2_OUT_ORDER = "LEH";
	SSize half_window_size = (SSize)floor((Real)WINDOW_SIZE/2.0);
	Size lpCt = 0;
	vector1< Svm_node_rosettaOP > features;
	for(SSize ii=position-half_window_size; ii<=position+half_window_size; ++ii){
		if(ii<0 || ii>=SSize(fasta.size())){
			Size tmp_index = lpCt*(rd1_preds[1].size()+1)+(rd1_preds[1].size()+1);//For HLE position 4
			Svm_node_rosettaOP tmpNode = Svm_node_rosettaOP(new Svm_node_rosetta(tmp_index,1));
			features.push_back(tmpNode);
		}
		else
			for(Size jj=0; jj<rd1_preds[1].size(); ++jj){
				Size rd1_pred_typeCorrection = 99999;
				if(ss_type == "HLE")
					rd1_pred_typeCorrection = HLE_RD1_OUT_ORDER.find(HLE_RD2_IN_ORDER.at(jj))+1;
				else
					rd1_pred_typeCorrection = ABEGO_RD1_OUT_ORDER.find(ABEGO_RD2_IN_ORDER.at(jj))+1;
				Size tmp_index = lpCt*(rd1_preds[1].size()+1)+jj+1;
				Real tmp_value = rd1_preds[ii+1][rd1_pred_typeCorrection];
				Svm_node_rosettaOP tmpNode= Svm_node_rosettaOP(new Svm_node_rosetta(tmp_index,tmp_value));
				features.push_back(tmpNode);
			}
		lpCt++;
	}
	vector1 <Real> prob_estimates = rd2_model->predict_probability(features);
	vector1 <Real> probs_to_return;
	for(Size ii=0;ii<ss_type.size(); ++ii){
		Size rd1_pred_typeCorrection = 99999;
		if(ss_type == "HLE") //The order that the probabilities come out is not ABEGO HLE, this corrects.
			rd1_pred_typeCorrection = HLE_RD2_OUT_ORDER.find(HLE_RD2_IN_ORDER.at(ii));
		else
			rd1_pred_typeCorrection = ABEGO_RD2_OUT_ORDER.find(ABEGO_RD2_IN_ORDER.at(ii));
		probs_to_return.push_back(prob_estimates[rd1_pred_typeCorrection+1]);
	}
	return(probs_to_return);
}
/////////////////////////////////////////////////////////////////////////////////
//@brief: Rd1 prediction
/////////////////////////////////////////////////////////////////////////////////
vector1< vector1 <Real> > SS_predictor::predict_rd1(string fasta){
	vector1< vector1 <Real> > probs;
	for(Size ii=0; ii<fasta.size(); ++ii){
		string window_aa = get_window_aa(fasta,ii);
		probs.push_back(predict_pos_rd1(window_aa));
	}
	return probs;
}
/////////////////////////////////////////////////////////////////////////////////
//@brief: Rd2 prediction
/////////////////////////////////////////////////////////////////////////////////
vector1< vector1 <Real> > SS_predictor::predict_rd2(vector1< vector1< Real > > rd1_preds, string fasta){
	vector1< vector1 <Real> > probs;
	Size correction = 0;
	if(ss_type == "ABEGO")//Note the last position in ABEGO is ALWAYS a O
		correction = 1;
	for(Size ii=0; ii<fasta.size()-correction; ++ii){
		probs.push_back(predict_pos_rd2(rd1_preds,ii,fasta));
	}
	if(correction == 1){
		vector1 <Real> correction_probs;
		for(int ii=0;ii<4;++ii)
			correction_probs.push_back(0.0);
		correction_probs.push_back(1.0);
		probs.push_back(correction_probs);
	}
	return probs;
}
/////////////////////////////////////////////////////////////////////////////////
//@brief: Predicts the SS in a two round fashion. Probabilities are output ABEGO order or HLE order depending on type
/////////////////////////////////////////////////////////////////////////////////
vector1<vector1 <Real> > SS_predictor::predict_ss(string fasta){
	vector1<vector1 <Real> > rd1_preds = predict_rd1(fasta);
	vector1<vector1 <Real> > rd2_preds = predict_rd2(rd1_preds,fasta);
	return(rd2_preds);
}

/// @brief helper function to get SS char at a position from a vector of reals
char get_label( utility::vector1< core::Real > const & ss_pred_pos )
{
	char label = 'X';
	if((ss_pred_pos[1] >= ss_pred_pos[2]) && (ss_pred_pos[1] >= ss_pred_pos[3]))
		label = 'H';
	else
		if((ss_pred_pos[2] >= ss_pred_pos[1]) && (ss_pred_pos[2] >= ss_pred_pos[3]))
			label = 'L';
		else
			if((ss_pred_pos[3] >= ss_pred_pos[1]) && (ss_pred_pos[3] >= ss_pred_pos[2]))
				label = 'E';
	return(label);
}

/// @brief helper function to get SS char at a position from a vector of reals
core::Real get_prob( char wanted_ss, utility::vector1< core::Real > const & ss_pred_pos )
{
	// order is "H" "L" "E"
	if ( wanted_ss == 'H' ) {
		return ss_pred_pos[1];
	} else if ( wanted_ss == 'L' ) {
		return ss_pred_pos[2];
	} else if ( wanted_ss == 'E' ) {
		return ss_pred_pos[3];
	} else {
		utility_exit_with_message( "Error: secondary structure is something other than 'H', 'L', or 'E': " +
				std::string(1, wanted_ss) );
		return 99;
	}
}

} //ss_prediction
} //protocols
