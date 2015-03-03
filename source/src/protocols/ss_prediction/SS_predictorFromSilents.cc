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

#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/util.hh>
#include <core/fragment/SecondaryStructure.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
// C++ headers
#include <iostream>
#include <cstring>
#include <cmath>
#include <string>
#include <utility/libsvm/Svm_rosetta.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <iostream>

#include <protocols/ss_prediction/SS_predictor.hh>

namespace protocols{
namespace ss_prediction{

using std::string;
using utility::vector1;
using core::Real;
using core::Size;
using core::SSize;
using namespace utility::libsvm;

static thread_local basic::Tracer tr( "ss_prediction.ss_predictor" );

/////////////////////////////////////////////////////////////////////////////////
//@brief: This class both extracts the data for training if needed and does the prediction after that is done.
/////////////////////////////////////////////////////////////////////////////////
SS_predictorFromSilents::SS_predictorFromSilents(string type){//type = "train" or "run"
		if(type == "run"){
				std::cout << "not implemented" << std::endl;
		}
}

/////////////////////////////////////////////////////////////////////////////////
//@brief: gets the psipred data into C,H,E data 
///////////////////////////////////////////////////////////////////////////////////
core::fragment::SecondaryStructureOP ss_predictorFromSilents::get_psipred_ss2(){
		core::fragment::SecondaryStructureOP pred_psipred_ = core::fragment::SecondaryStructureOP( new SecondaryStructure() );
		pred_psipred_->read_psipred_ss2(option[in::file::psipred_ss2]());
		return(pred_psipred_);
}
/////////////////////////////////////////////////////////////////////////////////
//@brief: gets the silent files to secondarystructureOP object 
//////////////////////////////////////////////////////////////////////////////////
core::fragment::SecondaryStructureOP ss_predictorFromSilents::get_silents_ss2(){
		using namespace core::chemical;
		using namespace core::import_pose::pose_stream;
		core::fragment::SecondaryStructureOP pred_psipred_ = core::fragment::SecondaryStructureOP( new SecondaryStructure() );
		MetaPoseInputStream input = streams_from_cmd_line();
    ResidueTypeSetCOP rsd_set( rsd_set_from_cmd_line() );
		ConstantLengthFragSetOP fragset_ = new ConstantLengthFragSet( 1 );
		while(input.has_another_pose()){
				core::pose::PoseOP input_poseOP;
				input_poseOP = core::pose::PoseOP( new core::pose::Pose() );
				input.fill_pose(*input_poseOP,*rsd_set);
				steal_constant_length_frag_set_from_pose (*input_poseOP,fragSet_);
		}
		core::fragment::SecondaryStructureOP pred_structs_ = core::fragment::SecondaryStructureOP( new SecondaryStructure(*fragments_) );
		return(pred_structs_);
}
/////////////////////////////////////////////////////////////////////////////////
//@brief: See if there is only 1 sheet
////////////////////////////////////////////////////////////////////////////////////
bool ss_predictorFromSilents::only_1_sheet_chk(core::fragment::SecondaryStructureOP pred_psipred_){
		Size numb_potential_sheets = 0;
		Real sheet_threshold = 0.3;
		bool in_loop=false;
		for(Size ii=0; ii<=pred_psipred_->total_residue(); ++ii){
				if((pred_psipred_->sheet_fraction(ii)>sheet_threshold )&&(!in_sheet)){ 
						in_sheet=true;
						numb_potential_sheets++;
				}
				if((pred_psipred_->sheet_fraction(ii)<sheet_threshold )&&(in_sheet)){
						in_loop=false;
				}
		}
		if(in_loop)
				numb_potential_sheets++;
		if(numb_potential_sheets == 1)
				return(true);
		else
				return(false);
}

/////////////////////////////////////////////////////////////////////////////////
//@brief: This class both extracts the data for training if needed and does the prediction after that is done.
/////////////////////////////////////////////////////////////////////////////////
void SS_predictorFromSilents::output_training_data(){
		int range = 3;
		//step1 get data from psipred
		core::fragment::SecondaryStructureOP pred_psipred_ = get_psipred_ss2();
		//step2 get the data from the silent files. 
		core::fragment::SecondaryStructureOP pred_silents_ = get_silents_ss2();			
		//step3 determine if only 1 loop region predicted above 30% by psipred
		bool oneSheet = only_1_sheet_chk(pred_psipred_);
		//step4 get native dssp
		Pose native_pose;
			pose_from_pdb(
			native_pose,
			*rsd_set,
			option[ in::file::native ]()
			);
		core::scoring::dssp::Dssp dssp_obj( native_pose );
		dssp_obj.insert_ss_into_pose( native_pose );
		//output all data into 1 and 2 data lines
		std::string tag = core::pose::tag_from_pose(*input_poseOP);
		std::string out1_fl = (tag + ".1.data");
		std::string out15_fl = (tag+ ".15.data");
		utility::io::ozstream out1(out1_fl);
		utility::io::ozstream out15(out15_fl);
		//for consistency with my previous SVM I have chosen the following "L" type 2, "H" type 1, "E" type 3
		//out1: 
		//type 1:oneSheet 2-4:psipred_pos[L,H,E] 3 5-7:silents_pos[L,H,E]
		//out15:
		//type 1:oneSheet 2-44:psipred_pos[L,H,E] 45-89:silent_pos[L,H,E]
		//2:psipred-7 10:psipred@pos 17:psipred+7 18:silent-7(for sheet) 26:silents@pos 33::silent+7
		for(Size ii=1; ii<=native_pose.total_residue(); ++ii){
				Size type = 0;
				if(native_pose.secstruct(ii) == "H")
						type = 1;
				if(native_pose.secstruct(ii) == "L" 
						type = 2;
				if(native_pose.secstruct(ii) == "L" 
						type = 3;
				out1 << type << " ";
				out15 << type << " ";

				for (int kk=ii-range/2; kk<=

		
		


		


		 if(i_ss_dssp=="L")
              hle_type = 2
            else
              if(i_ss_dssp=="H")
                hle_type = 1
              else
                if(i_ss_dssp=="E")
                  hle_type = 3




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
