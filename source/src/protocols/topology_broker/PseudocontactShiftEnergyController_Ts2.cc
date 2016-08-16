// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

//////////////////////////////////////////////
///
/// @file PseudocontactShiftEnergyController.cc
///
/// @brief
///
/// @details
///
/// @param
///
/// @return
///
/// @remarks
///
/// @references
///
/// @authorv Christophe Schmitz //kalabharath & Oliver Lange
///
////////////////////////////////////////////////

// Unit Headers
#include <protocols/topology_broker/PseudocontactShiftEnergyController_Ts2.hh>

// Package Headers
#include <protocols/scoring/methods/pcsTs2/PseudocontactShiftEnergy.hh>
// Project Headers
#include <core/scoring/ScoreType.hh>
#include <core/scoring/methods/EnergyMethod.hh>

//TEMP
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Utility headers
#include <basic/Tracer.hh>

// ObjexxFCL Headers

// C++ headers


static THREAD_LOCAL basic::Tracer tr_control( "protocols.topology_broker.PseudocontactShiftEnergyController_Ts2" );

namespace protocols {
namespace topology_broker {

using namespace core;
using namespace core::scoring::constraints;
PseudocontactShiftEnergyController_Ts2::PseudocontactShiftEnergyController_Ts2()
{
	grid_edge_stage1_ = 50.0;
	grid_edge_stage2_ = 50.0;
	grid_edge_stage3_ = 50.0;
	grid_edge_stage4_ = 50.0;

	grid_step_stage1_ = 3.0;
	grid_step_stage2_ = 3.0;
	grid_step_stage3_ = 3.0;
	grid_step_stage4_ = 3.0;

	grid_small_cutoff_stage1_ = 0.0;
	grid_small_cutoff_stage2_ = 0.0;
	grid_small_cutoff_stage3_ = 0.0;
	grid_small_cutoff_stage4_ = 0.0;

	grid_large_cutoff_stage1_ = 25.0;
	grid_large_cutoff_stage2_ = 25.0;
	grid_large_cutoff_stage3_ = 25.0;
	grid_large_cutoff_stage4_ = 25.0;

	grid_cone_angle_cutoff_stage1_ = 180.0;
	grid_cone_angle_cutoff_stage2_ = 180.0;
	grid_cone_angle_cutoff_stage3_ = 180.0;
	grid_cone_angle_cutoff_stage4_ = 180.0;

	grid_residue_num_1_stage1_ = -1;
	grid_residue_num_1_stage2_ = -1;
	grid_residue_num_1_stage3_ = -1;
	grid_residue_num_1_stage4_ = -1;

	grid_residue_num_2_stage1_ = -1;
	grid_residue_num_2_stage2_ = -1;
	grid_residue_num_2_stage3_ = -1;
	grid_residue_num_2_stage4_ = -1;

	grid_atom_name_1_stage1_ = "BAD";
	grid_atom_name_1_stage2_ = "BAD";
	grid_atom_name_1_stage3_ = "BAD";
	grid_atom_name_1_stage4_ = "BAD";

	grid_atom_name_2_stage1_ = "BAD";
	grid_atom_name_2_stage2_ = "BAD";
	grid_atom_name_2_stage3_ = "BAD";
	grid_atom_name_2_stage4_ = "BAD";

	grid_k_vector_stage1_ = 0.0;
	grid_k_vector_stage2_ = 0.0;
	grid_k_vector_stage3_ = 0.0;
	grid_k_vector_stage4_ = 0.0;

	minimize_best_tensor_stage1_ = false;
	minimize_best_tensor_stage2_ = false;
	minimize_best_tensor_stage3_ = false;
	minimize_best_tensor_stage4_ = false;

	pcs_weight_stage1_ = 10.0;
	pcs_weight_stage2_ = 10.0;
	pcs_weight_stage3_ = 10.0;
	pcs_weight_stage4_ = 10.0;
}
void
PseudocontactShiftEnergyController_Ts2::set_defaults( ){

}

bool PseudocontactShiftEnergyController_Ts2::read_tag( std::string tag, std::istream& is ) {
	using namespace protocols::scoring::methods::pcsTs2;

	if ( tag == "TS2_GRID_EDGE_SIZE" ) {
		if ( (is >>
				grid_edge_stage1_ >>
				grid_edge_stage2_ >>
				grid_edge_stage3_ >>
				grid_edge_stage4_).fail() ) {
			std::cerr << "problem while parsing the tag '" << tag << "' . Expecting 4 numbers (integer or float)" << std::endl;
			utility_exit();
		}
		return true;
	}

	if ( tag == "TS2_GRID_STEP_SIZE" ) {
		if ( (is >>
				grid_step_stage1_ >>
				grid_step_stage2_ >>
				grid_step_stage3_ >>
				grid_step_stage4_).fail() ) {
			std::cerr << "problem while parsing the tag '" << tag << "' . Expecting 4 numbers (integer or float)" << std::endl;
			utility_exit();
		}
		return true;
	}

	if ( tag == "TS2_GRID_SMALL_CUTOFF" ) {
		if ( (is >>
				grid_small_cutoff_stage1_ >>
				grid_small_cutoff_stage2_ >>
				grid_small_cutoff_stage3_ >>
				grid_small_cutoff_stage4_).fail() ) {
			std::cerr << "problem while parsing the tag '" << tag << "' . Expecting 4 numbers (integer or float)" << std::endl;
			utility_exit();
		}
		return true;
	}

	if ( tag == "TS2_GRID_LARGE_CUTOFF" ) {
		if ( (is >>
				grid_large_cutoff_stage1_ >>
				grid_large_cutoff_stage2_ >>
				grid_large_cutoff_stage3_ >>
				grid_large_cutoff_stage4_).fail() ) {
			std::cerr << "problem while parsing the tag '" << tag << "' . Expecting 4 numbers (integer or float)" << std::endl;
			utility_exit();
		}
		return true;
	}

	if ( tag == "TS2_GRID_CONE_ANGLE_CUTOFF" ) {
		if ( (is >>
				grid_cone_angle_cutoff_stage1_ >>
				grid_cone_angle_cutoff_stage2_ >>
				grid_cone_angle_cutoff_stage3_ >>
				grid_cone_angle_cutoff_stage4_).fail() ) {
			std::cerr << "problem while parsing the tag '" << tag << "' . Expecting 4 numbers (integer or float)" << std::endl;
			utility_exit();
		}
		return true;
	}

	if ( tag == "TS2_VECTOR_ATOM1_NAME" ) {
		if ( (is >>
				grid_atom_name_1_stage1_ >>
				grid_atom_name_1_stage2_ >>
				grid_atom_name_1_stage3_ >>
				grid_atom_name_1_stage4_).fail() ) {
			std::cerr << "problem while parsing the tag '" << tag << "' . Expecting 4 strings" << std::endl;
			utility_exit();
		}
		return true;
	}

	if ( tag == "TS2_VECTOR_ATOM2_NAME" ) {
		if ( (is >>
				grid_atom_name_2_stage1_ >>
				grid_atom_name_2_stage2_ >>
				grid_atom_name_2_stage3_ >>
				grid_atom_name_2_stage4_).fail() ) {
			std::cerr << "problem while parsing the tag '" << tag << "' . Expecting 4 strings" << std::endl;
			utility_exit();
		}
		return true;
	}

	if ( tag == "TS2_VECTOR_ATOM1_RESIDUE_NUM" ) {
		if ( (is >>
				grid_residue_num_1_stage1_ >>
				grid_residue_num_1_stage2_ >>
				grid_residue_num_1_stage3_ >>
				grid_residue_num_1_stage4_).fail() ) {
			std::cerr << "problem while parsing the tag '" << tag << "' . Expecting 4 integers" << std::endl;
			utility_exit();
		}
		return true;
	}

	if ( tag == "TS2_VECTOR_ATOM2_RESIDUE_NUM" ) {
		if ( (is >>
				grid_residue_num_2_stage1_ >>
				grid_residue_num_2_stage2_ >>
				grid_residue_num_2_stage3_ >>
				grid_residue_num_2_stage4_).fail() ) {
			std::cerr << "problem while parsing the tag '" << tag << "' . Expecting 4 integers" << std::endl;
			utility_exit();
		}
		return true;
	}

	if ( tag == "TS2_K_VECTOR" ) {
		if ( (is >>
				grid_k_vector_stage1_ >>
				grid_k_vector_stage2_ >>
				grid_k_vector_stage3_ >>
				grid_k_vector_stage4_).fail() ) {
			std::cerr << "problem while parsing the tag '" << tag << "' . Expecting 4 numbers (integer or float)" << std::endl;
			utility_exit();
		}
		return true;
	}

	if ( tag == "TS2_MINIMIZE_BEST_TENSOR" ) {
		core::Size a, b, c, d;
		if ( (is >> a >> b >>c >> d).fail() ) {
			std::cerr << "problem while parsing the tag '" << tag << "' . Expecting 4 integers between 0 and 1" << std::endl;
			utility_exit();
		}
		if ( ((a != 0) && (a != 1)) ||
				((b != 0) && (b != 1)) ||
				((c != 0) && (c != 1)) ||
				((d != 0) && (d != 1)) ) {
			std::cerr << "problem while parsing the tag '" << tag << "' . Expecting 4 integers between 0 and 1" << std::endl;
			utility_exit();
		}
		minimize_best_tensor_stage1_ = (a == 1);
		minimize_best_tensor_stage2_ = (b == 1);
		minimize_best_tensor_stage3_ = (c == 1);
		minimize_best_tensor_stage4_ = (d == 1);
		return true;
	}

	if ( tag == "TS2_PCS_WEIGHT" ) {
		if ( (is >>
				pcs_weight_stage1_ >>
				pcs_weight_stage2_ >>
				pcs_weight_stage3_ >>
				pcs_weight_stage4_).fail() ) {
			std::cerr << "problem while parsing the tag '" << tag << "' . Expecting 4 numbers (integer or float)" << std::endl;
			utility_exit();
		}
		return true;
	}


	/*
	the weighting scheme is not working properly for the moment.
	if ( tag == "PCS_INPUT_FILE" ) {
	core::Real weight;
	std::string filename;
	if ((is >>
	filename >>
	weight).fail()){
	std::cerr << "problem while parsing the tag '" << tag << "' . Expecting 1 string followed by 1 number (integer or float) positif" << std::endl;
	utility_exit();
	}
	if(weight < 0){
	std::cerr << "problem while parsing the tag '" << tag << "' . Expecting 1 string followed by 1 number (integer or float) positif" << std::endl;
	utility_exit();
	}
	filenames_.push_back(filename);
	individual_weights_.push_back(weight);
	return true;
	}
	*/


	if ( tag == "TS2_PCS_INPUT_FILE" ) {
		core::Real weight;
		std::string filename;
		if ( (is >>
				filename
				).fail() ) {
			std::cerr << "problem while parsing the tag '" << tag << "' . Expecting 1 string" << std::endl;
			utility_exit();
		}
		//TODO: work with weigth scheme for individual lanthanide? this is temporary.
		weight = 1.0;

		if ( weight < 0 ) {
			std::cerr << "problem while parsing the tag '" << tag << "' . Expecting 1 string" << std::endl;
			utility_exit();
		}
		filenames_.push_back(filename);
		individual_weights_.push_back(weight);
		return true;
	}

	std::cerr << "The following tag is unknown " << tag << std::endl;


	return Parent::read_tag( tag, is );
}

void PseudocontactShiftEnergyController_Ts2::add_mover(
	moves::RandomMover& /* random_mover */,
	core::pose::Pose const& /*pose*/,
	abinitio::StageID stageID , /* abinitio sampler stage */
	core::scoring::ScoreFunction const& /*scorefxn*/,
	core::Real /*progress  progress within stage */
)
{

	using namespace protocols::scoring::methods::pcsTs2;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	//utility::vector1< Size > vec_exclude;
	//if ( option[ in::file::native_exclude_res].user() ) {
	// vec_exclude = option[ in::file::native_exclude_res ]();
	// PCS_Energy_parameters_manager_Ts2::get_instance()->set_vector_exclude_residues(vec_exclude);
	//}


	if ( stageID == abinitio::STAGE_1 ) {

		PCS_Energy_parameters_manager_Ts2::get_instance()->set_grid_param(grid_edge_stage1_,
			grid_step_stage1_,
			grid_small_cutoff_stage1_,
			grid_large_cutoff_stage1_,
			grid_cone_angle_cutoff_stage1_,
			grid_atom_name_1_stage1_,
			grid_atom_name_2_stage1_,
			grid_residue_num_1_stage1_,
			grid_residue_num_2_stage1_,
			grid_k_vector_stage1_,
			minimize_best_tensor_stage1_,
			pcs_weight_stage1_
		);
	}

	if ( stageID == abinitio::STAGE_2 ) {


		PCS_Energy_parameters_manager_Ts2::get_instance()->set_grid_param(grid_edge_stage2_,
			grid_step_stage2_,
			grid_small_cutoff_stage2_,
			grid_large_cutoff_stage2_,
			grid_cone_angle_cutoff_stage2_,
			grid_atom_name_1_stage2_,
			grid_atom_name_2_stage2_,
			grid_residue_num_1_stage2_,
			grid_residue_num_2_stage2_,
			grid_k_vector_stage2_,
			minimize_best_tensor_stage2_,
			pcs_weight_stage2_
		);
	}

	if ( stageID == abinitio::STAGE_3a ) {

		PCS_Energy_parameters_manager_Ts2::get_instance()->set_grid_param(grid_edge_stage3_,
			grid_step_stage3_,
			grid_small_cutoff_stage3_,
			grid_large_cutoff_stage3_,
			grid_cone_angle_cutoff_stage3_,
			grid_atom_name_1_stage3_,
			grid_atom_name_2_stage3_,
			grid_residue_num_1_stage3_,
			grid_residue_num_2_stage3_,
			grid_k_vector_stage3_,
			minimize_best_tensor_stage3_,
			pcs_weight_stage3_
		);
	}

	if ( stageID == abinitio::STAGE_3b ) {

		PCS_Energy_parameters_manager_Ts2::get_instance()->set_grid_param(grid_edge_stage3_,
			grid_step_stage3_,
			grid_small_cutoff_stage3_,
			grid_large_cutoff_stage3_,
			grid_cone_angle_cutoff_stage3_,
			grid_atom_name_1_stage3_,
			grid_atom_name_2_stage3_,
			grid_residue_num_1_stage3_,
			grid_residue_num_2_stage3_,
			grid_k_vector_stage3_,
			minimize_best_tensor_stage3_,
			pcs_weight_stage3_
		);
	}

	if ( stageID == abinitio::STAGE_4 ) {

		PCS_Energy_parameters_manager_Ts2::get_instance()->set_grid_param(grid_edge_stage4_,
			grid_step_stage4_,
			grid_small_cutoff_stage4_,
			grid_large_cutoff_stage4_,
			grid_cone_angle_cutoff_stage4_,
			grid_atom_name_1_stage4_,
			grid_atom_name_2_stage4_,
			grid_residue_num_1_stage4_,
			grid_residue_num_2_stage4_,
			grid_k_vector_stage4_,
			minimize_best_tensor_stage4_,
			pcs_weight_stage4_
		);
	}
}


void
PseudocontactShiftEnergyController_Ts2::init_after_reading(){
	using namespace protocols::scoring::methods::pcsTs2;
	control_grid_param();

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	//utility::vector1< Size > vec_exclude;
	//if ( option[ in::file::native_exclude_res ].user() ) {
	// vec_exclude = option[ in::file::native_exclude_res ]();
	//PCS_Energy_parameters_manager_Ts2::get_instance()->set_vector_exclude_residues(vec_exclude);
	//}


	PCS_Energy_parameters_manager_Ts2::get_instance()->set_grid_param(grid_edge_stage4_,
		grid_step_stage4_,
		grid_small_cutoff_stage4_,
		grid_large_cutoff_stage4_,
		grid_cone_angle_cutoff_stage4_,
		grid_atom_name_1_stage4_,
		grid_atom_name_2_stage4_,
		grid_residue_num_1_stage4_,
		grid_residue_num_2_stage4_,
		grid_k_vector_stage4_,
		minimize_best_tensor_stage4_,
		pcs_weight_stage4_
	);

	PCS_Energy_parameters_manager_Ts2::get_instance()->set_vector_name_and_weight(filenames_, individual_weights_);

}

void
PseudocontactShiftEnergyController_Ts2::control_grid_param(){

	if ( (grid_edge_stage1_ <= 0)||
			(grid_edge_stage2_ <= 0)||
			(grid_edge_stage3_ <= 0)||
			(grid_edge_stage4_ <= 0)||

			(grid_step_stage1_ <= 0)||
			(grid_step_stage2_ <= 0)||
			(grid_step_stage3_ <= 0)||
			(grid_step_stage4_ <= 0)||

			(grid_edge_stage1_ < grid_step_stage1_)||
			(grid_edge_stage2_ < grid_step_stage2_)||
			(grid_edge_stage3_ < grid_step_stage3_)||
			(grid_edge_stage4_ < grid_step_stage4_)||

			(grid_small_cutoff_stage1_  < 0)||
			(grid_small_cutoff_stage2_  < 0)||
			(grid_small_cutoff_stage3_  < 0)||
			(grid_small_cutoff_stage4_  < 0)||

			(grid_large_cutoff_stage1_ <= 0)||
			(grid_large_cutoff_stage2_ <= 0)||
			(grid_large_cutoff_stage3_ <= 0)||
			(grid_large_cutoff_stage4_ <= 0)||

			(grid_large_cutoff_stage1_ <= grid_small_cutoff_stage1_)||
			(grid_large_cutoff_stage2_ <= grid_small_cutoff_stage2_)||
			(grid_large_cutoff_stage3_ <= grid_small_cutoff_stage3_)||
			(grid_large_cutoff_stage4_ <= grid_small_cutoff_stage4_)||

			(grid_cone_angle_cutoff_stage1_ < 0)||
			(grid_cone_angle_cutoff_stage2_ < 0)||
			(grid_cone_angle_cutoff_stage3_ < 0)||
			(grid_cone_angle_cutoff_stage4_ < 0)||
			(grid_cone_angle_cutoff_stage1_ > 180)||
			(grid_cone_angle_cutoff_stage2_ > 180)||
			(grid_cone_angle_cutoff_stage3_ > 180)||
			(grid_cone_angle_cutoff_stage4_ > 180)||

			(pcs_weight_stage1_ < 0)||
			(pcs_weight_stage2_ < 0)||
			(pcs_weight_stage3_ < 0)||
			(pcs_weight_stage4_ < 0)
			) {

		std::cerr << "In the definition of the grid search parameters for the PseudocontactShift, you should ensure that for all stages:" << std::endl;
		std::cerr << "1. the GRID_EDGE_SIZE is strictly positif" << std::endl;
		std::cerr << "2. the GRID_STEP_SIZE is strictly positif" << std::endl;
		std::cerr << "3. the GRID_STEP_SIZE is smaller than the GRID_EDGE_SIZE" << std::endl;
		std::cerr << "4. the GRID_SMALL_CUTOFF is positif" << std::endl;
		std::cerr << "5. the GRID_LARGE_CUTOFF is strictly positif" << std::endl;
		std::cerr << "6. the GRID_LARGE_CUTOFF is strictly larger than the GRID_SMALL_CUTOFF" << std::endl;
		std::cerr << "7. the GRID_CONE_ANGLE_CUTOFF is between 0 and 180" << std::endl;
		std::cerr << "8. the PCS_WEIGHT is positif" << std::endl;
		utility_exit_with_message("Please, review your PseudocontactShift grid search parameter file and correct it");
	}

	if ( (grid_atom_name_1_stage1_ == "BAD")||
			(grid_atom_name_1_stage2_ == "BAD")||
			(grid_atom_name_1_stage3_ == "BAD")||
			(grid_atom_name_1_stage4_ == "BAD")||

			(grid_atom_name_2_stage1_  == "BAD")||
			(grid_atom_name_2_stage2_  == "BAD")||
			(grid_atom_name_2_stage3_  == "BAD")||
			(grid_atom_name_2_stage4_  == "BAD")||

			(grid_residue_num_1_stage1_ == -1)||
			(grid_residue_num_1_stage2_ == -1)||
			(grid_residue_num_1_stage3_ == -1)||
			(grid_residue_num_1_stage4_ == -1)||

			(grid_residue_num_2_stage1_ == -1)||
			(grid_residue_num_2_stage2_ == -1)||
			(grid_residue_num_2_stage3_ == -1)||
			(grid_residue_num_2_stage4_ == -1) ) {
		std::cerr << "In the definition of the grid search parameters for the PseudocontactShift, you must define the four following flags:" << std::endl;
		std::cerr << "VECTOR_ATOM1_NAME" << std::endl;
		std::cerr << "VECTOR_ATOM2_NAME" << std::endl;
		std::cerr << "VECTOR_ATOM1_RESIDUE_NUM" << std::endl;
		std::cerr << "VECTOR_ATOM1_RESIDUE_NUM" << std::endl;
		utility_exit_with_message("Please, review your PseudocontactShift grid search parameter file and correct it");
	}
}
} //topology_broker
} //protocols
