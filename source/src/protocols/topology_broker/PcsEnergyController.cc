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
/// @file protocols/topology_broker/PcsEnergyController.cc
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
/// @authorv Christophe Schmitz & Oliver Lange
///
////////////////////////////////////////////////

// Unit Headers
#include <protocols/topology_broker/PcsEnergyController.hh>

// Package Headers
#include <protocols/scoring/methods/pcs2/PcsEnergyParameterManager.hh>
#include <protocols/scoring/methods/pcs2/PcsGridSearchParameter.hh>
#include <protocols/scoring/methods/pcs2/PcsGridSearchParameterManager.hh>

// Project Headers

// Utility headers
#include <basic/Tracer.hh>

#include <utility/vector1.hh>
#include <basic/options/keys/OptionKeys.hh>


// ObjexxFCL Headers

// C++ headers

static THREAD_LOCAL basic::Tracer TR_PcsEnergyController( "protocols.topology_broker.PcsEnergyController" );

namespace protocols {
namespace topology_broker {

using namespace core;
//using namespace scoring::constraints;


PcsEnergyController::PcsEnergyController()
{
}

PcsEnergyController::PcsEnergyController(PcsEnergyController const &):
	TopologyClaimer()
{
}


PcsEnergyController::~PcsEnergyController()= default;

void
PcsEnergyController::set_defaults( ){
	using namespace protocols::scoring::methods::pcs2;
	PcsGridSearchParameterManager::get_instance()->incremente_n_multi_data();
}

bool PcsEnergyController::read_tag( std::string tag, std::istream& is ) {
	using namespace protocols::scoring::methods::pcs2;

	core::Size n_m_d;
	n_m_d = PcsGridSearchParameterManager::get_instance()->get_n_multi_data();

	PcsGridSearchParameter & g = PcsGridSearchParameterManager::get_instance()->get_grid_search_parameters(n_m_d);


	if ( tag == "PCS_WEIGHT" ) {
		if ( (is >>
				g.pcs_weight_stage1_ >>
				g.pcs_weight_stage2_ >>
				g.pcs_weight_stage3_ >>
				g.pcs_weight_stage4_).fail() ) {
			std::cerr << "problem while parsing the tag '" << tag << "' . Expecting 4 numbers (integer or float)" << std::endl;
			utility_exit();
		}
		return true;
	}


	if ( tag == "INCLUDE_ONLY" ) {
		if ( (is >>
				g.include_only_start_stage1_ >>
				g.include_only_end_stage1_
				).fail() ) {
			std::cerr << "problem while parsing the tag '" << tag << "' . Expecting 2 positive integer" << std::endl;
			utility_exit();
		}
		g.include_only_start_stage2_ = g.include_only_start_stage1_;
		g.include_only_start_stage3_ = g.include_only_start_stage1_;
		g.include_only_start_stage4_ = g.include_only_start_stage1_;

		g.include_only_end_stage2_ = g.include_only_end_stage1_;
		g.include_only_end_stage3_ = g.include_only_end_stage1_;
		g.include_only_end_stage4_ = g.include_only_end_stage1_;

		return true;
	}

	if ( tag == "INDIVIDUAL_SCALE" ) {
		if ( (is >>
				g.individual_scale_stage1_
				).fail() ) {
			std::cerr << "problem while parsing the tag '" << tag << "' . Expecting 1 positive number" << std::endl;
			utility_exit();
		}
		g.individual_scale_stage2_ = g.individual_scale_stage1_;
		g.individual_scale_stage3_ = g.individual_scale_stage1_;
		g.individual_scale_stage4_ = g.individual_scale_stage1_;

		return true;
	}

	if ( tag == "N_TRIALS_MINIMIZATION" ) {
		if ( (is >>
				g.n_trial_min_stage1_
				).fail() ) {
			std::cerr << "problem while parsing the tag '" << tag << "' . Expecting 1 positive integer" << std::endl;
			utility_exit();
		}
		g.n_trial_min_stage2_ = g.n_trial_min_stage1_;
		g.n_trial_min_stage3_ = g.n_trial_min_stage1_;
		g.n_trial_min_stage4_ = g.n_trial_min_stage1_;
		return true;
	}


	if ( tag == "PCS_INPUT_FILE" ) {
		core::Real weight;
		std::string filename;
		if ( (is >>
				filename >>
				weight
				).fail() ) {
			std::cerr << "problem while parsing the tag '" << tag << "' . Expecting 1 string followed by 1 number (integer or float)" << std::endl;
			utility_exit();
		}

		if ( weight < 0 ) {
			std::cerr << "problem while parsing the tag '" << tag << "' . Expecting 1 string followed by one POSITIVE weight" << std::endl;
			utility_exit();
		}
		g.filenames_.push_back(filename);
		g.individual_weights_.push_back(weight);
		return true;
	}

	std::cerr << "The following tag is unknown " << tag << std::endl;


	return Parent::read_tag( tag, is );
}

//This is called each time the stageID changed
void PcsEnergyController::add_mover(
	moves::RandomMover& /* random_mover */,
	core::pose::Pose const& /*pose*/,
	abinitio::StageID stageID , /* abinitio sampler stage */
	core::scoring::ScoreFunction const& /*scorefxn*/,
	core::Real /*progress  progress within stage */
)
{

	using namespace protocols::scoring::methods::pcs2;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	PcsGridSearchParameterManager * pcs_g_s_p_m = PcsGridSearchParameterManager::get_instance();

	//utility::vector1< Size > vec_exclude;
	//if ( option[ in::file::native_exclude_res].user() ) {
	// vec_exclude = option[ in::file::native_exclude_res ]();
	// PcsEnergyParameterManager::get_instance()->get_PCS_data_input_reference().set_vector_exclude_residues(vec_exclude);
	//}

	core::Size i_multi_data;
	core::Size n_m_d;

	n_m_d = pcs_g_s_p_m->get_n_multi_data();

	PcsEnergyParameterManager * pcs_e_p_m = PcsEnergyParameterManager::get_instance();

	for ( i_multi_data = 1;  i_multi_data <= n_m_d; ++i_multi_data ) {

		PcsGridSearchParameter & g = pcs_g_s_p_m->get_grid_search_parameters(i_multi_data);

		PcsEnergyParameter & pcs_e_p = pcs_e_p_m->get_PcsEnergyParameter_for(i_multi_data);

		if ( stageID == abinitio::STAGE_1 ) {

			pcs_e_p.set_grid_param(
				g.include_only_start_stage1_,
				g.include_only_end_stage1_,
				g.n_trial_min_stage1_,
				g.pcs_weight_stage1_,
				g.individual_scale_stage1_
			);
			TR_PcsEnergyController << "Switching STAGE_1; Grid parameters of center " << i_multi_data << " / " << n_m_d << std::endl;
			TR_PcsEnergyController << pcs_e_p;
		}

		if ( stageID == abinitio::STAGE_2 ) {


			pcs_e_p.set_grid_param(
				g.include_only_start_stage2_,
				g.include_only_end_stage2_,
				g.n_trial_min_stage2_,
				g.pcs_weight_stage2_,
				g.individual_scale_stage2_
			);
			TR_PcsEnergyController << "Switching STAGE_2; Grid parameters of center " << i_multi_data << " / " << n_m_d << std::endl;
			TR_PcsEnergyController << pcs_e_p;

		}

		if ( stageID == abinitio::STAGE_3a ) {

			pcs_e_p.set_grid_param(
				g.include_only_start_stage3_,
				g.include_only_end_stage3_,
				g.n_trial_min_stage3_,
				g.pcs_weight_stage3_,
				g.individual_scale_stage3_
			);
			TR_PcsEnergyController << "Switching STAGE_3a; Grid parameters of center " << i_multi_data << " / " << n_m_d << std::endl;
			TR_PcsEnergyController << pcs_e_p;

		}

		if ( stageID == abinitio::STAGE_3b ) {

			pcs_e_p.set_grid_param(
				g.include_only_start_stage3_,
				g.include_only_end_stage3_,
				g.n_trial_min_stage3_,
				g.pcs_weight_stage3_,
				g.individual_scale_stage3_
			);
			TR_PcsEnergyController << "Switching STAGE_3b; Grid parameters of center " << i_multi_data << " / " << n_m_d << std::endl;
			TR_PcsEnergyController << pcs_e_p;
		}

		if ( stageID == abinitio::STAGE_4 ) {

			pcs_e_p.set_grid_param(
				g.include_only_start_stage4_,
				g.include_only_end_stage4_,
				g.n_trial_min_stage4_,
				g.pcs_weight_stage4_,
				g.individual_scale_stage4_
			);
			TR_PcsEnergyController << "Switching STAGE_4; Grid parameters of center " << i_multi_data << " / " << n_m_d << std::endl;
			TR_PcsEnergyController << pcs_e_p;
		}
	} //for loop
}


//This is called after each PcsEnergyController CLAIMER is being read
void
PcsEnergyController::init_after_reading(){
	using namespace protocols::scoring::methods::pcs2;

	core::Size n_m_d;

	PcsEnergyParameterManager * pcs_e_p_m = PcsEnergyParameterManager::get_instance();

	n_m_d =  PcsGridSearchParameterManager::get_instance()->get_n_multi_data();

	PcsGridSearchParameter & g = PcsGridSearchParameterManager::get_instance()->get_grid_search_parameters(n_m_d);

	g.control_grid_param();

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	//utility::vector1< Size > vec_exclude;
	//if ( option[ in::file::native_exclude_res ].user() ) {
	// vec_exclude = option[ in::file::native_exclude_res ]();
	//PcsEnergyParameterManager::get_instance()get_PCS_data_input_reference().->set_vector_exclude_residues(vec_exclude);
	//}

	pcs_e_p_m->incremente_n_multi_data();

	pcs_e_p_m->get_PcsEnergyParameter_for(n_m_d).set_grid_param(
		g.include_only_start_stage4_,
		g.include_only_end_stage4_,
		g.n_trial_min_stage4_,
		g.pcs_weight_stage4_,
		g.individual_scale_stage4_
	);

	pcs_e_p_m->get_PcsEnergyParameter_for(n_m_d).set_vector_name_and_weight(g.filenames_, g.individual_weights_);

	/* core::Size k;
	for(k = 1; k <= n_m_d; k++){
	PcsGridSearchParameter & g_junk = PcsGridSearchParameterManager::get_instance()->get_grid_search_parameters(k);
	std::cerr << g_junk << std::endl;
	}
	*/
}

} //topology_broker
} //protocols
