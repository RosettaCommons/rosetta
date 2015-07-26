// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/init/score_function_corrections.cc
/// @brief  Initialize Score Funciton Corrections
/// @author Matthew O'Meara


// Unit headers
#include <core/init/score_function_corrections.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/cryst.OptionKeys.gen.hh>
#include <basic/options/keys/ddg.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/mistakes.OptionKeys.gen.hh>
#include <basic/options/keys/qsar.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/optimization.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>

// Core headers
#include <core/types.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>

// C++ headers
#include <istream>
#include <string>

using namespace basic::options;
using namespace basic::options::OptionKeys;

//MaximCode:
static thread_local basic::Tracer TR( "core.init.score_function_corrections" );

namespace core {
namespace init {

struct pre_talaris_2013_behavior_settings {
	pre_talaris_2013_behavior_settings();

	bool pre_talaris2013_geometries;
	bool hb_sp2_chipen;
	bool hb_fade_energy;
	bool hbond_measure_sp3acc_BAH_from_hvy;
	Real lj_hbond_hdis;
	Real lj_hbond_OH_donor_dis;
	Real hb_sp2_outer_width;
	Real expand_st_chi2sampling;
	std::string score_weights;
	std::string score_patch;
	bool analytic_etable_evaluation;
	std::string hbond_params;
	bool smooth_fa_elec;
	Real elec_min_dis;
	bool dun10;
	bool use_bicubic_interpolation;
};

pre_talaris_2013_behavior_settings::pre_talaris_2013_behavior_settings() :
	pre_talaris2013_geometries( true ),
	hb_sp2_chipen( false ),
	hb_fade_energy( false ),
	hbond_measure_sp3acc_BAH_from_hvy( false ),
	lj_hbond_hdis( 1.95 ),
	lj_hbond_OH_donor_dis( 3.0 ),
	hb_sp2_outer_width( 0.33333 ),
	expand_st_chi2sampling( false ),
	score_weights( "pre_talaris_2013_standard.wts" ),
	score_patch( "" ), // the default patch to score12 is handled in the core::scoring::get_score_function
	analytic_etable_evaluation( false ),
	hbond_params( "score12_params" ),
	smooth_fa_elec( false ),
	elec_min_dis( 1.5 ),
	dun10( false ),
	use_bicubic_interpolation( false )
{}

pre_talaris_2013_behavior_settings const restore_sc12_settings;

static thread_local basic::Tracer TR( "core.init.score_function_corrections" );

void
init_revert_to_pre_talaris_2013_mistake() {
	if ( option[ mistakes::restore_pre_talaris_2013_behavior ] ) {
		revert_to_pre_talaris_2013_defaults();
	}
}

void revert_to_pre_talaris_2013_defaults() {
	if ( ! option[ mistakes::chemical::pre_talaris2013_geometries ].user() ) {
		option[ mistakes::chemical::pre_talaris2013_geometries ].default_value( restore_sc12_settings.pre_talaris2013_geometries );
	}
	if ( ! option[ corrections::score::hb_sp2_chipen ].user() ) {
		option[ corrections::score::hb_sp2_chipen ].default_value( restore_sc12_settings.hb_sp2_chipen );
	}
	if ( ! option[ corrections::score::hb_fade_energy ].user() ) {
		option[ corrections::score::hb_fade_energy ].default_value( restore_sc12_settings.hb_fade_energy );
	}
	if ( ! option[ corrections::score::hbond_measure_sp3acc_BAH_from_hvy ].user() ) {
		option[ corrections::score::hbond_measure_sp3acc_BAH_from_hvy ].default_value( restore_sc12_settings.hbond_measure_sp3acc_BAH_from_hvy );
	}
	if ( ! option[ corrections::score::lj_hbond_hdis ].user() ) {
		option[ corrections::score::lj_hbond_hdis ].default_value( restore_sc12_settings.lj_hbond_hdis );
	}
	if ( ! option[ corrections::score::lj_hbond_OH_donor_dis ].user() ) {
		option[ corrections::score::lj_hbond_OH_donor_dis ].default_value( restore_sc12_settings.lj_hbond_OH_donor_dis );
	}
	if ( ! option[ corrections::score::hb_sp2_outer_width ].user() ) {
		option[ corrections::score::hb_sp2_outer_width ].default_value( restore_sc12_settings.hb_sp2_outer_width );
	}
	if ( ! option[ corrections::chemical::expand_st_chi2sampling ].user() ) {
		option[ corrections::chemical::expand_st_chi2sampling ].default_value( restore_sc12_settings.expand_st_chi2sampling );
	}
	if ( ! option[ score::weights ].user() ) {
		option[ score::weights ].default_value( restore_sc12_settings.score_weights );
		/// Don't set the score12 weights patch if the user has provided a value for either score::weights or score::patch
		if ( ! option[ score::patch ].user() ) {
			utility::vector1< std::string > default_score_patch(1);
			default_score_patch[1] = restore_sc12_settings.score_patch;
			option[ score::patch ].default_value( default_score_patch );
		}
	}
	if ( ! option[ score::analytic_etable_evaluation ].user() ) {
		option[ score::analytic_etable_evaluation ].default_value( restore_sc12_settings.analytic_etable_evaluation );
	}
	if ( ! option[ score::hbond_params ].user() ) {
		option[ score::hbond_params ].default_value( restore_sc12_settings.hbond_params );
	}
	if ( ! option[ score::smooth_fa_elec ].user() ) {
		option[ score::smooth_fa_elec ].default_value( restore_sc12_settings.smooth_fa_elec );
	}
	if ( ! option[ score::elec_min_dis ].user() ) {
		option[ score::elec_min_dis ].default_value( restore_sc12_settings.elec_min_dis );
	}
	if( ! option[ basic::options::OptionKeys::chemical::set_atom_properties ].user() ) {
		utility::vector1< std::string > params;
		params.push_back("fa_standard:ONH2:LK_DGFREE:-10.0000");
		params.push_back("fa_standard:NH2O:LK_DGFREE:-10.0000");
		params.push_back("fa_standard:Narg:LK_DGFREE:-11.0000");
		params.push_back("fa_standard:OH:LK_DGFREE:-6.7700" );
		option[ basic::options::OptionKeys::chemical::set_atom_properties ].value(params);
	}

	if ( ! option[ corrections::score::dun10 ].user() ) {
		option[corrections::score::dun10].default_value( restore_sc12_settings.dun10 );
	}

	if ( ! option[ corrections::score::use_bicubic_interpolation ].user() ) {
		option[corrections::score::use_bicubic_interpolation].default_value( restore_sc12_settings.use_bicubic_interpolation );
	}
	if ( ! option[ qsar::weights ].user() ) {
		option[ qsar::weights ].default_value( restore_sc12_settings.score_weights );
	}
	if ( ! option[ ddg::min_cst_weights ].user() ) {
		option[ ddg::min_cst_weights ].default_value( restore_sc12_settings.score_weights );
	}
	if ( ! option[ score::pack_weights ].user() ) {
		option[ score::pack_weights ].default_value( restore_sc12_settings.score_weights );
  }
}


void
init_hbond_sp2_correction() {
	if( option[ corrections::hbond_sp2_correction ]) {
		TR.Warning << "The -hbond_sp2_correction is deprecated.  The default talaris2013 score function should be relied upon instead." << std::endl;


		// Note, the weight sets that go with the hbond_sp2_correction
		// flag are handled in core/scoring/ScorefunctionFactory.hh

		if( ! option[ corrections::score::hb_sp2_chipen ].user() ) {
			option[ corrections::score::hb_sp2_chipen ].value( true );
		}

		if( ! option[ corrections::score::hb_sp2_BAH180_rise ].user() ) {
			option[ corrections::score::hb_sp2_BAH180_rise ].value( 0.75 );
		}

		if( ! option[ corrections::score::hb_sp2_outer_width ].user() ) {
			option[ corrections::score::hb_sp2_outer_width ].value( 0.357 );
		}

		if( ! option[ corrections::score::hb_fade_energy ].user() ) {
			option[ corrections::score::hb_fade_energy ].value( true );
		}

		if( ! option[ corrections::score::hbond_measure_sp3acc_BAH_from_hvy ].user() ) {
			option[ corrections::score::hbond_measure_sp3acc_BAH_from_hvy ].value( true );
		}

		if( ! option[ corrections::score::lj_hbond_hdis ].user() ) {
			option[ corrections::score::lj_hbond_hdis ].value( 1.75 );
		}

		if( ! option[ corrections::score::lj_hbond_OH_donor_dis ].user() ) {
			option[ corrections::score::lj_hbond_OH_donor_dis ].value( 2.6 );
		}

		if ( ! option[ corrections::chemical::expand_st_chi2sampling ].user() ) {
			option[ corrections::chemical::expand_st_chi2sampling ].value( true );
		}

		if( ! option[ score::hbond_params ].user() ) {
			option[ score::hbond_params ].value( "sp2_elec_params" );
		}

		if( ! option[ score::smooth_fa_elec ].user() ) {
			option[ score::smooth_fa_elec ].value( true );
		}

		if( ! option[ score::elec_min_dis ].user() ) {
			option[ score::elec_min_dis ].value( 1.6 );
		}

		if( ! option[ score::elec_r_option ].user() ) {
			option[ score::elec_r_option ].value( false );
		}

		if( ! option[ basic::options::OptionKeys::chemical::set_atom_properties ].user() || option[ mistakes::restore_pre_talaris_2013_behavior ] ) {
			std::cout << "overriding sc12 changes to the LK_DGFREEs" << std::endl;
			utility::vector1< std::string > params;
			params.push_back("fa_standard:ONH2:LK_DGFREE:-5.85");
			params.push_back("fa_standard:NH2O:LK_DGFREE:-7.8");
			params.push_back("fa_standard:Narg:LK_DGFREE:-10.0");
			params.push_back("fa_standard:OH:LK_DGFREE:-6.70" );
			option[ basic::options::OptionKeys::chemical::set_atom_properties ].value(params);
		}

		if ( ! option[ corrections::score::dun10 ].user() ) {
			option[corrections::score::dun10].value( true );
		}

		if ( ! option[ corrections::score::use_bicubic_interpolation ].user() ) {
			option[corrections::score::use_bicubic_interpolation].value( true );
		}

		if ( ! option[ score::weights ].user() ) {
			option[ score::weights ].value( "sp2_correction" );
		}

		if( option[corrections::correct] ){
			throw utility::excn::EXCN_BadInput( "The -corrections::hbond_sp2_correction is incompatible with the -corrections:correct flag.");
		}

 		if( option[ corrections::facts_default ]) {
 			throw utility::excn::EXCN_BadInput( "The -corrections::hbond_sp2_correction is incompatible with the -corrections:facts_default flag.");
		}
	}
}

void
init_facts_correction() {
 	// corrections specific for FACTS: this has to precede -correct
 	if( option[corrections::facts_default]) {
 		// Call all the corrections from -correct
 		option[corrections::correct].value(true);

 		/// Below are FACTS specific flags
 		// Hbond LJ parameters
		if ( ! option[ corrections::score::lj_hbond_hdis ].user() ) {
			option[corrections::score::lj_hbond_hdis].value(2.3);
		}

		if ( ! option[ corrections::score::lj_hbond_OH_donor_dis ].user() ) {
			option[corrections::score::lj_hbond_OH_donor_dis].value(3.4);
		}

 		// FACTS options - these are not fully optimized, can be changed in the future
 		if ( ! option[score::facts_GBpair_cut].user() )
			option[score::facts_GBpair_cut].value(10.0);

 		if ( ! option[score::facts_dshift].user() ){
			utility::vector1< Real > params;
			params.push_back( 0.0 );
			params.push_back( 1.5 );
			params.push_back( 1.5 );
			params.push_back( 1.5 );
			option[score::facts_dshift].value( params );
		}

 		if ( ! option[score::facts_die].user() )
			option[score::facts_die].value(1.0);

 		if ( ! option[score::facts_kappa].user() )
			option[score::facts_kappa].value(12.0);

 		if ( ! option[score::facts_plane_to_self].user() )
			option[score::facts_plane_to_self].value(true);

 		if ( ! option[score::facts_asp_patch].user() )
			option[score::facts_asp_patch].value(3);

 		// Below are default for sp2correction which are FACTS defaults too
 		if ( ! option[ corrections::score::use_bicubic_interpolation ].user() )
			option[corrections::score::use_bicubic_interpolation].value(true);

 		if ( ! option[score::hbond_params].user() )
			option[score::hbond_params ].value( "sp2_params" );

 		if ( ! option[corrections::score::hb_sp2_chipen].user() )
			option[corrections::score::hb_sp2_chipen].value(true);

 		if ( ! option[corrections::score::hbond_measure_sp3acc_BAH_from_hvy].user() )
			option[corrections::score::hbond_measure_sp3acc_BAH_from_hvy].value(true);

		if( ! option[ corrections::score::hb_fade_energy ].user() )
			option[ corrections::score::hb_fade_energy ].value( true );

 	} // end facts_default
}

//MaximCode:
void
init_correct_correction() {

	// set default corrections
	if( option[corrections::correct]) {

		/// Legacy warning
		if ( ! option[ mistakes::restore_pre_talaris_2013_behavior ] ) {
			TR.Warning << "-correct does not behave the way it did before talaris2013 became default; consider adding the -restore_pre_talaris_2013_behavior flag" << std::endl;
		}

		// Pair energy
		if ( ! option[ corrections::score::no_his_his_pairE ].user() ) {
			option[corrections::score::no_his_his_pairE].value( true );
		}

		// p_aa_pp
		if ( ! option[ corrections::score::p_aa_pp ].user() ) {
			option[corrections::score::p_aa_pp].value( "scoring/score_functions/P_AA_pp/P_AA_pp_08.2009" );
		}
		if ( ! option[ corrections::score::p_aa_pp_nogridshift ].user() ) {
			option[corrections::score::p_aa_pp_nogridshift].value( true );
		}

		//Ramachandran
		if ( ! option[ corrections::score::rama_not_squared ].user() ) {
			option[corrections::score::rama_not_squared].value( true );
		}
		if ( ! option[ corrections::score::rama_map ].user() ) {
			option[ corrections::score::rama_map ].value("scoring/score_functions/rama/Rama09_noEH_kernel25_it08.dat");
		}

		//rotamer library
		if ( ! option[ corrections::score::dun02_file ].user() ) {
			option[corrections::score::dun02_file].value( "rotamer/bbdep02.May.sortlib-correct.12.2010" );
		}

		//icoor
		if ( ! option[ corrections::chemical::icoor_05_2009 ].user() ) {
			option[corrections::chemical::icoor_05_2009].value( true );
		}
		if ( ! option[ score::hbond_params ].user() ) {
			option[score::hbond_params].value( "correct_params" );
		}

		if ( ! option[ corrections::score::ch_o_bond_potential ].user() ) {
			option[ corrections::score::ch_o_bond_potential ].value("scoring/score_functions/carbon_hbond/ch_o_bond_potential_near_min_yf.dat");
		}

		if ( ! option[ score::weights ].user() ) {
			option[ score::weights ].value( "score12_w_corrections" );
		}
	}
}

//MaximCode:
void
init_shapovalov_lib_fixes_enable_correction()
{
	// set default corrections
	if( option[corrections::shapovalov_lib_fixes_enable])
	{
		if (option[corrections::shapovalov_lib::shap_dun10_enable])
		{
			if ( ! option[ corrections::shapovalov_lib::shap_dun10_dir ].user() )
			{
				std::string _smoothingRequsted = option[ corrections::shapovalov_lib::shap_dun10_smooth_level ];
				if (_smoothingRequsted.compare("1")==0 || _smoothingRequsted.compare("lowest_smooth")==0)
				{
					option[ corrections::shapovalov_lib::shap_dun10_dir ].value("rotamer/shapovalov/StpDwn_0-0-0");
					option[ corrections::shapovalov_lib::shap_dun10_smooth_level ].value("1");
				}
				else if (_smoothingRequsted.compare("2")==0 || _smoothingRequsted.compare("lower_smooth")==0)
				{
					option[ corrections::shapovalov_lib::shap_dun10_dir ].value("rotamer/shapovalov/StpDwn_2-2-2");
					option[ corrections::shapovalov_lib::shap_dun10_smooth_level ].value("2");
				}
				else if (_smoothingRequsted.compare("3")==0 || _smoothingRequsted.compare("low_smooth")==0)
				{
					option[ corrections::shapovalov_lib::shap_dun10_dir ].value("rotamer/shapovalov/StpDwn_5-5-5");
					option[ corrections::shapovalov_lib::shap_dun10_smooth_level ].value("3");
				}
				else if (_smoothingRequsted.compare("4")==0 || _smoothingRequsted.compare("average_smooth")==0)
				{
					option[ corrections::shapovalov_lib::shap_dun10_dir ].value("rotamer/shapovalov/StpDwn_10-10-10");
					option[ corrections::shapovalov_lib::shap_dun10_smooth_level ].value("4");
				}
				else if (_smoothingRequsted.compare("5")==0 || _smoothingRequsted.compare("higher_smooth")==0)
				{
					option[ corrections::shapovalov_lib::shap_dun10_dir ].value("rotamer/shapovalov/StpDwn_20-20-20");
					option[ corrections::shapovalov_lib::shap_dun10_smooth_level ].value("5");
				}
				else if (_smoothingRequsted.compare("6")==0 || _smoothingRequsted.compare("highest_smooth")==0)
				{
					option[ corrections::shapovalov_lib::shap_dun10_dir ].value("rotamer/shapovalov/StpDwn_25-25-25");
					option[ corrections::shapovalov_lib::shap_dun10_smooth_level ].value("6");
				}
				else
				{
					option[ corrections::shapovalov_lib::shap_dun10_dir ].value("rotamer/shapovalov/StpDwn_5-5-5");
					option[ corrections::shapovalov_lib::shap_dun10_smooth_level ].value("3");
				}
			}
			else
			{
				option[ corrections::shapovalov_lib::shap_dun10_smooth_level ].value("overriden");
			}
		}


		if (option[corrections::shapovalov_lib::shap_rama_enable])
		{
			if ( ! option[ corrections::shapovalov_lib::shap_rama_map ].user() )
			{
				std::string _smoothingRequsted = option[ corrections::shapovalov_lib::shap_rama_smooth_level ];
				if (_smoothingRequsted.compare("1")==0 || _smoothingRequsted.compare("lowest_smooth")==0)
				{
					option[ corrections::shapovalov_lib::shap_rama_map ].value("scoring/score_functions/rama/shapovalov/kappa75/all.ramaProb");
					option[ corrections::shapovalov_lib::shap_rama_smooth_level ].value("1");
				}
				else if (_smoothingRequsted.compare("2")==0 || _smoothingRequsted.compare("lower_smooth")==0)
				{
					option[ corrections::shapovalov_lib::shap_rama_map ].value("scoring/score_functions/rama/shapovalov/kappa50/all.ramaProb");
					option[ corrections::shapovalov_lib::shap_rama_smooth_level ].value("2");
				}
				else if (_smoothingRequsted.compare("3")==0 || _smoothingRequsted.compare("higher_smooth")==0)
				{
					option[ corrections::shapovalov_lib::shap_rama_map ].value("scoring/score_functions/rama/shapovalov/kappa37.5/all.ramaProb");
					option[ corrections::shapovalov_lib::shap_rama_smooth_level ].value("3");
				}
				else if (_smoothingRequsted.compare("4")==0 || _smoothingRequsted.compare("highest_smooth")==0)
				{
					option[ corrections::shapovalov_lib::shap_rama_map ].value("scoring/score_functions/rama/shapovalov/kappa25/all.ramaProb");
					option[ corrections::shapovalov_lib::shap_rama_smooth_level ].value("4");
				}
				else
				{
					option[ corrections::shapovalov_lib::shap_rama_map ].value("scoring/score_functions/rama/shapovalov/kappa25/all.ramaProb");
					option[ corrections::shapovalov_lib::shap_rama_smooth_level ].value("4");
				}
			}
			else
			{
				option[ corrections::shapovalov_lib::shap_rama_smooth_level ].value("overriden");
			}
		}


		if (option[corrections::shapovalov_lib::shap_p_aa_pp_enable])
		{
			if ( ! option[ corrections::shapovalov_lib::shap_p_aa_pp ].user() )
			{
				std::string _smoothingRequsted = option[ corrections::shapovalov_lib::shap_p_aa_pp_smooth_level ];
				if (_smoothingRequsted.compare("1")==0 || _smoothingRequsted.compare("low_smooth")==0)
				{
					option[ corrections::shapovalov_lib::shap_p_aa_pp ].value("scoring/score_functions/P_AA_pp/shapovalov/10deg/kappa131/a20.prop");
					option[ corrections::shapovalov_lib::shap_p_aa_pp_smooth_level ].value("1");
				}
				else if (_smoothingRequsted.compare("2")==0 || _smoothingRequsted.compare("high_smooth")==0)
				{
					option[ corrections::shapovalov_lib::shap_p_aa_pp ].value("scoring/score_functions/P_AA_pp/shapovalov/10deg/kappa50/a20.prop");
					option[ corrections::shapovalov_lib::shap_p_aa_pp_smooth_level ].value("2");
				}
				else
				{
					option[ corrections::shapovalov_lib::shap_p_aa_pp ].value("scoring/score_functions/P_AA_pp/shapovalov/10deg/kappa50/a20.prop");
					option[ corrections::shapovalov_lib::shap_p_aa_pp_smooth_level ].value("2");
				}
			}
			else
			{
				option[ corrections::shapovalov_lib::shap_p_aa_pp_smooth_level ].value("overriden");
			}
		}


	}
}


void
init_crystal_refinement_correction() {
	//fpd  crystal-refinement specific changes
	if( option[cryst::crystal_refine]) {
		// use -correct rama fixes
		if ( ! option[ corrections::score::rama_not_squared ].user() ) {
			option[corrections::score::rama_not_squared].value( true );
		}
		if ( ! option[ corrections::score::rama_map ].user() ) {
			option[ corrections::score::rama_map ].value("scoring/score_functions/rama/Rama09_noEH_kernel25_it08.dat");
		}

		// read pdbs properly
		if ( ! option[ in::missing_density_to_jump ].user() ) {
			option[in::missing_density_to_jump].value( true );
		}
		if ( ! option[ in::preserve_crystinfo ].user() ) {
			option[in::preserve_crystinfo].value( true );
		}

		// hacky way of letting centroid structures be scored against crystal data
		if ( ! option[ out::file::no_output_cen ].user() ) {
			option[out::file::no_output_cen].value( true );
		}

		// internal coordinate minimization change: don't dampen RB movement
		if ( ! option[ optimization::scale_rb ].user() ) {
			option[optimization::scale_rb].value( 1.0 );
		}

	}
}


void
init_beta_correction() {
  if( option[corrections::beta_july15]() || option[corrections::beta_july15_cart]() ) {
		// atom clone
		if( ! option[ basic::options::OptionKeys::chemical::clone_atom_types ].user() ) {
			utility::vector1< std::string > params;
			params.push_back( "fa_standard:CH1:CH0" );
			params.push_back( "fa_standard:CH1:CH0R" );
			params.push_back( "fa_standard:Hpol:HS" );
			params.push_back( "fa_standard:NH2O:NHOQ" );
			params.push_back( "fa_standard:Ntrp:NtrH" );
			params.push_back( "fa_standard:Ntrp:NtrR" );
			params.push_back( "fa_standard:OH:OHT" );
			params.push_back( "fa_standard:OH:OHY" );
			params.push_back( "fa_standard:ONH2:ONHQ" );
			params.push_back( "fa_standard:OOC:OOCE" );
			params.push_back( "fa_standard:S:SH1" );
			option[ basic::options::OptionKeys::chemical::clone_atom_types ].value(params);
		} else {
			TR.Warning << "Flag -beta_july15 is set but -clone_atom_types are also specified.  Not changing atom properties!" << std::endl;
		}

		// atom reassign
		if( ! option[ basic::options::OptionKeys::chemical::reassign_atom_types ].user() ) {
			utility::vector1< std::string > params;
			params.push_back( "fa_standard:ARG:CZ:CH0R" );
			params.push_back( "fa_standard:ARG:NE:NtrR" );
			params.push_back( "fa_standard:CYS:HG:HS" );
			params.push_back( "fa_standard:CYS:SG:SH1" );
			params.push_back( "fa_standard:GLN:NE2:NHOQ" );
			params.push_back( "fa_standard:GLN:OE1:ONHQ" );
			params.push_back( "fa_standard:GLU:OE1:OOCE" );
			params.push_back( "fa_standard:GLU:OE2:OOCE" );
			params.push_back( "fa_standard:HIS_D:ND1:NtrH" );
			params.push_back( "fa_standard:HIS:CG:CH0" );
			params.push_back( "fa_standard:HIS:NE2:NtrH" );
			params.push_back( "fa_standard:MET:HG:HS" );
			params.push_back( "fa_standard:PHE:CG:CH0" );
			params.push_back( "fa_standard:THR:OG1:OHT" );
			params.push_back( "fa_standard:TRP:CD2:CH0" );
			params.push_back( "fa_standard:TRP:CE2:CH0" );
			params.push_back( "fa_standard:TRP:CG:CH0" );
			params.push_back( "fa_standard:TYR:CG:CH0" );
			params.push_back( "fa_standard:TYR:CZ:CH0" );
			params.push_back( "fa_standard:TYR:OH:OHY" );
			option[ basic::options::OptionKeys::chemical::reassign_atom_types ].value(params);
		} else {
			TR.Warning << "Flag -beta_july15 is set but -reassign_atom_types are also specified.  Not changing atom properties!" << std::endl;
		}

		// atom properties
		if( ! option[ basic::options::OptionKeys::chemical::set_atom_properties ].user() ) {
			utility::vector1< std::string > params;
			params.push_back( "fa_standard:aroC:LJ_RADIUS:1.99917604399492");
			params.push_back( "fa_standard:aroC:LJ_WDEPTH:0.119406055088945");
			params.push_back( "fa_standard:aroC:LK_DGFREE:1.33045639161877");
			params.push_back( "fa_standard:aroC:LK_VOLUME:16.704");
			params.push_back( "fa_standard:CAbb:LJ_RADIUS:2.02180593808171");
			params.push_back( "fa_standard:CAbb:LJ_WDEPTH:0.055337476050435");
			params.push_back( "fa_standard:CAbb:LK_DGFREE:1.67297469277228");
			params.push_back( "fa_standard:CAbb:LK_VOLUME:12.137");
			params.push_back( "fa_standard:CH0:LJ_RADIUS:2.02180593808171");
			params.push_back( "fa_standard:CH0:LJ_WDEPTH:0.055337476050435");
			params.push_back( "fa_standard:CH0:LK_DGFREE:0.278285165627887");
			params.push_back( "fa_standard:CH0:LK_VOLUME:8.9980");
			params.push_back( "fa_standard:CH0R:LJ_RADIUS:2.02180593808171");
			params.push_back( "fa_standard:CH0R:LJ_WDEPTH:0.055337476050435");
			params.push_back( "fa_standard:CH0R:LK_DGFREE:-0.00121092877425643");
			params.push_back( "fa_standard:CH0R:LK_VOLUME:8.9980");
			params.push_back( "fa_standard:CH1:LJ_RADIUS:2.02180593808171");
			params.push_back( "fa_standard:CH1:LJ_WDEPTH:0.055337476050435");
			params.push_back( "fa_standard:CH1:LK_DGFREE:-1.30308454629699");
			params.push_back( "fa_standard:CH1:LK_VOLUME:10.686");
			params.push_back( "fa_standard:CH2:LJ_RADIUS:2.07118221995967");
			params.push_back( "fa_standard:CH2:LJ_WDEPTH:0.143813358773253");
			params.push_back( "fa_standard:CH2:LK_DGFREE:-0.280286405370075");
			params.push_back( "fa_standard:CH2:LK_VOLUME:18.331");
			params.push_back( "fa_standard:CH3:LJ_RADIUS:2.02346742314654");
			params.push_back( "fa_standard:CH3:LJ_WDEPTH:0.240860458620729");
			params.push_back( "fa_standard:CH3:LK_DGFREE:6.11867159786138");
			params.push_back( "fa_standard:CH3:LK_VOLUME:25.855");
			params.push_back( "fa_standard:CNH2:LJ_RADIUS:2.01580593808171");
			params.push_back( "fa_standard:CNH2:LJ_WDEPTH:0.0802014845204291");
			params.push_back( "fa_standard:CNH2:LK_DGFREE:1.46174215885178");
			params.push_back( "fa_standard:CNH2:LK_VOLUME:13.500");
			params.push_back( "fa_standard:CObb:LJ_RADIUS:1.92513578221377");
			params.push_back( "fa_standard:CObb:LJ_WDEPTH:0.117295282631773");
			params.push_back( "fa_standard:CObb:LK_DGFREE:1.99927073996285");
			params.push_back( "fa_standard:CObb:LK_VOLUME:13.221");
			params.push_back( "fa_standard:COO:LJ_RADIUS:1.92513578221377");
			params.push_back( "fa_standard:COO:LJ_WDEPTH:0.117295282631773");
			params.push_back( "fa_standard:COO:LK_DGFREE:-3.34461092877425");
			params.push_back( "fa_standard:COO:LK_VOLUME:14.653");
			params.push_back( "fa_standard:Hapo:LJ_RADIUS:1.26522545740706");
			params.push_back( "fa_standard:Hapo:LJ_WDEPTH:0.0404580528958375");
			params.push_back( "fa_standard:Haro:LJ_RADIUS:1.26522545740706");
			params.push_back( "fa_standard:Haro:LJ_WDEPTH:0.0404580528958375");
			params.push_back( "fa_standard:HNbb:LJ_RADIUS:0.946805938081717");
			params.push_back( "fa_standard:HNbb:LJ_WDEPTH:0.0199253639540697");
			params.push_back( "fa_standard:Hpol:LJ_RADIUS:0.946805938081717");
			params.push_back( "fa_standard:Hpol:LJ_WDEPTH:0.0199253639540697");
			params.push_back( "fa_standard:HS:LJ_RADIUS:0.630805938081717");
			params.push_back( "fa_standard:HS:LJ_WDEPTH:0.0492014845204291");
			params.push_back( "fa_standard:Narg:LJ_RADIUS:1.75348059380817");
			params.push_back( "fa_standard:Narg:LJ_WDEPTH:0.170120148452043");
			params.push_back( "fa_standard:Narg:LK_DGFREE:-6.50568602515778");
			params.push_back( "fa_standard:Narg:LK_LAMBDA:3.500");
			params.push_back( "fa_standard:Narg:LK_VOLUME:15.717");
			params.push_back( "fa_standard:Nbb:LJ_RADIUS:1.75348059380817");
			params.push_back( "fa_standard:Nbb:LJ_WDEPTH:0.170120148452043");
			params.push_back( "fa_standard:Nbb:LK_DGFREE:-7.29136191571457");
			params.push_back( "fa_standard:Nbb:LK_VOLUME:15.992");
			params.push_back( "fa_standard:NH2O:LJ_RADIUS:1.75348059380817");
			params.push_back( "fa_standard:NH2O:LJ_WDEPTH:0.170120148452043");
			params.push_back( "fa_standard:NH2O:LK_DGFREE:-5.36711092877426");
			params.push_back( "fa_standard:NH2O:LK_VOLUME:15.689");
			params.push_back( "fa_standard:Nhis:LJ_RADIUS:1.75348059380817");
			params.push_back( "fa_standard:Nhis:LJ_WDEPTH:0.170120148452043");
			params.push_back( "fa_standard:Nhis:LK_DGFREE:-10.2857109287743");
			params.push_back( "fa_standard:Nhis:LK_VOLUME:9.3177");
			params.push_back( "fa_standard:NHOQ:LJ_RADIUS:1.75348059380817");
			params.push_back( "fa_standard:NHOQ:LJ_WDEPTH:0.170120148452043");
			params.push_back( "fa_standard:NHOQ:LK_DGFREE:-5.20281092877426");
			params.push_back( "fa_standard:Nlys:LJ_RADIUS:1.75348059380817");
			params.push_back( "fa_standard:Nlys:LJ_WDEPTH:0.170120148452043");
			params.push_back( "fa_standard:Nlys:LK_DGFREE:-16.2981109287742");
			params.push_back( "fa_standard:Nlys:LK_LAMBDA:3.500");
			params.push_back( "fa_standard:Nlys:LK_VOLUME:16.514");
			params.push_back( "fa_standard:Npro:LJ_RADIUS:1.75348059380817");
			params.push_back( "fa_standard:Npro:LJ_WDEPTH:0.170120148452043");
			params.push_back( "fa_standard:Npro:LK_DGFREE:-3.71951092877425");
			params.push_back( "fa_standard:Npro:LK_VOLUME:3.7181");
			params.push_back( "fa_standard:NtrH:LJ_RADIUS:1.75348059380817");
			params.push_back( "fa_standard:NtrH:LJ_WDEPTH:0.170120148452043");
			params.push_back( "fa_standard:NtrH:LK_DGFREE:-6.35271092877426");
			params.push_back( "fa_standard:NtrH:LK_VOLUME:9.2829");
			params.push_back( "fa_standard:Ntrp:LJ_RADIUS:1.75348059380817");
			params.push_back( "fa_standard:Ntrp:LJ_WDEPTH:0.170120148452043");
			params.push_back( "fa_standard:Ntrp:LK_DGFREE:-5.95091092877426");
			params.push_back( "fa_standard:Ntrp:LK_VOLUME:9.5221");
			params.push_back( "fa_standard:NtrR:LJ_RADIUS:1.75348059380817");
			params.push_back( "fa_standard:NtrR:LJ_WDEPTH:0.170120148452043");
			params.push_back( "fa_standard:NtrR:LK_DGFREE:-9.77811092877423");
			params.push_back( "fa_standard:NtrR:LK_VOLUME:9.7792");
			params.push_back( "fa_standard:OCbb:LJ_RADIUS:1.55048059380817");
			params.push_back( "fa_standard:OCbb:LJ_WDEPTH:0.161120148452043");
			params.push_back( "fa_standard:OCbb:LK_DGFREE:-5.20513751964092");
			params.push_back( "fa_standard:OCbb:LK_VOLUME:12.196");
			params.push_back( "fa_standard:OH:LJ_RADIUS:1.55048059380817");
			params.push_back( "fa_standard:OH:LJ_WDEPTH:0.160120148452043");
			params.push_back( "fa_standard:OH:LK_DGFREE:-6.38881994292366");
			params.push_back( "fa_standard:OH:LK_VOLUME:10.722");
			params.push_back( "fa_standard:OHT:LJ_RADIUS:1.55048059380817");
			params.push_back( "fa_standard:OHT:LJ_WDEPTH:0.160120148452043");
			params.push_back( "fa_standard:OHT:LK_DGFREE:-8.67363943218534");
			params.push_back( "fa_standard:OHT:LK_VOLUME:10.843");
			params.push_back( "fa_standard:OHY:LJ_RADIUS:1.55048059380817");
			params.push_back( "fa_standard:OHY:LJ_WDEPTH:0.160120148452043");
			params.push_back( "fa_standard:OHY:LK_DGFREE:-4.16723176579274");
			params.push_back( "fa_standard:OHY:LK_VOLUME:10.628");
			params.push_back( "fa_standard:ONH2:LJ_RADIUS:1.55048059380817");
			params.push_back( "fa_standard:ONH2:LJ_WDEPTH:0.161120148452043");
			params.push_back( "fa_standard:ONH2:LK_DGFREE:-5.87105842194542");
			params.push_back( "fa_standard:ONH2:LK_VOLUME:10.102");
			params.push_back( "fa_standard:ONHQ:LJ_RADIUS:1.55048059380817");
			params.push_back( "fa_standard:ONHQ:LJ_WDEPTH:0.161120148452043");
			params.push_back( "fa_standard:ONHQ:LK_DGFREE:-7.09171092877425");
			params.push_back( "fa_standard:OOC:LJ_RADIUS:1.49381727066464");
			params.push_back( "fa_standard:OOC:LJ_WDEPTH:0.213201484520429");
			params.push_back( "fa_standard:OOC:LK_DGFREE:-7.86291092877426");
			params.push_back( "fa_standard:OOC:LK_LAMBDA:3.500");
			params.push_back( "fa_standard:OOC:LK_VOLUME:9.9956");
			params.push_back( "fa_standard:OOCE:LJ_RADIUS:1.49381727066464");
			params.push_back( "fa_standard:OOCE:LJ_WDEPTH:0.213201484520429");
			params.push_back( "fa_standard:OOCE:LK_DGFREE:-7.86281092877427");
			params.push_back( "fa_standard:OOCE:LK_LAMBDA:3.500");
			params.push_back( "fa_standard:OOCE:LK_VOLUME:9.9956");
			params.push_back( "fa_standard:S:LJ_RADIUS:1.88248061423444");
			params.push_back( "fa_standard:S:LJ_WDEPTH:0.433201484520429");
			params.push_back( "fa_standard:S:LK_DGFREE:-0.791901963657784");
			params.push_back( "fa_standard:S:LK_VOLUME:17.640");
			params.push_back( "fa_standard:SH1:LJ_RADIUS:1.88248061423444");
			params.push_back( "fa_standard:SH1:LJ_WDEPTH:0.433201484520429");
			params.push_back( "fa_standard:SH1:LK_DGFREE:3.84568907122574");
			params.push_back( "fa_standard:SH1:LK_VOLUME:23.240");
			option[ basic::options::OptionKeys::chemical::set_atom_properties ].value(params);
		} else {
			TR.Warning << "Flag -beta_july15 is set but -set_atom_properties are also specified.  Not changing atom properties!" << std::endl;
		}

		// hbonds
		if( ! option[ basic::options::OptionKeys::score::hb_acc_strength ].user() ) {
			utility::vector1< std::string > params;
			params.push_back( "hbacc_AHX:0.94");
			params.push_back( "hbacc_CXA:1.08");
			params.push_back( "hbacc_CXL:1.00");
			params.push_back( "hbacc_HXL:1.00");
			params.push_back( "hbacc_IMD:1.05");
			params.push_back( "hbacc_IME:0.94");
			params.push_back( "hbacc_PBA:0.94");
			option[ basic::options::OptionKeys::score::hb_acc_strength ].value(params);
		} else {
			TR.Warning << "Flag -beta_july15 is set but -hb_acc_strength are also specified.  Not changing atom properties!" << std::endl;
		}
		if( ! option[ basic::options::OptionKeys::score::hb_don_strength ].user() ) {
			utility::vector1< std::string > params;
			params.push_back( "hbdon_AHX:0.92");
			params.push_back( "hbdon_AMO:1.01");
			params.push_back( "hbdon_CXA:1.08");
			params.push_back( "hbdon_GDE:0.95");
			params.push_back( "hbdon_GDH:0.96");
			params.push_back( "hbdon_HXL:0.85");
			params.push_back( "hbdon_IMD:1.08");
			params.push_back( "hbdon_IME:1.31");
			params.push_back( "hbdon_IND:1.14");
			params.push_back( "hbdon_PBA:1.26");
			option[ basic::options::OptionKeys::score::hb_don_strength ].value(params);
		} else {
			TR.Warning << "Flag -beta_july15 is set but -hb_don_strength are also specified.  Not changing atom properties!" << std::endl;
		}

		// unmodifypot
		if( ! option[ basic::options::OptionKeys::score::unmodifypot ].user() ) {
			option[ basic::options::OptionKeys::score::unmodifypot ].value(true);
		}

		// sigmoidal dielectric
		if( ! option[ basic::options::OptionKeys::score::elec_sigmoidal_die ].user() ) {
			option[ basic::options::OptionKeys::score::elec_sigmoidal_die ].value(true);
		}
		if( ! option[ basic::options::OptionKeys::score::elec_sigmoidal_die_D ].user() ) {
			option[ basic::options::OptionKeys::score::elec_sigmoidal_die_D ].value(137.0);
		}
		if( ! option[ basic::options::OptionKeys::score::elec_sigmoidal_die_D0 ].user() ) {
			option[ basic::options::OptionKeys::score::elec_sigmoidal_die_D0 ].value(9.75);
		}
		if( ! option[ basic::options::OptionKeys::score::elec_sigmoidal_die_S ].user() ) {
			option[ basic::options::OptionKeys::score::elec_sigmoidal_die_S ].value(0.44);
		}

		// weights file
    if ( ! option[ score::weights ].user() ) {
			if( option[corrections::beta_july15_cart]() || option[optimization::nonideal]() ) {
				option[ score::weights ].value( "beta_july15_cart.wts" );
			} else {
				option[ score::weights ].value( "beta_july15.wts" );
			}
    } else {
			TR.Warning << "Flag -beta_july15 is set but -weights are also specified.  Not changing input weights file!" << std::endl;
		}

	}
}


void
init_nonideal_correction() {
	//PTC set up for flexible bond geometries
	if( option[optimization::nonideal]) {
		// add jumps to missing densities so cart_bonded_length behaves correctly
		if ( ! option[ in::missing_density_to_jump ].user() ) {
			option[ in::missing_density_to_jump ].value( true );
		}

		// use dualspace relax protocol - start with internal coordinate, finish with Cartesian
		if ( ! option[ relax::dualspace ].user() ) {
			option[ relax::dualspace ].value( true );
		}

		// use flexible bond angles during internal coordinate stage of dualspace relax
		if ( ! option[ relax::minimize_bond_angles ].user() ) {
			option[ relax::minimize_bond_angles ].value( true );
		}

		// limit default_max_cycles
		if ( ! option[ optimization::default_max_cycles ].user() ) {
			option[ optimization::default_max_cycles ].value( 200 );
		}

		// enable for rtmin/minpack - use nonideal unless nonideal option is defined or cartesian is specified as true
		if ( (! option[ optimization::scmin_nonideal ].user()) || option[ optimization::scmin_cartesian ] ) {
			option[ optimization::scmin_nonideal ].value( true );
		}

		// talaris2013_cart.wts - enables cart_bonded, disables pro_close, new reference weights
		if ( ! option[ score::weights ].user() ) {
			option[ score::weights ].value( "talaris2013_cart.wts" );
		}
		else {
			scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
			if ( !sfxn->ready_for_nonideal_scoring() )
				TR.Warning << "option[ optimization::nonideal ] - scorefunction not set up for nonideal/Cartesian scoring (safe to ignore this warning if setting weights via RosettaScripts)" << std::endl;
		}
	}
}

//void
//revert_hbond_sp2_correction_to_score12() {
//	if( ! option[ corrections::score::hb_sp2_chipen ].user() ) {
//		option[ corrections::score::hb_sp2_chipen ].value( false );
//	}
//
//	if( ! option[ corrections::score::hb_sp2_outer_width ].user() ) {
//		option[ corrections::score::hb_sp2_outer_width ].value( 0.33333 );
//	}
//
//	if( ! option[ corrections::score::hb_fade_energy ].user() ) {
//		option[ corrections::score::hb_fade_energy ].value( false );
//	}
//
//	if( ! option[ corrections::score::hbond_measure_sp3acc_BAH_from_hvy ].user() ) {
//		option[ corrections::score::hbond_measure_sp3acc_BAH_from_hvy ].value( false );
//	}
//
//	if ( ! option[ corrections::chemical::expand_st_chi2sampling ].user() ) {
//		option[ corrections::chemical::expand_st_chi2sampling ].value( false );
//	}
//
//	if( ! option[ score::hbond_params ].user() ) {
//		option[ score::hbond_params ].value( "score12_params" );
//	}
//
//	if( ! option[ score::smooth_fa_elec ].user() ) {
//		option[ score::smooth_fa_elec ].value( false );
//	}
//
//	if( ! option[ score::elec_min_dis ].user() ) {
//		option[ score::elec_min_dis ].value( 1.5 );
//	}
//
//}


void
init_restore_score12prime() {
	if (option[ corrections::score::score12prime ] ) {
		revert_to_pre_talaris_2013_defaults();
		option[ score::weights ].value( "score12prime.wts" );
		option[ score::patch ].value( "" );
	}
}

//MaximCode:
void
init_score_function_corrections(){
	init_revert_to_pre_talaris_2013_mistake();
	init_hbond_sp2_correction();
	init_facts_correction();
	init_correct_correction();
	init_crystal_refinement_correction();
	init_beta_correction();
	init_nonideal_correction(); // keep after beta
	init_dna_correction();

	init_shapovalov_lib_fixes_enable_correction();
}

void
check_score_function_sanity(
	std::string const & scorefxn_key,
	bool only_warn
) {
	if( scorefxn_key == "score12" ){
		if( ! option[ mistakes::restore_pre_talaris_2013_behavior ] ) {
			TR.Error
				<< "**********************************************" << std::endl
				<< "To use the '" << scorefxn_key << "' score function" << std::endl
				<< "please using the " << std::endl
				<< std::endl
				<< "        -restore_pre_talaris_2013_behavior" << std::endl
				<< std::endl
				<< "flag on the command line." << std::endl
				<< "**********************************************" << std::endl;
			if(!only_warn){
				throw utility::excn::EXCN_BadInput("Missing -restore_pre_talaris_2013_behavior flag.");
			}
		}
	} else if( scorefxn_key == "talaris2013" ){
		if( option[ mistakes::restore_pre_talaris_2013_behavior ] ) {
			TR.Error
				<< "**********************************************" << std::endl
				<< "To use the '" << scorefxn_key << "' score function" << std::endl
				<< "please don't using the " << std::endl
				<< std::endl
				<< "        -restore_pre_talaris_2013_behavior" << std::endl
				<< std::endl
				<< "flag on the command line." << std::endl
				<< "**********************************************" << std::endl;
			if(!only_warn){
				throw utility::excn::EXCN_BadInput("Using -restore_pre_talaris_2013_behavior flag with talaris2013 energy function");
			}
		}
	}

}

void
init_dna_correction()
{

	utility::vector1< std::string > dna_residue_types;
	dna_residue_types.push_back("ADE");
	dna_residue_types.push_back("CYT");
	dna_residue_types.push_back("GUA");
	dna_residue_types.push_back("THY");

	if ( option[ corrections::newdna ] ) {

		{ // cloning atomtypes
			utility::vector1< std::string > mods;

			/// reassign the N9/N1 atoms from Ntrp atomtype to Npro atomtype
			// mods.push_back( "fa_standard:OOC:OPdd" ); // phosphate oxygens
			// mods.push_back( "fa_standard:OH:OEdd" ); // phosphate ester oxygens, was ONH2, doesnt make much difference
			// mods.push_back( "fa_standard:OH:ORdd" ); // ribose oxygen
			mods.push_back( "fa_standard:Phos:Pdna" ); // Phos in phosphate backbone
			mods.push_back( "fa_standard:OCbb:OBdd" ); // carbonyl oxygens in the bases

			if ( basic::options::option[ basic::options::OptionKeys::chemical::clone_atom_types ].user() ) {
				/// add these at the end so they take precedence
				utility::vector1< std::string > const user_mods
					( basic::options::option[ basic::options::OptionKeys::chemical::clone_atom_types ]() );
				for ( Size i=1; i<= user_mods.size(); ++i ) {
					if ( std::find( mods.begin(), mods.end(), user_mods[i] ) == mods.end() ) {
						mods.push_back( user_mods[i] );
					}
				}
			}
			basic::options::option[ basic::options::OptionKeys::chemical::clone_atom_types ].value( mods );

		}

		{ // reassigning atomtype properties
			utility::vector1< std::string > mods;
			mods.push_back( "fa_standard:Pdna:LK_DGFREE:-4.1" ); // from way back, Jim H., based on S
			mods.push_back( "fa_standard:Pdna:LK_VOLUME:14.7" ); // -- ditto --

			if ( basic::options::option[ basic::options::OptionKeys::chemical::set_atom_properties ].user() ) {
				/// add these at the end so they take precedence
				utility::vector1< std::string > const user_mods
					( basic::options::option[ basic::options::OptionKeys::chemical::set_atom_properties ]() );
				for ( Size i=1; i<= user_mods.size(); ++i ) {
					if ( std::find( mods.begin(), mods.end(), user_mods[i] ) == mods.end() ) {
						mods.push_back( user_mods[i] );
					}
				}
			}
			basic::options::option[ basic::options::OptionKeys::chemical::set_atom_properties ].value( mods );
		}


		{ // reassigning atomtypes

			utility::vector1< std::string > mods;

			/// reassign the N9/N1 atoms from Ntrp atomtype to Npro atomtype
			mods.push_back( "fa_standard:ADE:N9:Npro" );
			mods.push_back( "fa_standard:GUA:N9:Npro" );
			mods.push_back( "fa_standard:THY:N1:Npro" );
			mods.push_back( "fa_standard:CYT:N1:Npro" );

			for ( Size ii=1; ii<= dna_residue_types.size(); ++ii ) {
				std::string const & resname( dna_residue_types[ii] );
				/// phosphate Phos
				mods.push_back( "fa_standard:"+resname+":P:Pdna" );
			}

			/// the OCbb atomtypes in the bases; Etable.cc:modify_pot_one_pair does funny things to OCbb atomtype-- trouble?
			mods.push_back( "fa_standard:GUA:O6:OBdd" );
			mods.push_back( "fa_standard:CYT:O2:OBdd" );
			mods.push_back( "fa_standard:THY:O2:OBdd" );
			mods.push_back( "fa_standard:THY:O4:OBdd" );

			if ( basic::options::option[ basic::options::OptionKeys::chemical::reassign_atom_types ].user() ) {
				/// add these at the end so they take precedence
				utility::vector1< std::string > const user_mods
					( basic::options::option[ basic::options::OptionKeys::chemical::reassign_atom_types ]() );
				for ( Size i=1; i<= user_mods.size(); ++i ) {
					if ( std::find( mods.begin(), mods.end(), user_mods[i] ) == mods.end() ) {
						mods.push_back( user_mods[i] );
					}
				}
			}
			basic::options::option[ basic::options::OptionKeys::chemical::reassign_atom_types ].value( mods );
		}

		{ // reassigning icoor values

			utility::vector1< std::string > mods;

			for( Size ii=1; ii<= dna_residue_types.size(); ++ii ) {
				std::string const & resname( dna_residue_types[ii]);
				/// the default way is to have O4 as the torsion atom for H2'' but that's not as good as this:
				/// not sure why H2' follows H2'' but this is how it is in the params files right now...
				mods.push_back( "fa_standard:"+resname+":H2'':120.0,68.0,1.1,C2',C1',C3'" );
				mods.push_back( "fa_standard:"+resname+ ":H2':120.0,68.0,1.1,C2',C1',H2''" );

				/// it might be useful to standardize the backbone geometry across residue types,
				/// for when we are making mutations
				/// these are some of the changes from the blab branch
				mods.push_back( "fa_standard:"+resname+":UPPER:-180.0,60.2,1.608,O3',C3',C4'");
				mods.push_back( "fa_standard:"+resname+":LOWER:-64.0,76.3,1.608,P,O5',C5'");
			}


			if ( basic::options::option[ basic::options::OptionKeys::chemical::reassign_icoor ].user() ) {
				/// add these at the end so they take precedence
				utility::vector1< std::string > const user_mods
					( basic::options::option[ basic::options::OptionKeys::chemical::reassign_icoor ]() );
				for ( Size i=1; i<= user_mods.size(); ++i ) {
					if ( std::find( mods.begin(), mods.end(), user_mods[i] ) == mods.end() ) {
						mods.push_back( user_mods[i] );
					}
				}
			}
			basic::options::option[ basic::options::OptionKeys::chemical::reassign_icoor ].value( mods );
		}


	}


}
} // namespace
} // namespace
