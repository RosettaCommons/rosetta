// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/init/score_function_corrections.cc
/// @brief  Initialize Score Funciton Corrections
/// @author Matthew O'Meara
/// @author Frank and/or Hahnbeom did the beta 15 stuff

// Unit headers
#include <core/init/score_function_corrections.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/cryst.OptionKeys.gen.hh>
#include <basic/options/keys/ddg.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/mistakes.OptionKeys.gen.hh>
#include <basic/options/keys/qsar.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/optimization.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>
#include <basic/options/keys/rings.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <utility/options/OptionCollection.hh>

// Core headers
#include <core/types.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// C++ headers
#include <istream>
#include <string>

using namespace basic::options;
using namespace basic::options::OptionKeys;

//MaximCode:
static basic::Tracer TR( "core.init.score_function_corrections" );

namespace core {
namespace init {

struct pre_talaris_2013_behavior_settings {
	pre_talaris_2013_behavior_settings();

	bool pre_talaris2013_geometries{ true };
	bool hb_sp2_chipen{ false };
	bool hb_fade_energy{ false };
	bool hbond_measure_sp3acc_BAH_from_hvy{ false };
	Real lj_hbond_hdis = 1.95;
	Real lj_hbond_OH_donor_dis = 3.0;
	Real hb_sp2_outer_width = 0.33333;
	Real expand_st_chi2sampling{ false };
	std::string score_weights;
	std::string score_patch;
	bool analytic_etable_evaluation{ false };
	std::string hbond_params;
	bool smooth_fa_elec{ false };
	Real elec_min_dis = 1.5;
	bool dun10{ false };
	bool use_bicubic_interpolation{ false };
};

pre_talaris_2013_behavior_settings::pre_talaris_2013_behavior_settings() :
	score_weights( "pre_talaris_2013_standard.wts" ),
	score_patch( "" ), // the default patch to score12 is handled in the core::scoring::get_score_function
	hbond_params( "score12_params" )
{}

pre_talaris_2013_behavior_settings const restore_sc12_settings;

static basic::Tracer TR( "core.init.score_function_corrections" );

// entry point: apply all flag-based corrections
void
init_score_function_corrections( utility::options::OptionCollection & options ){
	init_pre_beta_nov15_etable( options ); // call first

	init_revert_to_talaris( options ); //option[ corrections::restore_talaris_behavior ] keep before pre_talaris_2013
	init_revert_to_pre_talaris_2013_mistake( options );
	init_hbond_sp2_correction( options );
	init_facts_correction( options );
	init_correct_correction( options );
	init_crystal_refinement_correction( options );
	init_beta_correction( options );
	init_nonideal_correction( options ); // keep after beta
	init_shapovalov_lib_fixes_enable_correction( options ); // keep after beta
	init_dna_correction( options );
	init_spades_score_function_correction( options );
}

void
init_revert_to_pre_talaris_2013_mistake( utility::options::OptionCollection & options ) {
	if ( options[ mistakes::restore_pre_talaris_2013_behavior ] ) {
		revert_to_pre_talaris_2013_defaults( options );
	}
}

void revert_to_pre_talaris_2013_defaults( utility::options::OptionCollection & options ) {
	if ( ! options[ mistakes::chemical::pre_talaris2013_geometries ].user() ) {
		options[ mistakes::chemical::pre_talaris2013_geometries ].default_value( restore_sc12_settings.pre_talaris2013_geometries );
	}
	if ( ! options[ corrections::score::hb_sp2_chipen ].user() ) {
		options[ corrections::score::hb_sp2_chipen ].default_value( restore_sc12_settings.hb_sp2_chipen );
	}
	if ( ! options[ corrections::score::hb_fade_energy ].user() ) {
		options[ corrections::score::hb_fade_energy ].default_value( restore_sc12_settings.hb_fade_energy );
	}
	if ( ! options[ corrections::score::hbond_measure_sp3acc_BAH_from_hvy ].user() ) {
		options[ corrections::score::hbond_measure_sp3acc_BAH_from_hvy ].default_value( restore_sc12_settings.hbond_measure_sp3acc_BAH_from_hvy );
	}
	if ( ! options[ corrections::score::lj_hbond_hdis ].user() ) {
		options[ corrections::score::lj_hbond_hdis ].default_value( restore_sc12_settings.lj_hbond_hdis );
	}
	if ( ! options[ corrections::score::lj_hbond_OH_donor_dis ].user() ) {
		options[ corrections::score::lj_hbond_OH_donor_dis ].default_value( restore_sc12_settings.lj_hbond_OH_donor_dis );
	}
	if ( ! options[ corrections::score::hb_sp2_outer_width ].user() ) {
		options[ corrections::score::hb_sp2_outer_width ].default_value( restore_sc12_settings.hb_sp2_outer_width );
	}
	if ( ! options[ corrections::chemical::expand_st_chi2sampling ].user() ) {
		options[ corrections::chemical::expand_st_chi2sampling ].default_value( restore_sc12_settings.expand_st_chi2sampling );
	}
	if ( ! options[ score::weights ].user() ) {
		options[ score::weights ].default_value( restore_sc12_settings.score_weights );
		/// Don't set the score12 weights patch if the user has provided a value for either score::weights or score::patch
		if ( ! options[ score::patch ].user() ) {
			utility::vector1< std::string > default_score_patch(1);
			default_score_patch[1] = restore_sc12_settings.score_patch;
			options[ score::patch ].default_value( default_score_patch );
		}
	}
	if ( ! options[ score::analytic_etable_evaluation ].user() ) {
		options[ score::analytic_etable_evaluation ].default_value( restore_sc12_settings.analytic_etable_evaluation );
	}
	if ( ! options[ score::hbond_params ].user() ) {
		options[ score::hbond_params ].default_value( restore_sc12_settings.hbond_params );
	}
	if ( ! options[ score::smooth_fa_elec ].user() ) {
		options[ score::smooth_fa_elec ].default_value( restore_sc12_settings.smooth_fa_elec );
	}
	if ( ! options[ score::elec_min_dis ].user() ) {
		options[ score::elec_min_dis ].default_value( restore_sc12_settings.elec_min_dis );
	}
	if ( ! options[ basic::options::OptionKeys::chemical::set_atom_properties ].user() ) {
		utility::vector1< std::string > params;
		params.push_back("fa_standard:ONH2:LK_DGFREE:-10.0000");
		params.push_back("fa_standard:NH2O:LK_DGFREE:-10.0000");
		params.push_back("fa_standard:Narg:LK_DGFREE:-11.0000");
		params.push_back("fa_standard:OH:LK_DGFREE:-6.7700" );
		params.push_back("fa_standard:OW:LK_DGFREE:-6.7700" );
		options[ basic::options::OptionKeys::chemical::set_atom_properties ].value(params);
	}

	if ( ! options[ corrections::score::dun10 ].user() ) {
		options[corrections::score::dun10].default_value( restore_sc12_settings.dun10 );
	}

	if ( ! options[ corrections::score::use_bicubic_interpolation ].user() ) {
		options[corrections::score::use_bicubic_interpolation].default_value( restore_sc12_settings.use_bicubic_interpolation );
	}
	if ( ! options[ qsar::weights ].user() ) {
		options[ qsar::weights ].default_value( restore_sc12_settings.score_weights );
	}
	if ( ! options[ ddg::min_cst_weights ].user() ) {
		options[ ddg::min_cst_weights ].default_value( restore_sc12_settings.score_weights );
	}
	if ( ! options[ score::pack_weights ].user() ) {
		options[ score::pack_weights ].default_value( restore_sc12_settings.score_weights );
	}
}

void
init_pre_beta_nov15_etable( utility::options::OptionCollection & options ){
	// atom properties
	// should work on setting atom_type_set...  make sure not to call any beta...
	if ( options[ corrections::restore_talaris_behavior ] ||
			options[ mistakes::restore_pre_talaris_2013_behavior ] ||
			options[ corrections::hbond_sp2_correction ] ||
			options[corrections::facts_default] ) {
		if ( options[corrections::beta]() ||
				options[corrections::beta_jan25]() || options[corrections::beta_jan25_cart]() ||
				options[corrections::beta_nov16]() || options[corrections::beta_nov16_cart]() ||
				options[corrections::beta_july15]() || options[corrections::beta_july15_cart]() ) {
			utility_exit_with_message("pre_beta_nov15_etable is not compatible with any of -beta flags!");
		}
	} else {
		return;
	}

	if ( ! options[ basic::options::OptionKeys::chemical::set_atom_properties ].user() ) {
		utility::vector1< std::string > params;
		params.push_back("fa_standard:CNH2:LJ_RADIUS:2.0000");
		params.push_back("fa_standard:CNH2:LJ_WDEPTH:0.1200");
		params.push_back("fa_standard:CNH2:LK_DGFREE:0.0000");
		params.push_back("fa_standard:CNH2:LK_VOLUME:14.7000");
		params.push_back("fa_standard:COO:LJ_RADIUS:2.0000");
		params.push_back("fa_standard:COO:LJ_WDEPTH:0.1200");
		params.push_back("fa_standard:COO:LK_DGFREE:-1.4000");
		params.push_back("fa_standard:COO:LK_VOLUME:8.3000");
		params.push_back("fa_standard:CH1:LJ_RADIUS:2.0000");
		params.push_back("fa_standard:CH1:LJ_WDEPTH:0.0486");
		params.push_back("fa_standard:CH1:LK_DGFREE:-0.2500");
		params.push_back("fa_standard:CH1:LK_VOLUME:23.7000");
		params.push_back("fa_standard:CH2:LJ_RADIUS:2.0000");
		params.push_back("fa_standard:CH2:LJ_WDEPTH:0.1142");
		params.push_back("fa_standard:CH2:LK_DGFREE:0.5200");
		params.push_back("fa_standard:CH2:LK_VOLUME:22.4000");
		params.push_back("fa_standard:CH3:LJ_RADIUS:2.0000");
		params.push_back("fa_standard:CH3:LJ_WDEPTH:0.1811");
		params.push_back("fa_standard:CH3:LK_DGFREE:1.5000");
		params.push_back("fa_standard:CH3:LK_VOLUME:30.0000");
		params.push_back("fa_standard:aroC:LJ_RADIUS:2.0000");
		params.push_back("fa_standard:aroC:LJ_WDEPTH:0.1200");
		params.push_back("fa_standard:aroC:LK_DGFREE:0.0800");
		params.push_back("fa_standard:aroC:LK_VOLUME:18.4000");
		params.push_back("fa_standard:Ntrp:LJ_RADIUS:1.7500");
		params.push_back("fa_standard:Ntrp:LJ_WDEPTH:0.2384");
		params.push_back("fa_standard:Ntrp:LK_DGFREE:-8.9000");
		params.push_back("fa_standard:Ntrp:LK_VOLUME:4.4000");
		params.push_back("fa_standard:Nhis:LJ_RADIUS:1.7500");
		params.push_back("fa_standard:Nhis:LJ_WDEPTH:0.2384");
		params.push_back("fa_standard:Nhis:LK_DGFREE:-4.0000");
		params.push_back("fa_standard:Nhis:LK_VOLUME:4.4000");
		params.push_back("fa_standard:NH2O:LJ_RADIUS:1.7500");
		params.push_back("fa_standard:NH2O:LJ_WDEPTH:0.2384");
		//params.push_back("fa_standard:NH2O:LK_DGFREE:-7.8000");
		params.push_back("fa_standard:NH2O:LK_VOLUME:11.2000");
		params.push_back("fa_standard:Nlys:LJ_RADIUS:1.7500");
		params.push_back("fa_standard:Nlys:LJ_WDEPTH:0.2384");
		params.push_back("fa_standard:Nlys:LK_DGFREE:-20.0000");
		params.push_back("fa_standard:Nlys:LK_VOLUME:11.2000");
		params.push_back("fa_standard:Narg:LJ_RADIUS:1.7500");
		params.push_back("fa_standard:Narg:LJ_WDEPTH:0.2384");
		//params.push_back("fa_standard:Narg:LK_DGFREE:-10.0000");
		params.push_back("fa_standard:Narg:LK_VOLUME:11.2000");
		params.push_back("fa_standard:Npro:LJ_RADIUS:1.7500");
		params.push_back("fa_standard:Npro:LJ_WDEPTH:0.2384");
		params.push_back("fa_standard:Npro:LK_DGFREE:-1.5500");
		params.push_back("fa_standard:Npro:LK_VOLUME:0.0000");
		params.push_back("fa_standard:OH:LJ_RADIUS:1.5500");
		params.push_back("fa_standard:OH:LJ_WDEPTH:0.1591");
		//params.push_back("fa_standard:OH:LK_DGFREE:-6.7000");
		params.push_back("fa_standard:OH:LK_VOLUME:10.8000");
		params.push_back("fa_standard:OW:LJ_RADIUS:1.5500");
		params.push_back("fa_standard:OW:LJ_WDEPTH:0.1591");
		//params.push_back("fa_standard:OW:LK_DGFREE:-6.7000");
		params.push_back("fa_standard:OW:LK_VOLUME:10.8000");
		params.push_back("fa_standard:ONH2:LJ_RADIUS:1.5500");
		params.push_back("fa_standard:ONH2:LJ_WDEPTH:0.1591");
		//params.push_back("fa_standard:ONH2:LK_DGFREE:-5.8500");
		params.push_back("fa_standard:ONH2:LK_VOLUME:10.8000");
		params.push_back("fa_standard:OOC:LJ_RADIUS:1.5500");
		params.push_back("fa_standard:OOC:LJ_WDEPTH:0.2100");
		params.push_back("fa_standard:OOC:LK_DGFREE:-10.0000");
		params.push_back("fa_standard:OOC:LK_VOLUME:10.8000");
		params.push_back("fa_standard:S:LJ_RADIUS:1.9000");
		params.push_back("fa_standard:S:LJ_WDEPTH:0.1600");
		params.push_back("fa_standard:S:LK_DGFREE:-4.1000");
		params.push_back("fa_standard:S:LK_VOLUME:14.7000");
		params.push_back("fa_standard:Nbb:LJ_RADIUS:1.7500");
		params.push_back("fa_standard:Nbb:LJ_WDEPTH:0.2384");
		params.push_back("fa_standard:Nbb:LK_DGFREE:-5.0000");
		params.push_back("fa_standard:Nbb:LK_VOLUME:4.4000");
		params.push_back("fa_standard:CAbb:LJ_RADIUS:2.0000");
		params.push_back("fa_standard:CAbb:LJ_WDEPTH:0.0486");
		params.push_back("fa_standard:CAbb:LK_DGFREE:1.0000");
		params.push_back("fa_standard:CAbb:LK_VOLUME:23.7000");
		params.push_back("fa_standard:CObb:LJ_RADIUS:2.0000");
		params.push_back("fa_standard:CObb:LJ_WDEPTH:0.1400");
		params.push_back("fa_standard:CObb:LK_DGFREE:1.0000");
		params.push_back("fa_standard:CObb:LK_VOLUME:14.7000");
		params.push_back("fa_standard:OCbb:LJ_RADIUS:1.5500");
		params.push_back("fa_standard:OCbb:LJ_WDEPTH:0.1591");
		params.push_back("fa_standard:OCbb:LK_DGFREE:-5.0000");
		params.push_back("fa_standard:OCbb:LK_VOLUME:10.8000");
		params.push_back("fa_standard:Hpol:LJ_RADIUS:1.0000");
		params.push_back("fa_standard:Hpol:LJ_WDEPTH:0.0500");
		params.push_back("fa_standard:Hapo:LJ_RADIUS:1.2000");
		params.push_back("fa_standard:Hapo:LJ_WDEPTH:0.0500");
		params.push_back("fa_standard:Haro:LJ_RADIUS:1.2000");
		params.push_back("fa_standard:Haro:LJ_WDEPTH:0.0500");
		params.push_back("fa_standard:HNbb:LJ_RADIUS:1.0000");
		params.push_back("fa_standard:HNbb:LJ_WDEPTH:0.0500");

		// lambda
		params.push_back("fa_standard:Narg:LK_LAMBDA:6.0000");
		params.push_back("fa_standard:Nlys:LK_LAMBDA:6.0000");
		params.push_back("fa_standard:OOC:LK_LAMBDA:6.0000");

		// this is very hacky... put it here because splitting it into revert_to_pre_talaris_2013_default
		//will clear all the information above
		// Doing this because we know we will get rid of this .cc file in near future
		if ( options[ mistakes::restore_pre_talaris_2013_behavior ] ) { //pre-talaris
			params.push_back("fa_standard:ONH2:LK_DGFREE:-10.0000");
			params.push_back("fa_standard:NH2O:LK_DGFREE:-10.0000");
			params.push_back("fa_standard:Narg:LK_DGFREE:-11.0000");
			params.push_back("fa_standard:OH:LK_DGFREE:-6.7700" );
			params.push_back("fa_standard:OW:LK_DGFREE:-6.7700" );
		} else { // talaris
			params.push_back("fa_standard:ONH2:LK_DGFREE:-5.8500");
			params.push_back("fa_standard:NH2O:LK_DGFREE:-7.8000");
			params.push_back("fa_standard:Narg:LK_DGFREE:-10.0000");
			params.push_back("fa_standard:OH:LK_DGFREE:-6.7000");
			params.push_back("fa_standard:OW:LK_DGFREE:-6.7000");
		}

		options[ basic::options::OptionKeys::chemical::set_atom_properties ].value(params);
	}
}

void
init_revert_to_talaris( utility::options::OptionCollection & options ) {
	if ( !(options[ corrections::restore_talaris_behavior ] ||
			options[ mistakes::restore_pre_talaris_2013_behavior ]
			)
			) {
		return;
	}

	// specific for talaris / pre-talaris
	if ( !options[ mistakes::restore_pre_talaris_2013_behavior ] ) {
		// hbond_param_set
		if ( ! options[ score::hbond_params ].user() ) {
			options[ score::hbond_params ].value( "sp2_elec_params" );
		}

		// weights file
		if ( ! options[ score::weights ].user() ) {
			options[ score::weights ].value( "talaris2014.wts" );
		} else {
			TR.Warning << "Flag -restore_talaris_behavior is set but -weights are also specified.  Not changing input weights file!" << std::endl;
		}
	}

	// below are shared b/w talaris and pre-talaris
	// residue_type_set
	// let chemical manager take care of this part instead
	//if ( ! options[ basic::options::OptionKeys::in::file::residue_type_set ].user() ) {
	// options[ basic::options::OptionKeys::in::file::residue_type_set ].value("fa_standard_talaris");
	//}

	// unmodifypot
	if ( ! options[ basic::options::OptionKeys::score::unmodifypot ].user() ) {
		options[ basic::options::OptionKeys::score::unmodifypot ].value(false);
	}

	// attractive hydrogens
	if ( ! options[ basic::options::OptionKeys::score::fa_Hatr ].user() ) {
		options[ basic::options::OptionKeys::score::fa_Hatr ].value(false);
	}

	// sigmoidal dielectric
	if ( ! options[ basic::options::OptionKeys::score::elec_sigmoidal_die ].user() ) {
		options[ basic::options::OptionKeys::score::elec_sigmoidal_die ].value(false);
	}

	// roland's new smoothed p_aa_pp and fa_dun libraries
	if ( ! options[ basic::options::OptionKeys::corrections::shapovalov_lib_fixes_enable ].user() ) {
		options[ basic::options::OptionKeys::corrections::shapovalov_lib_fixes_enable ].value(false);
	}

	// bbdep omega
	if ( ! options[ basic::options::OptionKeys::corrections::score::bbdep_omega ].user() ) {
		options[ basic::options::OptionKeys::corrections::score::bbdep_omega ].value(false);
	}

	// updated fa_elec countpair logic
	if ( ! options[ basic::options::OptionKeys::score::elec_representative_cp ].user() ) {
		options[ basic::options::OptionKeys::score::elec_representative_cp ].value(false);
	}
	if ( ! options[ basic::options::OptionKeys::score::elec_representative_cp_flip ].user() ) {
		options[ basic::options::OptionKeys::score::elec_representative_cp_flip ].value(false);
	}

	// updated cartbonded parameters
	if ( ! options[ basic::options::OptionKeys::score::bonded_params_dir ].user() ) {
		options[ basic::options::OptionKeys::score::bonded_params_dir ].value("scoring/score_functions/bondlength_bondangle_talaris/");
	}
}

void
init_hbond_sp2_correction( utility::options::OptionCollection & options ) {
	if ( options[ corrections::hbond_sp2_correction ] ) {
		TR.Warning << "The -hbond_sp2_correction is deprecated.  The default talaris2013 score function should be relied upon instead." << std::endl;


		// Note, the weight sets that go with the hbond_sp2_correction
		// flag are handled in core/scoring/ScorefunctionFactory.hh

		if ( ! options[ corrections::score::hb_sp2_chipen ].user() ) {
			options[ corrections::score::hb_sp2_chipen ].value( true );
		}

		if ( ! options[ corrections::score::hb_sp2_BAH180_rise ].user() ) {
			options[ corrections::score::hb_sp2_BAH180_rise ].value( 0.75 );
		}

		if ( ! options[ corrections::score::hb_sp2_outer_width ].user() ) {
			options[ corrections::score::hb_sp2_outer_width ].value( 0.357 );
		}

		if ( ! options[ corrections::score::hb_fade_energy ].user() ) {
			options[ corrections::score::hb_fade_energy ].value( true );
		}

		if ( ! options[ corrections::score::hbond_measure_sp3acc_BAH_from_hvy ].user() ) {
			options[ corrections::score::hbond_measure_sp3acc_BAH_from_hvy ].value( true );
		}

		if ( ! options[ corrections::score::lj_hbond_hdis ].user() ) {
			options[ corrections::score::lj_hbond_hdis ].value( 1.75 );
		}

		if ( ! options[ corrections::score::lj_hbond_OH_donor_dis ].user() ) {
			options[ corrections::score::lj_hbond_OH_donor_dis ].value( 2.6 );
		}

		if ( ! options[ corrections::chemical::expand_st_chi2sampling ].user() ) {
			options[ corrections::chemical::expand_st_chi2sampling ].value( true );
		}

		if ( ! options[ score::hbond_params ].user() ) {
			options[ score::hbond_params ].value( "sp2_elec_params" );
		}

		if ( ! options[ score::smooth_fa_elec ].user() ) {
			options[ score::smooth_fa_elec ].value( true );
		}

		if ( ! options[ score::elec_min_dis ].user() ) {
			options[ score::elec_min_dis ].value( 1.6 );
		}

		if ( ! options[ score::elec_r_option ].user() ) {
			options[ score::elec_r_option ].value( false );
		}

		/*
		if ( ! options[ basic::options::OptionKeys::chemical::set_atom_properties ].user() || options[ mistakes::restore_pre_talaris_2013_behavior ] ) {
		std::cout << "overriding sc12 changes to the LK_DGFREEs" << std::endl;
		utility::vector1< std::string > params;
		params.push_back("fa_standard:ONH2:LK_DGFREE:-5.85");
		params.push_back("fa_standard:NH2O:LK_DGFREE:-7.8");
		params.push_back("fa_standard:Narg:LK_DGFREE:-10.0");
		params.push_back("fa_standard:OH:LK_DGFREE:-6.70" );
		options[ basic::options::OptionKeys::chemical::set_atom_properties ].value(params);
		}
		*/

		if ( ! options[ corrections::score::dun10 ].user() ) {
			options[corrections::score::dun10].value( true );
		}

		if ( ! options[ corrections::score::use_bicubic_interpolation ].user() ) {
			options[corrections::score::use_bicubic_interpolation].value( true );
		}

		if ( ! options[ score::weights ].user() ) {
			options[ score::weights ].value( "sp2_correction" );
		}

		if ( options[corrections::correct] ) {
			throw CREATE_EXCEPTION(utility::excn::BadInput,  "The -corrections::hbond_sp2_correction is incompatible with the -corrections:correct flag.");
		}

		if ( options[ corrections::facts_default ] ) {
			throw CREATE_EXCEPTION(utility::excn::BadInput,  "The -corrections::hbond_sp2_correction is incompatible with the -corrections:facts_default flag.");
		}
	}
}

void
init_facts_correction( utility::options::OptionCollection & options ) {
	// corrections specific for FACTS: this has to precede -correct
	if ( options[corrections::facts_default] ) {
		// Call all the corrections from -correct
		options[corrections::correct].value(true);

		/// Below are FACTS specific flags
		// Hbond LJ parameters
		if ( ! options[ corrections::score::lj_hbond_hdis ].user() ) {
			options[corrections::score::lj_hbond_hdis].value(2.3);
		}

		if ( ! options[ corrections::score::lj_hbond_OH_donor_dis ].user() ) {
			options[corrections::score::lj_hbond_OH_donor_dis].value(3.4);
		}

		// FACTS options - these are not fully optimized, can be changed in the future
		if ( ! options[score::facts_GBpair_cut].user() ) {
			options[score::facts_GBpair_cut].value(10.0);
		}

		if ( ! options[score::facts_dshift].user() ) {
			utility::vector1< Real > params;
			params.push_back( 0.0 );
			params.push_back( 1.5 );
			params.push_back( 1.5 );
			params.push_back( 1.5 );
			options[score::facts_dshift].value( params );
		}

		if ( ! options[score::facts_die].user() ) {
			options[score::facts_die].value(1.0);
		}

		if ( ! options[score::facts_kappa].user() ) {
			options[score::facts_kappa].value(12.0);
		}

		if ( ! options[score::facts_plane_to_self].user() ) {
			options[score::facts_plane_to_self].value(true);
		}

		if ( ! options[score::facts_asp_patch].user() ) {
			options[score::facts_asp_patch].value(3);
		}

		// Below are default for sp2correction which are FACTS defaults too
		if ( ! options[ corrections::score::use_bicubic_interpolation ].user() ) {
			options[corrections::score::use_bicubic_interpolation].value(true);
		}

		if ( ! options[score::hbond_params].user() ) {
			options[score::hbond_params ].value( "sp2_params" );
		}

		if ( ! options[corrections::score::hb_sp2_chipen].user() ) {
			options[corrections::score::hb_sp2_chipen].value(true);
		}

		if ( ! options[corrections::score::hbond_measure_sp3acc_BAH_from_hvy].user() ) {
			options[corrections::score::hbond_measure_sp3acc_BAH_from_hvy].value(true);
		}

		if ( ! options[ corrections::score::hb_fade_energy ].user() ) {
			options[ corrections::score::hb_fade_energy ].value( true );
		}

	} // end facts_default
}

//MaximCode:
void
init_correct_correction( utility::options::OptionCollection & options ) {

	// set default corrections
	if ( options[corrections::correct] ) {

		/// Legacy warning
		runtime_assert_string_msg( options[ mistakes::restore_pre_talaris_2013_behavior ](), "\"-correct\" does not behave the way it did before talaris2013 became default, and is not intended for scorefunctions later than score12.  It is not recommended for the talaris2013, talaris2014, ref2015, or beta scorefunctions.  If you wish to use score12 with the \"-correct\" flag, use the \"-restore_pre_talaris_2013_behavior\" flag as well.");

		// Pair energy
		if ( ! options[ corrections::score::no_his_his_pairE ].user() ) {
			options[corrections::score::no_his_his_pairE].value( true );
		}

		// p_aa_pp
		if ( ! options[ corrections::score::p_aa_pp ].user() ) {
			options[corrections::score::p_aa_pp].value( "scoring/score_functions/P_AA_pp/P_AA_pp_08.2009" );
		}
		if ( ! options[ corrections::score::p_aa_pp_nogridshift ].user() ) {
			options[corrections::score::p_aa_pp_nogridshift].value( true );
		}

		//Ramachandran
		if ( ! options[ corrections::score::rama_not_squared ].user() ) {
			options[corrections::score::rama_not_squared].value( true );
		}
		if ( ! options[ corrections::score::rama_map ].user() ) {
			options[ corrections::score::rama_map ].value("scoring/score_functions/rama/Rama09_noEH_kernel25_it08.dat");
		}

		//rotamer library
		if ( ! options[ corrections::score::dun02_file ].user() ) {
			options[corrections::score::dun02_file].value( "rotamer/bbdep02.May.sortlib-correct.12.2010" );
		}

		//icoor
		if ( ! options[ corrections::chemical::icoor_05_2009 ].user() ) {
			options[corrections::chemical::icoor_05_2009].value( true );
		}
		if ( ! options[ score::hbond_params ].user() ) {
			options[score::hbond_params].value( "correct_params" );
		}

		if ( ! options[ corrections::score::ch_o_bond_potential ].user() ) {
			options[ corrections::score::ch_o_bond_potential ].value("scoring/score_functions/carbon_hbond/ch_o_bond_potential_near_min_yf.dat");
		}

		if ( ! options[ score::weights ].user() ) {
			options[ score::weights ].value( "score12_w_corrections" );
		}
	}
}

//MaximCode:
void
init_shapovalov_lib_fixes_enable_correction( utility::options::OptionCollection & options )
{
	// set default corrections
	if ( options[corrections::shapovalov_lib_fixes_enable] ) {
		if ( options[corrections::shapovalov_lib::shap_dun10_enable] ) {
			if ( ! options[ corrections::shapovalov_lib::shap_dun10_dir ].user() ) {
				std::string _smoothingRequsted = options[ corrections::shapovalov_lib::shap_dun10_smooth_level ];
				if ( _smoothingRequsted=="1" || _smoothingRequsted=="lowest_smooth" ) {
					options[ corrections::shapovalov_lib::shap_dun10_dir ].value("rotamer/shapovalov/StpDwn_0-0-0");
					options[ corrections::shapovalov_lib::shap_dun10_smooth_level ].value("1");
				} else if ( _smoothingRequsted=="2" || _smoothingRequsted=="lower_smooth" ) {
					options[ corrections::shapovalov_lib::shap_dun10_dir ].value("rotamer/shapovalov/StpDwn_2-2-2");
					options[ corrections::shapovalov_lib::shap_dun10_smooth_level ].value("2");
				} else if ( _smoothingRequsted=="3" || _smoothingRequsted=="low_smooth" ) {
					options[ corrections::shapovalov_lib::shap_dun10_dir ].value("rotamer/shapovalov/StpDwn_5-5-5");
					options[ corrections::shapovalov_lib::shap_dun10_smooth_level ].value("3");
				} else if ( _smoothingRequsted=="4" || _smoothingRequsted=="average_smooth" ) {
					options[ corrections::shapovalov_lib::shap_dun10_dir ].value("rotamer/shapovalov/StpDwn_10-10-10");
					options[ corrections::shapovalov_lib::shap_dun10_smooth_level ].value("4");
				} else if ( _smoothingRequsted=="5" || _smoothingRequsted=="higher_smooth" ) {
					options[ corrections::shapovalov_lib::shap_dun10_dir ].value("rotamer/shapovalov/StpDwn_20-20-20");
					options[ corrections::shapovalov_lib::shap_dun10_smooth_level ].value("5");
				} else if ( _smoothingRequsted=="6" || _smoothingRequsted=="highest_smooth" ) {
					options[ corrections::shapovalov_lib::shap_dun10_dir ].value("rotamer/shapovalov/StpDwn_25-25-25");
					options[ corrections::shapovalov_lib::shap_dun10_smooth_level ].value("6");
				} else {
					options[ corrections::shapovalov_lib::shap_dun10_dir ].value("rotamer/shapovalov/StpDwn_5-5-5");
					options[ corrections::shapovalov_lib::shap_dun10_smooth_level ].value("3");
				}
			} else {
				options[ corrections::shapovalov_lib::shap_dun10_smooth_level ].value("overriden");
			}
		}


		if ( options[corrections::shapovalov_lib::shap_rama_enable] ) {
			if ( ! options[ corrections::shapovalov_lib::shap_rama_map ].user() ) {
				std::string _smoothingRequsted = options[ corrections::shapovalov_lib::shap_rama_smooth_level ];
				if ( _smoothingRequsted=="1" || _smoothingRequsted=="lowest_smooth" ) {
					options[ corrections::shapovalov_lib::shap_rama_map ].value("scoring/score_functions/rama/shapovalov/kappa75/all.ramaProb");
					options[ corrections::shapovalov_lib::shap_rama_smooth_level ].value("1");
				} else if ( _smoothingRequsted=="2" || _smoothingRequsted=="lower_smooth" ) {
					options[ corrections::shapovalov_lib::shap_rama_map ].value("scoring/score_functions/rama/shapovalov/kappa50/all.ramaProb");
					options[ corrections::shapovalov_lib::shap_rama_smooth_level ].value("2");
				} else if ( _smoothingRequsted=="3" || _smoothingRequsted=="higher_smooth" ) {
					options[ corrections::shapovalov_lib::shap_rama_map ].value("scoring/score_functions/rama/shapovalov/kappa37.5/all.ramaProb");
					options[ corrections::shapovalov_lib::shap_rama_smooth_level ].value("3");
				} else if ( _smoothingRequsted=="4" || _smoothingRequsted=="highest_smooth" ) {
					options[ corrections::shapovalov_lib::shap_rama_map ].value("scoring/score_functions/rama/shapovalov/kappa25/all.ramaProb");
					options[ corrections::shapovalov_lib::shap_rama_smooth_level ].value("4");
				} else {
					options[ corrections::shapovalov_lib::shap_rama_map ].value("scoring/score_functions/rama/shapovalov/kappa25/all.ramaProb");
					options[ corrections::shapovalov_lib::shap_rama_smooth_level ].value("4");
				}
			} else {
				options[ corrections::shapovalov_lib::shap_rama_smooth_level ].value("overriden");
			}
		}


		if ( options[corrections::shapovalov_lib::shap_p_aa_pp_enable] ) {
			if ( ! options[ corrections::shapovalov_lib::shap_p_aa_pp ].user() ) {
				std::string _smoothingRequsted = options[ corrections::shapovalov_lib::shap_p_aa_pp_smooth_level ];
				if ( _smoothingRequsted=="1" || _smoothingRequsted=="low_smooth" ) {
					options[ corrections::shapovalov_lib::shap_p_aa_pp ].value("scoring/score_functions/P_AA_pp/shapovalov/10deg/kappa131/a20.prop");
					options[ corrections::shapovalov_lib::shap_p_aa_pp_smooth_level ].value("1");
				} else if ( _smoothingRequsted=="2" || _smoothingRequsted=="high_smooth" ) {
					options[ corrections::shapovalov_lib::shap_p_aa_pp ].value("scoring/score_functions/P_AA_pp/shapovalov/10deg/kappa50/a20.prop");
					options[ corrections::shapovalov_lib::shap_p_aa_pp_smooth_level ].value("2");
				} else {
					options[ corrections::shapovalov_lib::shap_p_aa_pp ].value("scoring/score_functions/P_AA_pp/shapovalov/10deg/kappa50/a20.prop");
					options[ corrections::shapovalov_lib::shap_p_aa_pp_smooth_level ].value("2");
				}
			} else {
				options[ corrections::shapovalov_lib::shap_p_aa_pp_smooth_level ].value("overriden");
			}
		}


	}
}


void
init_crystal_refinement_correction( utility::options::OptionCollection & options ) {
	//fpd  crystal-refinement specific changes
	if ( options[cryst::crystal_refine] ) {
		// read pdbs properly
		if ( ! options[ in::missing_density_to_jump ].user() ) {
			options[in::missing_density_to_jump].value( true );
		}
		if ( ! options[ in::preserve_crystinfo ].user() ) {
			options[in::preserve_crystinfo].value( true );
		}

		// let centroid structures be scored against crystal data
		if ( ! options[ out::file::no_output_cen ].user() ) {
			options[out::file::no_output_cen].value( true );
		}

		// internal coordinate minimization change: don't dampen RB movement
		if ( ! options[ optimization::scale_rb ].user() ) {
			options[optimization::scale_rb].value( 1.0 );
		}

		// properly read & write glycans
		if ( ! options[ in::include_sugars ].user() ) {
			options[in::include_sugars].value( true );
		}
		if ( ! options[ in::alternate_3_letter_codes ].user() ) {
			options[in::alternate_3_letter_codes].value( "pdb_sugar" );
		}
		if ( ! options[ in::auto_detect_glycan_connections ].user() ) {
			options[in::auto_detect_glycan_connections].value( true );
		}
		if ( ! options[ out::file::write_glycan_pdb_codes ].user() ) {
			options[out::file::write_glycan_pdb_codes].value( true );
		}
		if ( ! options[ score::ideal_sugars ].user() ) {
			options[score::ideal_sugars].value( true );
		}
	}

	// another nonspecific one
	if ( options[edensity::score_symm_complex]() ) {
		options[optimization::old_sym_min].value( true );
	}
}

void
init_beta_correction( utility::options::OptionCollection & options ) {
	// -beta activates most recent changes
	if ( options[corrections::beta]() ) {
		options[ corrections::beta_jan25 ].value(true);
	}
	if ( options[corrections::beta_cart]() ) {
		options[ corrections::beta_jan25_cart ].value(true);
	}

	// if multiple flags specified use the most recent
	if ( options[corrections::beta_jan25]() || options[corrections::beta_jan25_cart]() ) {
		init_beta_jan25_correction( options );
		init_gen_potential_settings( options, false );
	} else if ( options[ corrections::gen_potential ]() ) {
		init_beta_nov16_correction( options, false );  // genpot is "on top" of beta_nov16
		init_gen_potential_settings( options );
	} else if ( options[corrections::beta_nov16]() || options[corrections::beta_nov16_cart]() ) {
		init_beta_nov16_correction( options );
	} else if ( options[corrections::beta_july15]() || options[corrections::beta_july15_cart]() ) {
		init_beta_july15_correction( options );
	}
}


void
init_beta_july15_correction( utility::options::OptionCollection & options ) {
	// atom clone
	if ( ! options[ basic::options::OptionKeys::chemical::clone_atom_types ].user() ) {
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
		options[ basic::options::OptionKeys::chemical::clone_atom_types ].value(params);
	} else {
		TR.Warning << "Flag -beta_july15 is set but -clone_atom_types are also specified.  Not changing atom properties!" << std::endl;
	}

	// atom reassign
	if ( ! options[ basic::options::OptionKeys::chemical::reassign_atom_types ].user() ) {
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
		options[ basic::options::OptionKeys::chemical::reassign_atom_types ].value(params);
	} else {
		TR.Warning << "Flag -beta_july15 is set but -reassign_atom_types are also specified.  Not changing atom properties!" << std::endl;
	}

	// atom properties
	if ( ! options[ basic::options::OptionKeys::chemical::set_atom_properties ].user() ) {
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
		params.push_back( "fa_standard:OW:LJ_RADIUS:1.55048059380817");
		params.push_back( "fa_standard:OW:LJ_WDEPTH:0.160120148452043");
		params.push_back( "fa_standard:OW:LK_DGFREE:-6.38881994292366");
		params.push_back( "fa_standard:OW:LK_VOLUME:10.722");
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
		options[ basic::options::OptionKeys::chemical::set_atom_properties ].value(params);
	} else {
		TR.Warning << "Flag -beta_july15 is set but -set_atom_properties are also specified.  Not changing atom properties!" << std::endl;
	}

	// cart_improper autogen is not necessary with beta_genpot
	if ( ! options[ basic::options::OptionKeys::corrections::score::no_autogen_cart_improper ].user() ) {
		options[ basic::options::OptionKeys::corrections::score::no_autogen_cart_improper ].value( true );
	}

	// hbonds
	if ( ! options[ basic::options::OptionKeys::score::hb_acc_strength ].user() ) {
		utility::vector1< std::string > params;
		params.push_back( "hbacc_AHX:0.94");
		params.push_back( "hbacc_CXA:1.08");
		params.push_back( "hbacc_CXL:1.00");
		params.push_back( "hbacc_HXL:1.00");
		params.push_back( "hbacc_IMD:1.05");
		params.push_back( "hbacc_IME:0.94");
		params.push_back( "hbacc_PBA:0.94");
		params.push_back( "hbacc_GENERIC_SP2SC:1.19");
		params.push_back( "hbacc_GENERIC_SP3SC:1.19");
		params.push_back( "hbacc_GENERIC_RINGSC:1.19");

		options[ basic::options::OptionKeys::score::hb_acc_strength ].value(params);
	} else {
		TR.Warning << "Flag -beta_july15 is set but -hb_acc_strength are also specified.  Not changing atom properties!" << std::endl;
	}
	if ( ! options[ basic::options::OptionKeys::score::hb_don_strength ].user() ) {
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
		params.push_back( "hbdon_GENERIC_SC:1.45");
		options[ basic::options::OptionKeys::score::hb_don_strength ].value(params);
	} else {
		TR.Warning << "Flag -beta_july15 is set but -hb_don_strength are also specified.  Not changing atom properties!" << std::endl;
	}

	// unmodifypot
	if ( ! options[ basic::options::OptionKeys::score::unmodifypot ].user() ) {
		options[ basic::options::OptionKeys::score::unmodifypot ].value(true);
	}

	// sigmoidal dielectric
	if ( ! options[ basic::options::OptionKeys::score::elec_sigmoidal_die ].user() ) {
		options[ basic::options::OptionKeys::score::elec_sigmoidal_die ].value(true);
	}
	if ( ! options[ basic::options::OptionKeys::score::elec_sigmoidal_die_D ].user() ) {
		options[ basic::options::OptionKeys::score::elec_sigmoidal_die_D ].value(137.0);
	}
	if ( ! options[ basic::options::OptionKeys::score::elec_sigmoidal_die_D0 ].user() ) {
		options[ basic::options::OptionKeys::score::elec_sigmoidal_die_D0 ].value(9.75);
	}
	if ( ! options[ basic::options::OptionKeys::score::elec_sigmoidal_die_S ].user() ) {
		options[ basic::options::OptionKeys::score::elec_sigmoidal_die_S ].value(0.44);
	}

	// weights file
	if ( ! options[ score::weights ].user() ) {
		if ( options[corrections::beta_july15_cart]() || options[optimization::nonideal]() ) {
			options[ score::weights ].value( "beta_july15_cart.wts" );
		} else {
			options[ score::weights ].value( "beta_july15.wts" );
		}
	} else {
		TR.Warning << "Flag -beta_july15 is set but -weights are also specified.  Not changing input weights file!" << std::endl;
	}
}

void
init_beta_nov16_correction( utility::options::OptionCollection & options, bool setweight ) {
	// incompatibility check
	if ( options[corrections::score::rama_prepro_steep]() ) {
		utility_exit_with_message("-rama_prepro_steep is incompatible with -beta_nov16.  Use beta_nov15 or remove this flag");
	}

	// atom properties
	if ( ! options[ basic::options::OptionKeys::chemical::set_atom_properties ].user() ) {
		utility::vector1< std::string > params;
		// keep DGFREEs only
		params.push_back("fa_standard:aroC:LK_DGFREE:2.22283");
		params.push_back("fa_standard:CAbb:LK_DGFREE:4.44945");
		params.push_back("fa_standard:CH0:LK_DGFREE:1.24868");
		params.push_back("fa_standard:CH1:LK_DGFREE:-6.49218");
		params.push_back("fa_standard:CH2:LK_DGFREE:-2.55184");
		params.push_back("fa_standard:CH3:LK_DGFREE:7.72716");
		params.push_back("fa_standard:CNH2:LK_DGFREE:3.70334");
		params.push_back("fa_standard:CObb:LK_DGFREE:3.57899");
		params.push_back("fa_standard:COO:LK_DGFREE:-2.50876");
		params.push_back("fa_standard:Narg:LK_DGFREE:-8.69602");
		params.push_back("fa_standard:Nbb:LK_DGFREE:-12.84665");
		params.push_back("fa_standard:NH2O:LK_DGFREE:-7.66671");
		params.push_back("fa_standard:Nhis:LK_DGFREE:-9.72534");
		params.push_back("fa_standard:Nlys:LK_DGFREE:-18.74326");
		params.push_back("fa_standard:Npro:LK_DGFREE:-1.51111");
		params.push_back("fa_standard:Ntrp:LK_DGFREE:-10.64481");
		params.push_back("fa_standard:NtrR:LK_DGFREE:-4.92802");
		params.push_back("fa_standard:OCbb:LK_DGFREE:-9.52921");
		params.push_back("fa_standard:OH:LK_DGFREE:-5.46060");
		params.push_back("fa_standard:OW:LK_DGFREE:-5.46060");
		params.push_back("fa_standard:ONH2:LK_DGFREE:-5.03501");
		params.push_back("fa_standard:OOC:LK_DGFREE:-10.20822");
		params.push_back("fa_standard:S:LK_DGFREE:-4.89802");
		params.push_back("fa_standard:SH1:LK_DGFREE:2.07945");

		//fd & rpav: water
		params.push_back("fa_standard:Opoint:LJ_RADIUS:1.542743");
		params.push_back("fa_standard:Opoint:LJ_WDEPTH:0.161947");
		params.push_back("fa_standard:Opoint:LK_VOLUME:10.8000");
		params.push_back("fa_standard:Owat:LJ_RADIUS:1.542743");
		params.push_back("fa_standard:Owat:LJ_WDEPTH:0.161947");
		params.push_back("fa_standard:Owat:LK_DGFREE:-4.5480");
		params.push_back("fa_standard:Hwat:LJ_RADIUS:0.901681");
		params.push_back("fa_standard:Hwat:LJ_WDEPTH:0.01");

		options[ basic::options::OptionKeys::chemical::set_atom_properties ].value(params);
	} else {
		TR.Warning << "Flag -beta_nov16 is set but -set_atom_properties are also specified.  Not changing atom properties!" << std::endl;
	}

	// atomic charge
	if ( ! options[ basic::options::OptionKeys::chemical::set_atomic_charge ].user() ) {
		utility::vector1< std::string > params;
		params.push_back("fa_standard:ALA:1HB:0.12155");
		params.push_back("fa_standard:ALA:2HB:0.12155");
		params.push_back("fa_standard:ALA:3HB:0.12155");
		params.push_back("fa_standard:ALA:C:0.76511");
		params.push_back("fa_standard:ALA:CA:0.11475");
		params.push_back("fa_standard:ALA:CB:-0.36465");
		params.push_back("fa_standard:ALA:H:0.50819");
		params.push_back("fa_standard:ALA:HA:0.14754");
		params.push_back("fa_standard:ALA:N:-0.77049");
		params.push_back("fa_standard:ALA:O:-0.76511");
		params.push_back("fa_standard:ARG:1HB:0.07619");
		params.push_back("fa_standard:ARG:1HD:0.07619");
		params.push_back("fa_standard:ARG:1HG:0.07619");
		params.push_back("fa_standard:ARG:1HH1:0.29787");
		params.push_back("fa_standard:ARG:1HH2:0.29787");
		params.push_back("fa_standard:ARG:2HB:0.07619");
		params.push_back("fa_standard:ARG:2HD:0.07619");
		params.push_back("fa_standard:ARG:2HG:0.07619");
		params.push_back("fa_standard:ARG:2HH1:0.29787");
		params.push_back("fa_standard:ARG:2HH2:0.29787");
		params.push_back("fa_standard:ARG:C:0.76511");
		params.push_back("fa_standard:ARG:CA:0.11475");
		params.push_back("fa_standard:ARG:CB:-0.08557");
		params.push_back("fa_standard:ARG:CD:0.14210");
		params.push_back("fa_standard:ARG:CG:-0.08557");
		params.push_back("fa_standard:ARG:CZ:0.40572");
		params.push_back("fa_standard:ARG:H:0.50819");
		params.push_back("fa_standard:ARG:HA:0.14754");
		params.push_back("fa_standard:ARG:HE:0.28589");
		params.push_back("fa_standard:ARG:N:-0.77049");
		params.push_back("fa_standard:ARG:NE:-0.39713");
		params.push_back("fa_standard:ARG:NH1:-0.45704");
		params.push_back("fa_standard:ARG:NH2:-0.45704");
		params.push_back("fa_standard:ARG:O:-0.76511");
		params.push_back("fa_standard:ASN:1HB:0.07322");
		params.push_back("fa_standard:ASN:1HD2:0.26032");
		params.push_back("fa_standard:ASN:2HB:0.07322");
		params.push_back("fa_standard:ASN:2HD2:0.24405");
		params.push_back("fa_standard:ASN:C:0.76511");
		params.push_back("fa_standard:ASN:CA:0.11475");
		params.push_back("fa_standard:ASN:CB:-0.14643");
		params.push_back("fa_standard:ASN:CG:0.44743");
		params.push_back("fa_standard:ASN:H:0.50819");
		params.push_back("fa_standard:ASN:HA:0.14754");
		params.push_back("fa_standard:ASN:N:-0.77049");
		params.push_back("fa_standard:ASN:ND2:-0.50437");
		params.push_back("fa_standard:ASN:O:-0.76511");
		params.push_back("fa_standard:ASN:OD1:-0.44743");
		params.push_back("fa_standard:ASP:1HB:0.13276");
		params.push_back("fa_standard:ASP:2HB:0.13276");
		params.push_back("fa_standard:ASP:C:0.76511");
		params.push_back("fa_standard:ASP:CA:0.11475");
		params.push_back("fa_standard:ASP:CB:-0.29888");
		params.push_back("fa_standard:ASP:CG:0.75106");
		params.push_back("fa_standard:ASP:H:0.50819");
		params.push_back("fa_standard:ASP:HA:0.14754");
		params.push_back("fa_standard:ASP:N:-0.77049");
		params.push_back("fa_standard:ASP:O:-0.76511");
		params.push_back("fa_standard:ASP:OD1:-0.85885");
		params.push_back("fa_standard:ASP:OD2:-0.85885");
		params.push_back("fa_standard:CYD:1HB:0.12155");
		params.push_back("fa_standard:CYD:2HB:0.12155");
		params.push_back("fa_standard:CYD:C:0.76511");
		params.push_back("fa_standard:CYD:CA:0.11475");
		params.push_back("fa_standard:CYD:CB:-0.13506");
		params.push_back("fa_standard:CYD:H:0.50819");
		params.push_back("fa_standard:CYD:HA:0.14754");
		params.push_back("fa_standard:CYD:N:-0.77049");
		params.push_back("fa_standard:CYD:O:-0.76511");
		params.push_back("fa_standard:CYD:SG:-0.10804");
		params.push_back("fa_standard:CYS:1HB:0.12155");
		params.push_back("fa_standard:CYS:2HB:0.12155");
		params.push_back("fa_standard:CYS:C:0.76511");
		params.push_back("fa_standard:CYS:CA:0.11475");
		params.push_back("fa_standard:CYS:CB:-0.14856");
		params.push_back("fa_standard:CYS:H:0.50819");
		params.push_back("fa_standard:CYS:HA:0.14754");
		params.push_back("fa_standard:CYS:HG:0.21609");
		params.push_back("fa_standard:CYS:N:-0.77049");
		params.push_back("fa_standard:CYS:O:-0.76511");
		params.push_back("fa_standard:CYS:SG:-0.31063");
		params.push_back("fa_standard:GLN:1HB:0.09669");
		params.push_back("fa_standard:GLN:1HE2:0.34378");
		params.push_back("fa_standard:GLN:1HG:0.09669");
		params.push_back("fa_standard:GLN:2HB:0.09669");
		params.push_back("fa_standard:GLN:2HE2:0.32229");
		params.push_back("fa_standard:GLN:2HG:0.09669");
		params.push_back("fa_standard:GLN:C:0.76511");
		params.push_back("fa_standard:GLN:CA:0.11475");
		params.push_back("fa_standard:GLN:CB:-0.19338");
		params.push_back("fa_standard:GLN:CD:0.59087");
		params.push_back("fa_standard:GLN:CG:-0.19338");
		params.push_back("fa_standard:GLN:H:0.50819");
		params.push_back("fa_standard:GLN:HA:0.14754");
		params.push_back("fa_standard:GLN:N:-0.77049");
		params.push_back("fa_standard:GLN:NE2:-0.66607");
		params.push_back("fa_standard:GLN:O:-0.76511");
		params.push_back("fa_standard:GLN:OE1:-0.59087");
		params.push_back("fa_standard:GLU:1HB:0.05862");
		params.push_back("fa_standard:GLU:1HG:0.05862");
		params.push_back("fa_standard:GLU:2HB:0.05862");
		params.push_back("fa_standard:GLU:2HG:0.05862");
		params.push_back("fa_standard:GLU:C:0.76511");
		params.push_back("fa_standard:GLU:CA:0.11475");
		params.push_back("fa_standard:GLU:CB:-0.16925");
		params.push_back("fa_standard:GLU:CD:0.50590");
		params.push_back("fa_standard:GLU:CG:-0.25364");
		params.push_back("fa_standard:GLU:H:0.50819");
		params.push_back("fa_standard:GLU:HA:0.14754");
		params.push_back("fa_standard:GLU:N:-0.77049");
		params.push_back("fa_standard:GLU:O:-0.76511");
		params.push_back("fa_standard:GLU:OE1:-0.65874");
		params.push_back("fa_standard:GLU:OE2:-0.65874");
		params.push_back("fa_standard:GLY:1HA:0.14754");
		params.push_back("fa_standard:GLY:2HA:0.14754");
		params.push_back("fa_standard:GLY:C:0.76511");
		params.push_back("fa_standard:GLY:CA:-0.03279");
		params.push_back("fa_standard:GLY:H:0.50819");
		params.push_back("fa_standard:GLY:N:-0.77049");
		params.push_back("fa_standard:GLY:O:-0.76511");
		params.push_back("fa_standard:HIS_D:1HB:0.08456");
		params.push_back("fa_standard:HIS_D:2HB:0.08456");
		params.push_back("fa_standard:HIS_D:C:0.76511");
		params.push_back("fa_standard:HIS_D:CA:0.11475");
		params.push_back("fa_standard:HIS_D:CB:-0.08456");
		params.push_back("fa_standard:HIS_D:CD2:0.20671");
		params.push_back("fa_standard:HIS_D:CE1:0.23490");
		params.push_back("fa_standard:HIS_D:CG:-0.04698");
		params.push_back("fa_standard:HIS_D:H:0.50819");
		params.push_back("fa_standard:HIS_D:HA:0.14754");
		params.push_back("fa_standard:HIS_D:HD1:0.30067");
		params.push_back("fa_standard:HIS_D:HD2:0.09396");
		params.push_back("fa_standard:HIS_D:HE1:0.12215");
		params.push_back("fa_standard:HIS_D:N:-0.77049");
		params.push_back("fa_standard:HIS_D:ND1:-0.33825");
		params.push_back("fa_standard:HIS_D:NE2:-0.65772");
		params.push_back("fa_standard:HIS_D:O:-0.76511");
		params.push_back("fa_standard:HIS:1HB:0.08456");
		params.push_back("fa_standard:HIS:2HB:0.08456");
		params.push_back("fa_standard:HIS:C:0.76511");
		params.push_back("fa_standard:HIS:CA:0.11475");
		params.push_back("fa_standard:HIS:CB:-0.07517");
		params.push_back("fa_standard:HIS:CD2:-0.04698");
		params.push_back("fa_standard:HIS:CE1:0.23490");
		params.push_back("fa_standard:HIS:CG:0.20671");
		params.push_back("fa_standard:HIS:H:0.50819");
		params.push_back("fa_standard:HIS:HA:0.14754");
		params.push_back("fa_standard:HIS:HD2:0.08456");
		params.push_back("fa_standard:HIS:HE1:0.12215");
		params.push_back("fa_standard:HIS:HE2:0.30067");
		params.push_back("fa_standard:HIS:N:-0.77049");
		params.push_back("fa_standard:HIS:ND1:-0.65772");
		params.push_back("fa_standard:HIS:NE2:-0.33825");
		params.push_back("fa_standard:HIS:O:-0.76511");
		params.push_back("fa_standard:ILE:1HD1:0.12155");
		params.push_back("fa_standard:ILE:1HG1:0.12155");
		params.push_back("fa_standard:ILE:1HG2:0.12155");
		params.push_back("fa_standard:ILE:2HD1:0.12155");
		params.push_back("fa_standard:ILE:2HG1:0.12155");
		params.push_back("fa_standard:ILE:2HG2:0.12155");
		params.push_back("fa_standard:ILE:3HD1:0.12155");
		params.push_back("fa_standard:ILE:3HG2:0.12155");
		params.push_back("fa_standard:ILE:C:0.76511");
		params.push_back("fa_standard:ILE:CA:0.11475");
		params.push_back("fa_standard:ILE:CB:-0.12155");
		params.push_back("fa_standard:ILE:CD1:-0.36465");
		params.push_back("fa_standard:ILE:CG1:-0.24310");
		params.push_back("fa_standard:ILE:CG2:-0.36465");
		params.push_back("fa_standard:ILE:H:0.50819");
		params.push_back("fa_standard:ILE:HA:0.14754");
		params.push_back("fa_standard:ILE:HB:0.12155");
		params.push_back("fa_standard:ILE:N:-0.77049");
		params.push_back("fa_standard:ILE:O:-0.76511");
		params.push_back("fa_standard:LEU:1HB:0.12155");
		params.push_back("fa_standard:LEU:1HD1:0.12155");
		params.push_back("fa_standard:LEU:1HD2:0.12155");
		params.push_back("fa_standard:LEU:2HB:0.12155");
		params.push_back("fa_standard:LEU:2HD1:0.12155");
		params.push_back("fa_standard:LEU:2HD2:0.12155");
		params.push_back("fa_standard:LEU:3HD1:0.12155");
		params.push_back("fa_standard:LEU:3HD2:0.12155");
		params.push_back("fa_standard:LEU:C:0.76511");
		params.push_back("fa_standard:LEU:CA:0.11475");
		params.push_back("fa_standard:LEU:CB:-0.24310");
		params.push_back("fa_standard:LEU:CD1:-0.36465");
		params.push_back("fa_standard:LEU:CD2:-0.36465");
		params.push_back("fa_standard:LEU:CG:-0.12155");
		params.push_back("fa_standard:LEU:H:0.50819");
		params.push_back("fa_standard:LEU:HA:0.14754");
		params.push_back("fa_standard:LEU:HG:0.12155");
		params.push_back("fa_standard:LEU:N:-0.77049");
		params.push_back("fa_standard:LEU:O:-0.76511");
		params.push_back("fa_standard:LYS:1HB:0.08943");
		params.push_back("fa_standard:LYS:1HD:0.08943");
		params.push_back("fa_standard:LYS:1HE:0.05026");
		params.push_back("fa_standard:LYS:1HG:0.08943");
		params.push_back("fa_standard:LYS:1HZ:0.32446");
		params.push_back("fa_standard:LYS:2HB:0.08943");
		params.push_back("fa_standard:LYS:2HD:0.08943");
		params.push_back("fa_standard:LYS:2HE:0.05026");
		params.push_back("fa_standard:LYS:2HG:0.08943");
		params.push_back("fa_standard:LYS:2HZ:0.32446");
		params.push_back("fa_standard:LYS:3HZ:0.32446");
		params.push_back("fa_standard:LYS:C:0.76511");
		params.push_back("fa_standard:LYS:CA:0.11475");
		params.push_back("fa_standard:LYS:CB:-0.17497");
		params.push_back("fa_standard:LYS:CD:-0.17497");
		params.push_back("fa_standard:LYS:CE:0.20694");
		params.push_back("fa_standard:LYS:CG:-0.17497");
		params.push_back("fa_standard:LYS:H:0.50819");
		params.push_back("fa_standard:LYS:HA:0.14754");
		params.push_back("fa_standard:LYS:N:-0.77049");
		params.push_back("fa_standard:LYS:NZ:-0.29249");
		params.push_back("fa_standard:LYS:O:-0.76511");
		params.push_back("fa_standard:MET:1HB:0.12155");
		params.push_back("fa_standard:MET:1HE:0.12155");
		params.push_back("fa_standard:MET:1HG:0.12155");
		params.push_back("fa_standard:MET:2HB:0.12155");
		params.push_back("fa_standard:MET:2HE:0.12155");
		params.push_back("fa_standard:MET:2HG:0.12155");
		params.push_back("fa_standard:MET:3HE:0.12155");
		params.push_back("fa_standard:MET:C:0.76511");
		params.push_back("fa_standard:MET:CA:0.11475");
		params.push_back("fa_standard:MET:CB:-0.24310");
		params.push_back("fa_standard:MET:CE:-0.29712");
		params.push_back("fa_standard:MET:CG:-0.18908");
		params.push_back("fa_standard:MET:H:0.50819");
		params.push_back("fa_standard:MET:HA:0.14754");
		params.push_back("fa_standard:MET:N:-0.77049");
		params.push_back("fa_standard:MET:O:-0.76511");
		params.push_back("fa_standard:MET:SD:-0.12155");
		params.push_back("fa_standard:PHE:1HB:0.14323");
		params.push_back("fa_standard:PHE:2HB:0.14323");
		params.push_back("fa_standard:PHE:C:0.76511");
		params.push_back("fa_standard:PHE:CA:0.11475");
		params.push_back("fa_standard:PHE:CB:-0.28647");
		params.push_back("fa_standard:PHE:CD1:-0.18302");
		params.push_back("fa_standard:PHE:CD2:-0.18302");
		params.push_back("fa_standard:PHE:CE1:-0.18302");
		params.push_back("fa_standard:PHE:CE2:-0.18302");
		params.push_back("fa_standard:PHE:CG:0.00000");
		params.push_back("fa_standard:PHE:CZ:-0.18302");
		params.push_back("fa_standard:PHE:H:0.50819");
		params.push_back("fa_standard:PHE:HA:0.14754");
		params.push_back("fa_standard:PHE:HD1:0.18302");
		params.push_back("fa_standard:PHE:HD2:0.18302");
		params.push_back("fa_standard:PHE:HE1:0.18302");
		params.push_back("fa_standard:PHE:HE2:0.18302");
		params.push_back("fa_standard:PHE:HZ:0.18302");
		params.push_back("fa_standard:PHE:N:-0.77049");
		params.push_back("fa_standard:PHE:O:-0.76511");
		params.push_back("fa_standard:PRO:1HB:0.14754");
		params.push_back("fa_standard:PRO:1HD:0.14754");
		params.push_back("fa_standard:PRO:1HG:0.14754");
		params.push_back("fa_standard:PRO:2HB:0.14754");
		params.push_back("fa_standard:PRO:2HD:0.14754");
		params.push_back("fa_standard:PRO:2HG:0.14754");
		params.push_back("fa_standard:PRO:C:0.76511");
		params.push_back("fa_standard:PRO:CA:0.03279");
		params.push_back("fa_standard:PRO:CB:-0.29508");
		params.push_back("fa_standard:PRO:CG:-0.29508");
		params.push_back("fa_standard:PRO:HA:0.14754");
		params.push_back("fa_standard:PRO:N:-0.47541");
		params.push_back("fa_standard:PRO:O:-0.76511");
		params.push_back("fa_standard:SER:1HB:0.07337");
		params.push_back("fa_standard:SER:2HB:0.07337");
		params.push_back("fa_standard:SER:C:0.76511");
		params.push_back("fa_standard:SER:CA:0.11475");
		params.push_back("fa_standard:SER:CB:0.04076");
		params.push_back("fa_standard:SER:H:0.50819");
		params.push_back("fa_standard:SER:HA:0.14754");
		params.push_back("fa_standard:SER:HG:0.35053");
		params.push_back("fa_standard:SER:N:-0.77049");
		params.push_back("fa_standard:SER:O:-0.76511");
		params.push_back("fa_standard:SER:OG:-0.53802");
		params.push_back("fa_standard:THR:1HG2:0.08381");
		params.push_back("fa_standard:THR:2HG2:0.08381");
		params.push_back("fa_standard:THR:3HG2:0.08381");
		params.push_back("fa_standard:THR:C:0.76511");
		params.push_back("fa_standard:THR:CA:0.11475");
		params.push_back("fa_standard:THR:CB:0.13037");
		params.push_back("fa_standard:THR:CG2:-0.25143");
		params.push_back("fa_standard:THR:H:0.50819");
		params.push_back("fa_standard:THR:HA:0.14754");
		params.push_back("fa_standard:THR:HB:0.08381");
		params.push_back("fa_standard:THR:HG1:0.40042");
		params.push_back("fa_standard:THR:N:-0.77049");
		params.push_back("fa_standard:THR:O:-0.76511");
		params.push_back("fa_standard:THR:OG1:-0.61459");
		params.push_back("fa_standard:TRP:1HB:0.10831");
		params.push_back("fa_standard:TRP:2HB:0.10831");
		params.push_back("fa_standard:TRP:C:0.76511");
		params.push_back("fa_standard:TRP:CA:0.11475");
		params.push_back("fa_standard:TRP:CB:-0.21662");
		params.push_back("fa_standard:TRP:CD1:0.04212");
		params.push_back("fa_standard:TRP:CD2:-0.02407");
		params.push_back("fa_standard:TRP:CE2:0.15644");
		params.push_back("fa_standard:TRP:CE3:-0.13839");
		params.push_back("fa_standard:TRP:CG:-0.03610");
		params.push_back("fa_standard:TRP:CH2:-0.13839");
		params.push_back("fa_standard:TRP:CZ2:-0.13839");
		params.push_back("fa_standard:TRP:CZ3:-0.13839");
		params.push_back("fa_standard:TRP:H:0.50819");
		params.push_back("fa_standard:TRP:HA:0.14754");
		params.push_back("fa_standard:TRP:HD1:0.13839");
		params.push_back("fa_standard:TRP:HE1:0.45730");
		params.push_back("fa_standard:TRP:HE3:0.13839");
		params.push_back("fa_standard:TRP:HH2:0.13839");
		params.push_back("fa_standard:TRP:HZ2:0.13839");
		params.push_back("fa_standard:TRP:HZ3:0.13839");
		params.push_back("fa_standard:TRP:N:-0.77049");
		params.push_back("fa_standard:TRP:NE1:-0.73409");
		params.push_back("fa_standard:TRP:O:-0.76511");
		params.push_back("fa_standard:TYR:1HB:0.07353");
		params.push_back("fa_standard:TYR:2HB:0.07353");
		params.push_back("fa_standard:TYR:C:0.76511");
		params.push_back("fa_standard:TYR:CA:0.11475");
		params.push_back("fa_standard:TYR:CB:-0.14706");
		params.push_back("fa_standard:TYR:CD1:-0.09395");
		params.push_back("fa_standard:TYR:CD2:-0.09395");
		params.push_back("fa_standard:TYR:CE1:-0.09395");
		params.push_back("fa_standard:TYR:CE2:-0.09395");
		params.push_back("fa_standard:TYR:CG:0.00000");
		params.push_back("fa_standard:TYR:CZ:0.08987");
		params.push_back("fa_standard:TYR:H:0.50819");
		params.push_back("fa_standard:TYR:HA:0.14754");
		params.push_back("fa_standard:TYR:HD1:0.09395");
		params.push_back("fa_standard:TYR:HD2:0.09395");
		params.push_back("fa_standard:TYR:HE1:0.09395");
		params.push_back("fa_standard:TYR:HE2:0.09395");
		params.push_back("fa_standard:TYR:HH:0.35130");
		params.push_back("fa_standard:TYR:N:-0.77049");
		params.push_back("fa_standard:TYR:O:-0.76511");
		params.push_back("fa_standard:TYR:OH:-0.44117");
		params.push_back("fa_standard:VAL:1HG1:0.12155");
		params.push_back("fa_standard:VAL:1HG2:0.12155");
		params.push_back("fa_standard:VAL:2HG1:0.12155");
		params.push_back("fa_standard:VAL:2HG2:0.12155");
		params.push_back("fa_standard:VAL:3HG1:0.12155");
		params.push_back("fa_standard:VAL:3HG2:0.12155");
		params.push_back("fa_standard:VAL:C:0.76511");
		params.push_back("fa_standard:VAL:CA:0.11475");
		params.push_back("fa_standard:VAL:CB:-0.12155");
		params.push_back("fa_standard:VAL:CG1:-0.36465");
		params.push_back("fa_standard:VAL:CG2:-0.36465");
		params.push_back("fa_standard:VAL:H:0.50819");
		params.push_back("fa_standard:VAL:HA:0.14754");
		params.push_back("fa_standard:VAL:HB:0.12155");
		params.push_back("fa_standard:VAL:N:-0.77049");
		params.push_back("fa_standard:VAL:O:-0.76511");
		params.push_back("fa_standard:DAL:1HB:0.12155");
		params.push_back("fa_standard:DAL:2HB:0.12155");
		params.push_back("fa_standard:DAL:3HB:0.12155");
		params.push_back("fa_standard:DAL:C:0.76511");
		params.push_back("fa_standard:DAL:CA:0.11475");
		params.push_back("fa_standard:DAL:CB:-0.36465");
		params.push_back("fa_standard:DAL:H:0.50819");
		params.push_back("fa_standard:DAL:HA:0.14754");
		params.push_back("fa_standard:DAL:N:-0.77049");
		params.push_back("fa_standard:DAL:O:-0.76511");
		params.push_back("fa_standard:DAR:1HB:0.07619");
		params.push_back("fa_standard:DAR:1HD:0.07619");
		params.push_back("fa_standard:DAR:1HG:0.07619");
		params.push_back("fa_standard:DAR:1HH1:0.29787");
		params.push_back("fa_standard:DAR:1HH2:0.29787");
		params.push_back("fa_standard:DAR:2HB:0.07619");
		params.push_back("fa_standard:DAR:2HD:0.07619");
		params.push_back("fa_standard:DAR:2HG:0.07619");
		params.push_back("fa_standard:DAR:2HH1:0.29787");
		params.push_back("fa_standard:DAR:2HH2:0.29787");
		params.push_back("fa_standard:DAR:C:0.76511");
		params.push_back("fa_standard:DAR:CA:0.11475");
		params.push_back("fa_standard:DAR:CB:-0.08557");
		params.push_back("fa_standard:DAR:CD:0.14210");
		params.push_back("fa_standard:DAR:CG:-0.08557");
		params.push_back("fa_standard:DAR:CZ:0.40572");
		params.push_back("fa_standard:DAR:H:0.50819");
		params.push_back("fa_standard:DAR:HA:0.14754");
		params.push_back("fa_standard:DAR:HE:0.28589");
		params.push_back("fa_standard:DAR:N:-0.77049");
		params.push_back("fa_standard:DAR:NE:-0.39713");
		params.push_back("fa_standard:DAR:NH1:-0.45704");
		params.push_back("fa_standard:DAR:NH2:-0.45704");
		params.push_back("fa_standard:DAR:O:-0.76511");
		params.push_back("fa_standard:DAN:1HB:0.07322");
		params.push_back("fa_standard:DAN:1HD2:0.26032");
		params.push_back("fa_standard:DAN:2HB:0.07322");
		params.push_back("fa_standard:DAN:2HD2:0.24405");
		params.push_back("fa_standard:DAN:C:0.76511");
		params.push_back("fa_standard:DAN:CA:0.11475");
		params.push_back("fa_standard:DAN:CB:-0.14643");
		params.push_back("fa_standard:DAN:CG:0.44743");
		params.push_back("fa_standard:DAN:H:0.50819");
		params.push_back("fa_standard:DAN:HA:0.14754");
		params.push_back("fa_standard:DAN:N:-0.77049");
		params.push_back("fa_standard:DAN:ND2:-0.50437");
		params.push_back("fa_standard:DAN:O:-0.76511");
		params.push_back("fa_standard:DAN:OD1:-0.44743");
		params.push_back("fa_standard:DAS:1HB:0.13276");
		params.push_back("fa_standard:DAS:2HB:0.13276");
		params.push_back("fa_standard:DAS:C:0.76511");
		params.push_back("fa_standard:DAS:CA:0.11475");
		params.push_back("fa_standard:DAS:CB:-0.29888");
		params.push_back("fa_standard:DAS:CG:0.75106");
		params.push_back("fa_standard:DAS:H:0.50819");
		params.push_back("fa_standard:DAS:HA:0.14754");
		params.push_back("fa_standard:DAS:N:-0.77049");
		params.push_back("fa_standard:DAS:O:-0.76511");
		params.push_back("fa_standard:DAS:OD1:-0.85885");
		params.push_back("fa_standard:DAS:OD2:-0.85885");
		params.push_back("fa_standard:DCS:1HB:0.12155");
		params.push_back("fa_standard:DCS:2HB:0.12155");
		params.push_back("fa_standard:DCS:C:0.76511");
		params.push_back("fa_standard:DCS:CA:0.11475");
		params.push_back("fa_standard:DCS:CB:-0.14856");
		params.push_back("fa_standard:DCS:H:0.50819");
		params.push_back("fa_standard:DCS:HA:0.14754");
		params.push_back("fa_standard:DCS:HG:0.21609");
		params.push_back("fa_standard:DCS:N:-0.77049");
		params.push_back("fa_standard:DCS:O:-0.76511");
		params.push_back("fa_standard:DCS:SG:-0.31063");
		params.push_back("fa_standard:DGN:1HB:0.09669");
		params.push_back("fa_standard:DGN:1HE2:0.34378");
		params.push_back("fa_standard:DGN:1HG:0.09669");
		params.push_back("fa_standard:DGN:2HB:0.09669");
		params.push_back("fa_standard:DGN:2HE2:0.32229");
		params.push_back("fa_standard:DGN:2HG:0.09669");
		params.push_back("fa_standard:DGN:C:0.76511");
		params.push_back("fa_standard:DGN:CA:0.11475");
		params.push_back("fa_standard:DGN:CB:-0.19338");
		params.push_back("fa_standard:DGN:CD:0.59087");
		params.push_back("fa_standard:DGN:CG:-0.19338");
		params.push_back("fa_standard:DGN:H:0.50819");
		params.push_back("fa_standard:DGN:HA:0.14754");
		params.push_back("fa_standard:DGN:N:-0.77049");
		params.push_back("fa_standard:DGN:NE2:-0.66607");
		params.push_back("fa_standard:DGN:O:-0.76511");
		params.push_back("fa_standard:DGN:OE1:-0.59087");
		params.push_back("fa_standard:DGU:1HB:0.05862");
		params.push_back("fa_standard:DGU:1HG:0.05862");
		params.push_back("fa_standard:DGU:2HB:0.05862");
		params.push_back("fa_standard:DGU:2HG:0.05862");
		params.push_back("fa_standard:DGU:C:0.76511");
		params.push_back("fa_standard:DGU:CA:0.11475");
		params.push_back("fa_standard:DGU:CB:-0.16925");
		params.push_back("fa_standard:DGU:CD:0.50590");
		params.push_back("fa_standard:DGU:CG:-0.25364");
		params.push_back("fa_standard:DGU:H:0.50819");
		params.push_back("fa_standard:DGU:HA:0.14754");
		params.push_back("fa_standard:DGU:N:-0.77049");
		params.push_back("fa_standard:DGU:O:-0.76511");
		params.push_back("fa_standard:DGU:OE1:-0.65874");
		params.push_back("fa_standard:DGU:OE2:-0.65874");
		params.push_back("fa_standard:DHI_D:1HB:0.08456");
		params.push_back("fa_standard:DHI_D:2HB:0.08456");
		params.push_back("fa_standard:DHI_D:C:0.76511");
		params.push_back("fa_standard:DHI_D:CA:0.11475");
		params.push_back("fa_standard:DHI_D:CB:-0.08456");
		params.push_back("fa_standard:DHI_D:CD2:0.20671");
		params.push_back("fa_standard:DHI_D:CE1:0.23490");
		params.push_back("fa_standard:DHI_D:CG:-0.04698");
		params.push_back("fa_standard:DHI_D:H:0.50819");
		params.push_back("fa_standard:DHI_D:HA:0.14754");
		params.push_back("fa_standard:DHI_D:HD1:0.30067");
		params.push_back("fa_standard:DHI_D:HD2:0.09396");
		params.push_back("fa_standard:DHI_D:HE1:0.12215");
		params.push_back("fa_standard:DHI_D:N:-0.77049");
		params.push_back("fa_standard:DHI_D:ND1:-0.33825");
		params.push_back("fa_standard:DHI_D:NE2:-0.65772");
		params.push_back("fa_standard:DHI_D:O:-0.76511");
		params.push_back("fa_standard:DHI:1HB:0.08456");
		params.push_back("fa_standard:DHI:2HB:0.08456");
		params.push_back("fa_standard:DHI:C:0.76511");
		params.push_back("fa_standard:DHI:CA:0.11475");
		params.push_back("fa_standard:DHI:CB:-0.07517");
		params.push_back("fa_standard:DHI:CD2:-0.04698");
		params.push_back("fa_standard:DHI:CE1:0.23490");
		params.push_back("fa_standard:DHI:CG:0.20671");
		params.push_back("fa_standard:DHI:H:0.50819");
		params.push_back("fa_standard:DHI:HA:0.14754");
		params.push_back("fa_standard:DHI:HD2:0.08456");
		params.push_back("fa_standard:DHI:HE1:0.12215");
		params.push_back("fa_standard:DHI:HE2:0.30067");
		params.push_back("fa_standard:DHI:N:-0.77049");
		params.push_back("fa_standard:DHI:ND1:-0.65772");
		params.push_back("fa_standard:DHI:NE2:-0.33825");
		params.push_back("fa_standard:DHI:O:-0.76511");
		params.push_back("fa_standard:DIL:1HD1:0.12155");
		params.push_back("fa_standard:DIL:1HG1:0.12155");
		params.push_back("fa_standard:DIL:1HG2:0.12155");
		params.push_back("fa_standard:DIL:2HD1:0.12155");
		params.push_back("fa_standard:DIL:2HG1:0.12155");
		params.push_back("fa_standard:DIL:2HG2:0.12155");
		params.push_back("fa_standard:DIL:3HD1:0.12155");
		params.push_back("fa_standard:DIL:3HG2:0.12155");
		params.push_back("fa_standard:DIL:C:0.76511");
		params.push_back("fa_standard:DIL:CA:0.11475");
		params.push_back("fa_standard:DIL:CB:-0.12155");
		params.push_back("fa_standard:DIL:CD1:-0.36465");
		params.push_back("fa_standard:DIL:CG1:-0.24310");
		params.push_back("fa_standard:DIL:CG2:-0.36465");
		params.push_back("fa_standard:DIL:H:0.50819");
		params.push_back("fa_standard:DIL:HA:0.14754");
		params.push_back("fa_standard:DIL:HB:0.12155");
		params.push_back("fa_standard:DIL:N:-0.77049");
		params.push_back("fa_standard:DIL:O:-0.76511");
		params.push_back("fa_standard:DLE:1HB:0.12155");
		params.push_back("fa_standard:DLE:1HD1:0.12155");
		params.push_back("fa_standard:DLE:1HD2:0.12155");
		params.push_back("fa_standard:DLE:2HB:0.12155");
		params.push_back("fa_standard:DLE:2HD1:0.12155");
		params.push_back("fa_standard:DLE:2HD2:0.12155");
		params.push_back("fa_standard:DLE:3HD1:0.12155");
		params.push_back("fa_standard:DLE:3HD2:0.12155");
		params.push_back("fa_standard:DLE:C:0.76511");
		params.push_back("fa_standard:DLE:CA:0.11475");
		params.push_back("fa_standard:DLE:CB:-0.24310");
		params.push_back("fa_standard:DLE:CD1:-0.36465");
		params.push_back("fa_standard:DLE:CD2:-0.36465");
		params.push_back("fa_standard:DLE:CG:-0.12155");
		params.push_back("fa_standard:DLE:H:0.50819");
		params.push_back("fa_standard:DLE:HA:0.14754");
		params.push_back("fa_standard:DLE:HG:0.12155");
		params.push_back("fa_standard:DLE:N:-0.77049");
		params.push_back("fa_standard:DLE:O:-0.76511");
		params.push_back("fa_standard:DLY:1HB:0.08943");
		params.push_back("fa_standard:DLY:1HD:0.08943");
		params.push_back("fa_standard:DLY:1HE:0.05026");
		params.push_back("fa_standard:DLY:1HG:0.08943");
		params.push_back("fa_standard:DLY:1HZ:0.32446");
		params.push_back("fa_standard:DLY:2HB:0.08943");
		params.push_back("fa_standard:DLY:2HD:0.08943");
		params.push_back("fa_standard:DLY:2HE:0.05026");
		params.push_back("fa_standard:DLY:2HG:0.08943");
		params.push_back("fa_standard:DLY:2HZ:0.32446");
		params.push_back("fa_standard:DLY:3HZ:0.32446");
		params.push_back("fa_standard:DLY:C:0.76511");
		params.push_back("fa_standard:DLY:CA:0.11475");
		params.push_back("fa_standard:DLY:CB:-0.17497");
		params.push_back("fa_standard:DLY:CD:-0.17497");
		params.push_back("fa_standard:DLY:CE:0.20694");
		params.push_back("fa_standard:DLY:CG:-0.17497");
		params.push_back("fa_standard:DLY:H:0.50819");
		params.push_back("fa_standard:DLY:HA:0.14754");
		params.push_back("fa_standard:DLY:N:-0.77049");
		params.push_back("fa_standard:DLY:NZ:-0.29249");
		params.push_back("fa_standard:DLY:O:-0.76511");
		params.push_back("fa_standard:DME:1HB:0.12155");
		params.push_back("fa_standard:DME:1HE:0.12155");
		params.push_back("fa_standard:DME:1HG:0.12155");
		params.push_back("fa_standard:DME:2HB:0.12155");
		params.push_back("fa_standard:DME:2HE:0.12155");
		params.push_back("fa_standard:DME:2HG:0.12155");
		params.push_back("fa_standard:DME:3HE:0.12155");
		params.push_back("fa_standard:DME:C:0.76511");
		params.push_back("fa_standard:DME:CA:0.11475");
		params.push_back("fa_standard:DME:CB:-0.24310");
		params.push_back("fa_standard:DME:CE:-0.29712");
		params.push_back("fa_standard:DME:CG:-0.18908");
		params.push_back("fa_standard:DME:H:0.50819");
		params.push_back("fa_standard:DME:HA:0.14754");
		params.push_back("fa_standard:DME:N:-0.77049");
		params.push_back("fa_standard:DME:O:-0.76511");
		params.push_back("fa_standard:DME:SD:-0.12155");
		params.push_back("fa_standard:DPH:1HB:0.14323");
		params.push_back("fa_standard:DPH:2HB:0.14323");
		params.push_back("fa_standard:DPH:C:0.76511");
		params.push_back("fa_standard:DPH:CA:0.11475");
		params.push_back("fa_standard:DPH:CB:-0.28647");
		params.push_back("fa_standard:DPH:CD1:-0.18302");
		params.push_back("fa_standard:DPH:CD2:-0.18302");
		params.push_back("fa_standard:DPH:CE1:-0.18302");
		params.push_back("fa_standard:DPH:CE2:-0.18302");
		params.push_back("fa_standard:DPH:CG:0.00000");
		params.push_back("fa_standard:DPH:CZ:-0.18302");
		params.push_back("fa_standard:DPH:H:0.50819");
		params.push_back("fa_standard:DPH:HA:0.14754");
		params.push_back("fa_standard:DPH:HD1:0.18302");
		params.push_back("fa_standard:DPH:HD2:0.18302");
		params.push_back("fa_standard:DPH:HE1:0.18302");
		params.push_back("fa_standard:DPH:HE2:0.18302");
		params.push_back("fa_standard:DPH:HZ:0.18302");
		params.push_back("fa_standard:DPH:N:-0.77049");
		params.push_back("fa_standard:DPH:O:-0.76511");
		params.push_back("fa_standard:DPR:1HB:0.14754");
		params.push_back("fa_standard:DPR:1HD:0.14754");
		params.push_back("fa_standard:DPR:1HG:0.14754");
		params.push_back("fa_standard:DPR:2HB:0.14754");
		params.push_back("fa_standard:DPR:2HD:0.14754");
		params.push_back("fa_standard:DPR:2HG:0.14754");
		params.push_back("fa_standard:DPR:C:0.76511");
		params.push_back("fa_standard:DPR:CA:0.03279");
		params.push_back("fa_standard:DPR:CB:-0.29508");
		params.push_back("fa_standard:DPR:CG:-0.29508");
		params.push_back("fa_standard:DPR:HA:0.14754");
		params.push_back("fa_standard:DPR:N:-0.47541");
		params.push_back("fa_standard:DPR:O:-0.76511");
		params.push_back("fa_standard:DSE:1HB:0.07337");
		params.push_back("fa_standard:DSE:2HB:0.07337");
		params.push_back("fa_standard:DSE:C:0.76511");
		params.push_back("fa_standard:DSE:CA:0.11475");
		params.push_back("fa_standard:DSE:CB:0.04076");
		params.push_back("fa_standard:DSE:H:0.50819");
		params.push_back("fa_standard:DSE:HA:0.14754");
		params.push_back("fa_standard:DSE:HG:0.35053");
		params.push_back("fa_standard:DSE:N:-0.77049");
		params.push_back("fa_standard:DSE:O:-0.76511");
		params.push_back("fa_standard:DSE:OG:-0.53802");
		params.push_back("fa_standard:DTH:1HG2:0.08381");
		params.push_back("fa_standard:DTH:2HG2:0.08381");
		params.push_back("fa_standard:DTH:3HG2:0.08381");
		params.push_back("fa_standard:DTH:C:0.76511");
		params.push_back("fa_standard:DTH:CA:0.11475");
		params.push_back("fa_standard:DTH:CB:0.13037");
		params.push_back("fa_standard:DTH:CG2:-0.25143");
		params.push_back("fa_standard:DTH:H:0.50819");
		params.push_back("fa_standard:DTH:HA:0.14754");
		params.push_back("fa_standard:DTH:HB:0.08381");
		params.push_back("fa_standard:DTH:HG1:0.40042");
		params.push_back("fa_standard:DTH:N:-0.77049");
		params.push_back("fa_standard:DTH:O:-0.76511");
		params.push_back("fa_standard:DTH:OG1:-0.61459");
		params.push_back("fa_standard:DTR:1HB:0.10831");
		params.push_back("fa_standard:DTR:2HB:0.10831");
		params.push_back("fa_standard:DTR:C:0.76511");
		params.push_back("fa_standard:DTR:CA:0.11475");
		params.push_back("fa_standard:DTR:CB:-0.21662");
		params.push_back("fa_standard:DTR:CD1:0.04212");
		params.push_back("fa_standard:DTR:CD2:-0.02407");
		params.push_back("fa_standard:DTR:CE2:0.15644");
		params.push_back("fa_standard:DTR:CE3:-0.13839");
		params.push_back("fa_standard:DTR:CG:-0.03610");
		params.push_back("fa_standard:DTR:CH2:-0.13839");
		params.push_back("fa_standard:DTR:CZ2:-0.13839");
		params.push_back("fa_standard:DTR:CZ3:-0.13839");
		params.push_back("fa_standard:DTR:H:0.50819");
		params.push_back("fa_standard:DTR:HA:0.14754");
		params.push_back("fa_standard:DTR:HD1:0.13839");
		params.push_back("fa_standard:DTR:HE1:0.45730");
		params.push_back("fa_standard:DTR:HE3:0.13839");
		params.push_back("fa_standard:DTR:HH2:0.13839");
		params.push_back("fa_standard:DTR:HZ2:0.13839");
		params.push_back("fa_standard:DTR:HZ3:0.13839");
		params.push_back("fa_standard:DTR:N:-0.77049");
		params.push_back("fa_standard:DTR:NE1:-0.73409");
		params.push_back("fa_standard:DTR:O:-0.76511");
		params.push_back("fa_standard:DTY:1HB:0.07353");
		params.push_back("fa_standard:DTY:2HB:0.07353");
		params.push_back("fa_standard:DTY:C:0.76511");
		params.push_back("fa_standard:DTY:CA:0.11475");
		params.push_back("fa_standard:DTY:CB:-0.14706");
		params.push_back("fa_standard:DTY:CD1:-0.09395");
		params.push_back("fa_standard:DTY:CD2:-0.09395");
		params.push_back("fa_standard:DTY:CE1:-0.09395");
		params.push_back("fa_standard:DTY:CE2:-0.09395");
		params.push_back("fa_standard:DTY:CG:0.00000");
		params.push_back("fa_standard:DTY:CZ:0.08987");
		params.push_back("fa_standard:DTY:H:0.50819");
		params.push_back("fa_standard:DTY:HA:0.14754");
		params.push_back("fa_standard:DTY:HD1:0.09395");
		params.push_back("fa_standard:DTY:HD2:0.09395");
		params.push_back("fa_standard:DTY:HE1:0.09395");
		params.push_back("fa_standard:DTY:HE2:0.09395");
		params.push_back("fa_standard:DTY:HH:0.35130");
		params.push_back("fa_standard:DTY:N:-0.77049");
		params.push_back("fa_standard:DTY:O:-0.76511");
		params.push_back("fa_standard:DTY:OH:-0.44117");
		params.push_back("fa_standard:DVA:1HG1:0.12155");
		params.push_back("fa_standard:DVA:1HG2:0.12155");
		params.push_back("fa_standard:DVA:2HG1:0.12155");
		params.push_back("fa_standard:DVA:2HG2:0.12155");
		params.push_back("fa_standard:DVA:3HG1:0.12155");
		params.push_back("fa_standard:DVA:3HG2:0.12155");
		params.push_back("fa_standard:DVA:C:0.76511");
		params.push_back("fa_standard:DVA:CA:0.11475");
		params.push_back("fa_standard:DVA:CB:-0.12155");
		params.push_back("fa_standard:DVA:CG1:-0.36465");
		params.push_back("fa_standard:DVA:CG2:-0.36465");
		params.push_back("fa_standard:DVA:H:0.50819");
		params.push_back("fa_standard:DVA:HA:0.14754");
		params.push_back("fa_standard:DVA:HB:0.12155");
		params.push_back("fa_standard:DVA:N:-0.77049");
		params.push_back("fa_standard:DVA:O:-0.76511");

		// fd & rpav water
		params.push_back("fa_standard:HOH:O:-0.2896");
		params.push_back("fa_standard:HOH:H1:0.1448");
		params.push_back("fa_standard:HOH:H2:0.1448");

		options[ basic::options::OptionKeys::chemical::set_atomic_charge ].value(params);
	} else {
		TR.Warning << "Flag -beta_nov16 is set but -set_atomic_charge is also specified.  Not changing atom properties!" << std::endl;
	}

	// patch charges
	if ( ! options[ basic::options::OptionKeys::chemical::set_patch_atomic_charge ].user() ) {
		utility::vector1< std::string > params;
		params.push_back("fa_standard:ALA:CtermProteinFull:C:0.35448");
		params.push_back("fa_standard:ALA:CtermProteinFull:O:-0.67724");
		params.push_back("fa_standard:ALA:CtermProteinFull:OXT:-0.67724");
		params.push_back("fa_standard:ALA:NtermProteinFull:1H:0.26154");
		params.push_back("fa_standard:ALA:NtermProteinFull:2H:0.26154");
		params.push_back("fa_standard:ALA:NtermProteinFull:3H:0.26154");
		params.push_back("fa_standard:ALA:NtermProteinFull:CA:0.19184");
		params.push_back("fa_standard:ALA:NtermProteinFull:HA:0.12794");
		params.push_back("fa_standard:ALA:NtermProteinFull:N:-0.10440");
		params.push_back("fa_standard:DAL:CtermProteinFull:C:0.35448");
		params.push_back("fa_standard:DAL:CtermProteinFull:O:-0.67724");
		params.push_back("fa_standard:DAL:CtermProteinFull:OXT:-0.67724");
		params.push_back("fa_standard:DAL:NtermProteinFull:1H:0.26154");
		params.push_back("fa_standard:DAL:NtermProteinFull:2H:0.26154");
		params.push_back("fa_standard:DAL:NtermProteinFull:3H:0.26154");
		params.push_back("fa_standard:DAL:NtermProteinFull:CA:0.19184");
		params.push_back("fa_standard:DAL:NtermProteinFull:HA:0.12794");
		params.push_back("fa_standard:DAL:NtermProteinFull:N:-0.10440");
		params.push_back("fa_standard:GLY:NtermProteinFull:1H:0.25156");
		params.push_back("fa_standard:GLY:NtermProteinFull:1HA:0.11215");
		params.push_back("fa_standard:GLY:NtermProteinFull:2H:0.25156");
		params.push_back("fa_standard:GLY:NtermProteinFull:2HA:0.11215");
		params.push_back("fa_standard:GLY:NtermProteinFull:3H:0.25156");
		params.push_back("fa_standard:GLY:NtermProteinFull:CA:0.13539");
		params.push_back("fa_standard:GLY:NtermProteinFull:N:-0.11438");
		params.push_back("fa_standard:PRO:NtermProteinFull:1H:0.17880");
		params.push_back("fa_standard:PRO:NtermProteinFull:1HB:0.09167");
		params.push_back("fa_standard:PRO:NtermProteinFull:1HD:0.09167");
		params.push_back("fa_standard:PRO:NtermProteinFull:1HG:0.09167");
		params.push_back("fa_standard:PRO:NtermProteinFull:2H:0.17880");
		params.push_back("fa_standard:PRO:NtermProteinFull:2HB:0.09167");
		params.push_back("fa_standard:PRO:NtermProteinFull:2HD:0.09167");
		params.push_back("fa_standard:PRO:NtermProteinFull:2HG:0.09167");
		params.push_back("fa_standard:PRO:NtermProteinFull:CA:0.13233");
		params.push_back("fa_standard:PRO:NtermProteinFull:CB:-0.06516");
		params.push_back("fa_standard:PRO:NtermProteinFull:CG:-0.06516");
		params.push_back("fa_standard:PRO:NtermProteinFull:HA:0.09167");
		params.push_back("fa_standard:PRO:NtermProteinFull:N:-0.00127");
		params.push_back("fa_standard:DPR:NtermProteinFull:1H:0.17880");
		params.push_back("fa_standard:DPR:NtermProteinFull:1HB:0.09167");
		params.push_back("fa_standard:DPR:NtermProteinFull:1HD:0.09167");
		params.push_back("fa_standard:DPR:NtermProteinFull:1HG:0.09167");
		params.push_back("fa_standard:DPR:NtermProteinFull:2H:0.17880");
		params.push_back("fa_standard:DPR:NtermProteinFull:2HB:0.09167");
		params.push_back("fa_standard:DPR:NtermProteinFull:2HD:0.09167");
		params.push_back("fa_standard:DPR:NtermProteinFull:2HG:0.09167");
		params.push_back("fa_standard:DPR:NtermProteinFull:CA:0.13233");
		params.push_back("fa_standard:DPR:NtermProteinFull:CB:-0.06516");
		params.push_back("fa_standard:DPR:NtermProteinFull:CG:-0.06516");
		params.push_back("fa_standard:DPR:NtermProteinFull:HA:0.09167");
		params.push_back("fa_standard:DPR:NtermProteinFull:N:-0.00127");
		options[ basic::options::OptionKeys::chemical::set_patch_atomic_charge ].value(params);
	} else {
		TR.Warning << "Flag -beta_nov16 is set but -set_patch_atomic_charge is also specified.  Not changing atom properties!" << std::endl;
	}

	// hbonds
	if ( ! options[ basic::options::OptionKeys::score::hb_acc_strength ].user() ) {
		utility::vector1< std::string > params;
		params.push_back( "hbacc_AHX:1.15");
		params.push_back( "hbacc_CXA:1.21");
		params.push_back( "hbacc_CXL:1.10");
		params.push_back( "hbacc_HXL:1.15");
		params.push_back( "hbacc_IMD:1.13");
		params.push_back( "hbacc_IME:1.17");
		params.push_back( "hbacc_PBA:1.19");
		params.push_back( "hbacc_H2O:1.23");
		params.push_back( "hbacc_GENERIC_SP2SC:1.19");
		params.push_back( "hbacc_GENERIC_SP3SC:1.19");
		params.push_back( "hbacc_GENERIC_RINGSC:1.19");
		options[ basic::options::OptionKeys::score::hb_acc_strength ].value(params);
	} else {
		TR.Warning << "Flag -beta_nov16 is set but -hb_acc_strength are also specified.  Not changing atom properties!" << std::endl;
	}
	if ( ! options[ basic::options::OptionKeys::score::hb_don_strength ].user() ) {
		utility::vector1< std::string > params;
		params.push_back( "hbdon_AHX:1.00");
		params.push_back( "hbdon_AMO:1.17");
		params.push_back( "hbdon_CXA:1.29");
		params.push_back( "hbdon_GDE:1.11");
		params.push_back( "hbdon_GDH:1.11");
		params.push_back( "hbdon_HXL:0.99");
		params.push_back( "hbdon_IMD:1.18");
		params.push_back( "hbdon_IME:1.42");
		params.push_back( "hbdon_IND:1.15");
		params.push_back( "hbdon_PBA:1.45");
		params.push_back( "hbdon_H2O:1.26");
		params.push_back( "hbdon_GENERIC_SC:1.45");
		options[ basic::options::OptionKeys::score::hb_don_strength ].value(params);
	} else {
		TR.Warning << "Flag -beta_nov16 is set but -hb_don_strength are also specified.  Not changing atom properties!" << std::endl;
	}

	// unmodifypot
	if ( ! options[ basic::options::OptionKeys::score::unmodifypot ].user() ) {
		options[ basic::options::OptionKeys::score::unmodifypot ].value(true);
	}

	// attractive hydrogens
	if ( ! options[ basic::options::OptionKeys::score::fa_Hatr ].user() ) {
		options[ basic::options::OptionKeys::score::fa_Hatr ].value(true);
	}

	// new sp3 acceptor
	if ( ! options[ basic::options::OptionKeys::score::hbond_new_sp3_acc ].user() ) {
		options[ basic::options::OptionKeys::score::hbond_new_sp3_acc ].value(true);
	}

	// sigmoidal dielectric
	if ( ! options[ basic::options::OptionKeys::score::elec_sigmoidal_die ].user() ) {
		options[ basic::options::OptionKeys::score::elec_sigmoidal_die ].value(true);
	}
	if ( ! options[ basic::options::OptionKeys::score::elec_sigmoidal_die_D ].user() ) {
		options[ basic::options::OptionKeys::score::elec_sigmoidal_die_D ].value(79.931);
	}
	if ( ! options[ basic::options::OptionKeys::score::elec_sigmoidal_die_D0 ].user() ) {
		options[ basic::options::OptionKeys::score::elec_sigmoidal_die_D0 ].value(6.648);
	}
	if ( ! options[ basic::options::OptionKeys::score::elec_sigmoidal_die_S ].user() ) {
		options[ basic::options::OptionKeys::score::elec_sigmoidal_die_S ].value(0.441546);
	}

	// lkball
	if ( ! options[ basic::options::OptionKeys::dna::specificity::lk_ball_ramp_width_A2 ].user() ) {
		options[ basic::options::OptionKeys::dna::specificity::lk_ball_ramp_width_A2 ].value(3.709);
	}
	if ( ! options[ basic::options::OptionKeys::dna::specificity::lk_ball_overlap_width_A2 ].user() ) {
		options[ basic::options::OptionKeys::dna::specificity::lk_ball_overlap_width_A2 ].value(2.60);
	}
	if ( ! options[ basic::options::OptionKeys::dna::specificity::lk_ball_overlap_gap ].user() ) {
		options[ basic::options::OptionKeys::dna::specificity::lk_ball_overlap_gap ].value(0.50);
	}
	if ( ! options[ basic::options::OptionKeys::dna::specificity::lk_ball_water_fade ].user() ) {
		options[ basic::options::OptionKeys::dna::specificity::lk_ball_water_fade ].value(1.0);
	}
	if ( ! options[ basic::options::OptionKeys::dna::specificity::lk_ball_for_bb ].user() ) {
		options[ basic::options::OptionKeys::dna::specificity::lk_ball_for_bb ].value(true);
	}
	if ( ! options[ basic::options::OptionKeys::dna::specificity::lk_ball_bridge_angle_widthscale ].user() ) {
		options[ basic::options::OptionKeys::dna::specificity::lk_ball_bridge_angle_widthscale ].value(2.8);
	}
	if ( ! options[ basic::options::OptionKeys::dna::specificity::lk_ball_waters_sp2 ].user() ) {
		utility::vector1< core::Real > params;
		params.push_back(2.828);
		params.push_back(134.5);
		params.push_back(0.0);
		params.push_back(2.828);
		params.push_back(134.5);
		params.push_back(180.0);
		options[ basic::options::OptionKeys::dna::specificity::lk_ball_waters_sp2 ].value(params);
	}
	if ( ! options[ basic::options::OptionKeys::dna::specificity::lk_ball_waters_sp3 ].user() ) {
		utility::vector1< core::Real > params;
		params.push_back(2.828);
		params.push_back(109.3);
		params.push_back(120.0);
		params.push_back(2.828);
		params.push_back(109.3);
		params.push_back(240.0);
		options[ basic::options::OptionKeys::dna::specificity::lk_ball_waters_sp3 ].value(params);
	}
	if ( ! options[ basic::options::OptionKeys::dna::specificity::lk_ball_waters_ring ].user() ) {
		utility::vector1< core::Real > params;
		params.push_back(2.828);
		params.push_back(180);
		params.push_back(0);
		options[ basic::options::OptionKeys::dna::specificity::lk_ball_waters_ring ].value(params);
	}
	if ( ! options[ basic::options::OptionKeys::dna::specificity::lk_ball_waters_donor ].user() ) {
		options[ basic::options::OptionKeys::dna::specificity::lk_ball_waters_donor ].value(2.828);
	}


	// roland's new smoothed p_aa_pp and fa_dun libraries
	if ( ! options[ basic::options::OptionKeys::corrections::shapovalov_lib_fixes_enable ].user() ) {
		options[ basic::options::OptionKeys::corrections::shapovalov_lib_fixes_enable ].value(true);
	}
	if ( ! options[ basic::options::OptionKeys::corrections::shapovalov_lib::shap_dun10_enable ].user() ) {
		options[ basic::options::OptionKeys::corrections::shapovalov_lib::shap_dun10_enable ].value(false);
	}
	if ( ! options[ basic::options::OptionKeys::corrections::shapovalov_lib::shap_p_aa_pp_smooth_level ].user() ) {
		options[ basic::options::OptionKeys::corrections::shapovalov_lib::shap_p_aa_pp_smooth_level ].value("1");
	}

	// corrected torsions
	if ( ! options[ basic::options::OptionKeys::corrections::score::rama_pp_map ].user() ) {
		options[ basic::options::OptionKeys::corrections::score::rama_pp_map ].value("scoring/score_functions/rama/fd_beta_nov2016");
	}
	if ( ! options[ basic::options::OptionKeys::corrections::score::dun10_dir ].user() ) {
		options[ basic::options::OptionKeys::corrections::score::dun10_dir ].value("rotamer/beta_nov2016");
	}


	// bbdep omega
	if ( ! options[ basic::options::OptionKeys::corrections::score::bbdep_omega ].user() ) {
		options[ basic::options::OptionKeys::corrections::score::bbdep_omega ].value(true);
	}

	// updated fa_elec countpair logic
	if ( ! options[ basic::options::OptionKeys::score::elec_representative_cp_flip ].user() ) {
		options[ basic::options::OptionKeys::score::elec_representative_cp_flip ].value(true);
	}

	// updated cartbonded parameters
	if ( ! options[ basic::options::OptionKeys::score::bonded_params_dir ].user() ) {
		options[ basic::options::OptionKeys::score::bonded_params_dir ].value("scoring/score_functions/bondlength_bondangle");
	}

	// weights file
	if ( setweight && !options[ score::weights ].user() ) {
		if ( options[corrections::beta_nov16_cart]() || options[optimization::nonideal]() ) {
			options[ score::weights ].value( "beta_nov16_cart.wts" );
		} else {
			options[ score::weights ].value( "beta_nov16.wts" );
		}
	} else {
		TR.Warning << "Flag -beta_nov16 is set but either -weights are also specified or this is being used with the genpot score function.  Not changing input weights file!" << std::endl;
	}

	// protocol option but it should be default ...
	if ( !options[optimization::default_max_cycles]() ) {
		options[ optimization::default_max_cycles ].value( 200 );
	}
}

void
init_beta_jan25_correction( utility::options::OptionCollection & options, bool setweight ) {
	// incompatibility check
	if ( options[corrections::score::rama_prepro_steep]() ) {
		utility_exit_with_message("-rama_prepro_steep is incompatible with -beta_jan25.  Use beta_nov15 or remove this flag");
	}

	// atom properties
	if ( ! options[ basic::options::OptionKeys::chemical::set_atom_properties ].user() ) {
		utility::vector1< std::string > params;
		// keep DGFREEs only
		params.push_back("fa_standard:aroC:LK_DGFREE:2.22283");
		params.push_back("fa_standard:CAbb:LK_DGFREE:4.44945");
		params.push_back("fa_standard:CH0:LK_DGFREE:1.24868");
		params.push_back("fa_standard:CH1:LK_DGFREE:-6.49218");
		params.push_back("fa_standard:CH2:LK_DGFREE:-2.55184");
		params.push_back("fa_standard:CH3:LK_DGFREE:7.72716");
		params.push_back("fa_standard:CNH2:LK_DGFREE:3.70334");
		params.push_back("fa_standard:CObb:LK_DGFREE:3.57899");
		params.push_back("fa_standard:COO:LK_DGFREE:-2.50876");
		params.push_back("fa_standard:Narg:LK_DGFREE:-8.69602");
		params.push_back("fa_standard:Nbb:LK_DGFREE:-12.84665");
		params.push_back("fa_standard:NH2O:LK_DGFREE:-7.66671");
		params.push_back("fa_standard:Nhis:LK_DGFREE:-9.72534");
		params.push_back("fa_standard:Nlys:LK_DGFREE:-18.74326");
		params.push_back("fa_standard:Npro:LK_DGFREE:-1.51111");
		params.push_back("fa_standard:Ntrp:LK_DGFREE:-10.64481");
		params.push_back("fa_standard:NtrR:LK_DGFREE:-4.92802");
		params.push_back("fa_standard:OCbb:LK_DGFREE:-9.52921");
		params.push_back("fa_standard:OH:LK_DGFREE:-5.46060");
		params.push_back("fa_standard:OW:LK_DGFREE:-5.46060");
		params.push_back("fa_standard:ONH2:LK_DGFREE:-5.03501");
		params.push_back("fa_standard:OOC:LK_DGFREE:-10.20822");
		params.push_back("fa_standard:S:LK_DGFREE:-4.89802");
		params.push_back("fa_standard:SH1:LK_DGFREE:2.07945");

		//fd & rpav: water
		params.push_back("fa_standard:Opoint:LJ_RADIUS:1.542743");
		params.push_back("fa_standard:Opoint:LJ_WDEPTH:0.161947");
		params.push_back("fa_standard:Opoint:LK_VOLUME:10.8000");
		params.push_back("fa_standard:Owat:LJ_RADIUS:1.542743");
		params.push_back("fa_standard:Owat:LJ_WDEPTH:0.161947");
		params.push_back("fa_standard:Owat:LK_DGFREE:-4.5480");
		params.push_back("fa_standard:Hwat:LJ_RADIUS:0.901681");
		params.push_back("fa_standard:Hwat:LJ_WDEPTH:0.01");

		// LJ parameters fit in beta_jan25
		params.push_back("fa_standard:CH0:LJ_RADIUS:2.04");
		params.push_back("fa_standard:CH0:LJ_WDEPTH:0.075");
		params.push_back("fa_standard:CH1:LJ_RADIUS:2.01");
		params.push_back("fa_standard:CH1:LJ_WDEPTH:0.062");
		params.push_back("fa_standard:CH2:LJ_RADIUS:2.01");
		params.push_back("fa_standard:CH2:LJ_WDEPTH:0.064");
		params.push_back("fa_standard:CH3:LJ_RADIUS:2.06");
		params.push_back("fa_standard:CH3:LJ_WDEPTH:0.064");
		params.push_back("fa_standard:aroC:LJ_RADIUS:2.02");
		params.push_back("fa_standard:aroC:LJ_WDEPTH:0.070");
		params.push_back("fa_standard:Hapo:LJ_RADIUS:1.45");
		params.push_back("fa_standard:Hapo:LJ_WDEPTH:0.021");
		params.push_back("fa_standard:Haro:LJ_RADIUS:1.41");
		params.push_back("fa_standard:Haro:LJ_WDEPTH:0.021");

		options[ basic::options::OptionKeys::chemical::set_atom_properties ].value(params);
	} else {
		TR.Warning << "Flag -beta_jan25 is set but -set_atom_properties are also specified.  Not changing atom properties!" << std::endl;
	}

	// atomic charge
	if ( ! options[ basic::options::OptionKeys::chemical::set_atomic_charge ].user() ) {
		utility::vector1< std::string > params;
		params.push_back("fa_standard:ALA:1HB:0.12155");
		params.push_back("fa_standard:ALA:2HB:0.12155");
		params.push_back("fa_standard:ALA:3HB:0.12155");
		params.push_back("fa_standard:ALA:C:0.76511");
		params.push_back("fa_standard:ALA:CA:0.11475");
		params.push_back("fa_standard:ALA:CB:-0.36465");
		params.push_back("fa_standard:ALA:H:0.50819");
		params.push_back("fa_standard:ALA:HA:0.14754");
		params.push_back("fa_standard:ALA:N:-0.77049");
		params.push_back("fa_standard:ALA:O:-0.76511");
		params.push_back("fa_standard:ARG:1HB:0.07619");
		params.push_back("fa_standard:ARG:1HD:0.07619");
		params.push_back("fa_standard:ARG:1HG:0.07619");
		params.push_back("fa_standard:ARG:1HH1:0.29787");
		params.push_back("fa_standard:ARG:1HH2:0.29787");
		params.push_back("fa_standard:ARG:2HB:0.07619");
		params.push_back("fa_standard:ARG:2HD:0.07619");
		params.push_back("fa_standard:ARG:2HG:0.07619");
		params.push_back("fa_standard:ARG:2HH1:0.29787");
		params.push_back("fa_standard:ARG:2HH2:0.29787");
		params.push_back("fa_standard:ARG:C:0.76511");
		params.push_back("fa_standard:ARG:CA:0.11475");
		params.push_back("fa_standard:ARG:CB:-0.08557");
		params.push_back("fa_standard:ARG:CD:0.14210");
		params.push_back("fa_standard:ARG:CG:-0.08557");
		params.push_back("fa_standard:ARG:CZ:0.40572");
		params.push_back("fa_standard:ARG:H:0.50819");
		params.push_back("fa_standard:ARG:HA:0.14754");
		params.push_back("fa_standard:ARG:HE:0.28589");
		params.push_back("fa_standard:ARG:N:-0.77049");
		params.push_back("fa_standard:ARG:NE:-0.39713");
		params.push_back("fa_standard:ARG:NH1:-0.45704");
		params.push_back("fa_standard:ARG:NH2:-0.45704");
		params.push_back("fa_standard:ARG:O:-0.76511");
		params.push_back("fa_standard:ASN:1HB:0.07322");
		params.push_back("fa_standard:ASN:1HD2:0.26032");
		params.push_back("fa_standard:ASN:2HB:0.07322");
		params.push_back("fa_standard:ASN:2HD2:0.24405");
		params.push_back("fa_standard:ASN:C:0.76511");
		params.push_back("fa_standard:ASN:CA:0.11475");
		params.push_back("fa_standard:ASN:CB:-0.14643");
		params.push_back("fa_standard:ASN:CG:0.44743");
		params.push_back("fa_standard:ASN:H:0.50819");
		params.push_back("fa_standard:ASN:HA:0.14754");
		params.push_back("fa_standard:ASN:N:-0.77049");
		params.push_back("fa_standard:ASN:ND2:-0.50437");
		params.push_back("fa_standard:ASN:O:-0.76511");
		params.push_back("fa_standard:ASN:OD1:-0.44743");
		params.push_back("fa_standard:ASP:1HB:0.13276");
		params.push_back("fa_standard:ASP:2HB:0.13276");
		params.push_back("fa_standard:ASP:C:0.76511");
		params.push_back("fa_standard:ASP:CA:0.11475");
		params.push_back("fa_standard:ASP:CB:-0.29888");
		params.push_back("fa_standard:ASP:CG:0.75106");
		params.push_back("fa_standard:ASP:H:0.50819");
		params.push_back("fa_standard:ASP:HA:0.14754");
		params.push_back("fa_standard:ASP:N:-0.77049");
		params.push_back("fa_standard:ASP:O:-0.76511");
		params.push_back("fa_standard:ASP:OD1:-0.85885");
		params.push_back("fa_standard:ASP:OD2:-0.85885");
		params.push_back("fa_standard:CYD:1HB:0.12155");
		params.push_back("fa_standard:CYD:2HB:0.12155");
		params.push_back("fa_standard:CYD:C:0.76511");
		params.push_back("fa_standard:CYD:CA:0.11475");
		params.push_back("fa_standard:CYD:CB:-0.13506");
		params.push_back("fa_standard:CYD:H:0.50819");
		params.push_back("fa_standard:CYD:HA:0.14754");
		params.push_back("fa_standard:CYD:N:-0.77049");
		params.push_back("fa_standard:CYD:O:-0.76511");
		params.push_back("fa_standard:CYD:SG:-0.10804");
		params.push_back("fa_standard:CYS:1HB:0.12155");
		params.push_back("fa_standard:CYS:2HB:0.12155");
		params.push_back("fa_standard:CYS:C:0.76511");
		params.push_back("fa_standard:CYS:CA:0.11475");
		params.push_back("fa_standard:CYS:CB:-0.14856");
		params.push_back("fa_standard:CYS:H:0.50819");
		params.push_back("fa_standard:CYS:HA:0.14754");
		params.push_back("fa_standard:CYS:HG:0.21609");
		params.push_back("fa_standard:CYS:N:-0.77049");
		params.push_back("fa_standard:CYS:O:-0.76511");
		params.push_back("fa_standard:CYS:SG:-0.31063");
		params.push_back("fa_standard:GLN:1HB:0.09669");
		params.push_back("fa_standard:GLN:1HE2:0.34378");
		params.push_back("fa_standard:GLN:1HG:0.09669");
		params.push_back("fa_standard:GLN:2HB:0.09669");
		params.push_back("fa_standard:GLN:2HE2:0.32229");
		params.push_back("fa_standard:GLN:2HG:0.09669");
		params.push_back("fa_standard:GLN:C:0.76511");
		params.push_back("fa_standard:GLN:CA:0.11475");
		params.push_back("fa_standard:GLN:CB:-0.19338");
		params.push_back("fa_standard:GLN:CD:0.59087");
		params.push_back("fa_standard:GLN:CG:-0.19338");
		params.push_back("fa_standard:GLN:H:0.50819");
		params.push_back("fa_standard:GLN:HA:0.14754");
		params.push_back("fa_standard:GLN:N:-0.77049");
		params.push_back("fa_standard:GLN:NE2:-0.66607");
		params.push_back("fa_standard:GLN:O:-0.76511");
		params.push_back("fa_standard:GLN:OE1:-0.59087");
		params.push_back("fa_standard:GLU:1HB:0.05862");
		params.push_back("fa_standard:GLU:1HG:0.05862");
		params.push_back("fa_standard:GLU:2HB:0.05862");
		params.push_back("fa_standard:GLU:2HG:0.05862");
		params.push_back("fa_standard:GLU:C:0.76511");
		params.push_back("fa_standard:GLU:CA:0.11475");
		params.push_back("fa_standard:GLU:CB:-0.16925");
		params.push_back("fa_standard:GLU:CD:0.50590");
		params.push_back("fa_standard:GLU:CG:-0.25364");
		params.push_back("fa_standard:GLU:H:0.50819");
		params.push_back("fa_standard:GLU:HA:0.14754");
		params.push_back("fa_standard:GLU:N:-0.77049");
		params.push_back("fa_standard:GLU:O:-0.76511");
		params.push_back("fa_standard:GLU:OE1:-0.65874");
		params.push_back("fa_standard:GLU:OE2:-0.65874");
		params.push_back("fa_standard:GLY:1HA:0.14754");
		params.push_back("fa_standard:GLY:2HA:0.14754");
		params.push_back("fa_standard:GLY:C:0.76511");
		params.push_back("fa_standard:GLY:CA:-0.03279");
		params.push_back("fa_standard:GLY:H:0.50819");
		params.push_back("fa_standard:GLY:N:-0.77049");
		params.push_back("fa_standard:GLY:O:-0.76511");
		params.push_back("fa_standard:HIS_D:1HB:0.08456");
		params.push_back("fa_standard:HIS_D:2HB:0.08456");
		params.push_back("fa_standard:HIS_D:C:0.76511");
		params.push_back("fa_standard:HIS_D:CA:0.11475");
		params.push_back("fa_standard:HIS_D:CB:-0.08456");
		params.push_back("fa_standard:HIS_D:CD2:0.20671");
		params.push_back("fa_standard:HIS_D:CE1:0.23490");
		params.push_back("fa_standard:HIS_D:CG:-0.04698");
		params.push_back("fa_standard:HIS_D:H:0.50819");
		params.push_back("fa_standard:HIS_D:HA:0.14754");
		params.push_back("fa_standard:HIS_D:HD1:0.30067");
		params.push_back("fa_standard:HIS_D:HD2:0.09396");
		params.push_back("fa_standard:HIS_D:HE1:0.12215");
		params.push_back("fa_standard:HIS_D:N:-0.77049");
		params.push_back("fa_standard:HIS_D:ND1:-0.33825");
		params.push_back("fa_standard:HIS_D:NE2:-0.65772");
		params.push_back("fa_standard:HIS_D:O:-0.76511");
		params.push_back("fa_standard:HIS:1HB:0.08456");
		params.push_back("fa_standard:HIS:2HB:0.08456");
		params.push_back("fa_standard:HIS:C:0.76511");
		params.push_back("fa_standard:HIS:CA:0.11475");
		params.push_back("fa_standard:HIS:CB:-0.07517");
		params.push_back("fa_standard:HIS:CD2:-0.04698");
		params.push_back("fa_standard:HIS:CE1:0.23490");
		params.push_back("fa_standard:HIS:CG:0.20671");
		params.push_back("fa_standard:HIS:H:0.50819");
		params.push_back("fa_standard:HIS:HA:0.14754");
		params.push_back("fa_standard:HIS:HD2:0.08456");
		params.push_back("fa_standard:HIS:HE1:0.12215");
		params.push_back("fa_standard:HIS:HE2:0.30067");
		params.push_back("fa_standard:HIS:N:-0.77049");
		params.push_back("fa_standard:HIS:ND1:-0.65772");
		params.push_back("fa_standard:HIS:NE2:-0.33825");
		params.push_back("fa_standard:HIS:O:-0.76511");
		params.push_back("fa_standard:ILE:1HD1:0.12155");
		params.push_back("fa_standard:ILE:1HG1:0.12155");
		params.push_back("fa_standard:ILE:1HG2:0.12155");
		params.push_back("fa_standard:ILE:2HD1:0.12155");
		params.push_back("fa_standard:ILE:2HG1:0.12155");
		params.push_back("fa_standard:ILE:2HG2:0.12155");
		params.push_back("fa_standard:ILE:3HD1:0.12155");
		params.push_back("fa_standard:ILE:3HG2:0.12155");
		params.push_back("fa_standard:ILE:C:0.76511");
		params.push_back("fa_standard:ILE:CA:0.11475");
		params.push_back("fa_standard:ILE:CB:-0.12155");
		params.push_back("fa_standard:ILE:CD1:-0.36465");
		params.push_back("fa_standard:ILE:CG1:-0.24310");
		params.push_back("fa_standard:ILE:CG2:-0.36465");
		params.push_back("fa_standard:ILE:H:0.50819");
		params.push_back("fa_standard:ILE:HA:0.14754");
		params.push_back("fa_standard:ILE:HB:0.12155");
		params.push_back("fa_standard:ILE:N:-0.77049");
		params.push_back("fa_standard:ILE:O:-0.76511");
		params.push_back("fa_standard:LEU:1HB:0.12155");
		params.push_back("fa_standard:LEU:1HD1:0.12155");
		params.push_back("fa_standard:LEU:1HD2:0.12155");
		params.push_back("fa_standard:LEU:2HB:0.12155");
		params.push_back("fa_standard:LEU:2HD1:0.12155");
		params.push_back("fa_standard:LEU:2HD2:0.12155");
		params.push_back("fa_standard:LEU:3HD1:0.12155");
		params.push_back("fa_standard:LEU:3HD2:0.12155");
		params.push_back("fa_standard:LEU:C:0.76511");
		params.push_back("fa_standard:LEU:CA:0.11475");
		params.push_back("fa_standard:LEU:CB:-0.24310");
		params.push_back("fa_standard:LEU:CD1:-0.36465");
		params.push_back("fa_standard:LEU:CD2:-0.36465");
		params.push_back("fa_standard:LEU:CG:-0.12155");
		params.push_back("fa_standard:LEU:H:0.50819");
		params.push_back("fa_standard:LEU:HA:0.14754");
		params.push_back("fa_standard:LEU:HG:0.12155");
		params.push_back("fa_standard:LEU:N:-0.77049");
		params.push_back("fa_standard:LEU:O:-0.76511");
		params.push_back("fa_standard:LYS:1HB:0.08943");
		params.push_back("fa_standard:LYS:1HD:0.08943");
		params.push_back("fa_standard:LYS:1HE:0.05026");
		params.push_back("fa_standard:LYS:1HG:0.08943");
		params.push_back("fa_standard:LYS:1HZ:0.32446");
		params.push_back("fa_standard:LYS:2HB:0.08943");
		params.push_back("fa_standard:LYS:2HD:0.08943");
		params.push_back("fa_standard:LYS:2HE:0.05026");
		params.push_back("fa_standard:LYS:2HG:0.08943");
		params.push_back("fa_standard:LYS:2HZ:0.32446");
		params.push_back("fa_standard:LYS:3HZ:0.32446");
		params.push_back("fa_standard:LYS:C:0.76511");
		params.push_back("fa_standard:LYS:CA:0.11475");
		params.push_back("fa_standard:LYS:CB:-0.17497");
		params.push_back("fa_standard:LYS:CD:-0.17497");
		params.push_back("fa_standard:LYS:CE:0.20694");
		params.push_back("fa_standard:LYS:CG:-0.17497");
		params.push_back("fa_standard:LYS:H:0.50819");
		params.push_back("fa_standard:LYS:HA:0.14754");
		params.push_back("fa_standard:LYS:N:-0.77049");
		params.push_back("fa_standard:LYS:NZ:-0.29249");
		params.push_back("fa_standard:LYS:O:-0.76511");
		params.push_back("fa_standard:MET:1HB:0.12155");
		params.push_back("fa_standard:MET:1HE:0.12155");
		params.push_back("fa_standard:MET:1HG:0.12155");
		params.push_back("fa_standard:MET:2HB:0.12155");
		params.push_back("fa_standard:MET:2HE:0.12155");
		params.push_back("fa_standard:MET:2HG:0.12155");
		params.push_back("fa_standard:MET:3HE:0.12155");
		params.push_back("fa_standard:MET:C:0.76511");
		params.push_back("fa_standard:MET:CA:0.11475");
		params.push_back("fa_standard:MET:CB:-0.24310");
		params.push_back("fa_standard:MET:CE:-0.29712");
		params.push_back("fa_standard:MET:CG:-0.18908");
		params.push_back("fa_standard:MET:H:0.50819");
		params.push_back("fa_standard:MET:HA:0.14754");
		params.push_back("fa_standard:MET:N:-0.77049");
		params.push_back("fa_standard:MET:O:-0.76511");
		params.push_back("fa_standard:MET:SD:-0.12155");
		params.push_back("fa_standard:PHE:1HB:0.14323");
		params.push_back("fa_standard:PHE:2HB:0.14323");
		params.push_back("fa_standard:PHE:C:0.76511");
		params.push_back("fa_standard:PHE:CA:0.11475");
		params.push_back("fa_standard:PHE:CB:-0.28647");
		params.push_back("fa_standard:PHE:CD1:-0.18302");
		params.push_back("fa_standard:PHE:CD2:-0.18302");
		params.push_back("fa_standard:PHE:CE1:-0.18302");
		params.push_back("fa_standard:PHE:CE2:-0.18302");
		params.push_back("fa_standard:PHE:CG:0.00000");
		params.push_back("fa_standard:PHE:CZ:-0.18302");
		params.push_back("fa_standard:PHE:H:0.50819");
		params.push_back("fa_standard:PHE:HA:0.14754");
		params.push_back("fa_standard:PHE:HD1:0.18302");
		params.push_back("fa_standard:PHE:HD2:0.18302");
		params.push_back("fa_standard:PHE:HE1:0.18302");
		params.push_back("fa_standard:PHE:HE2:0.18302");
		params.push_back("fa_standard:PHE:HZ:0.18302");
		params.push_back("fa_standard:PHE:N:-0.77049");
		params.push_back("fa_standard:PHE:O:-0.76511");
		params.push_back("fa_standard:PRO:1HB:0.14754");
		params.push_back("fa_standard:PRO:1HD:0.14754");
		params.push_back("fa_standard:PRO:1HG:0.14754");
		params.push_back("fa_standard:PRO:2HB:0.14754");
		params.push_back("fa_standard:PRO:2HD:0.14754");
		params.push_back("fa_standard:PRO:2HG:0.14754");
		params.push_back("fa_standard:PRO:C:0.76511");
		params.push_back("fa_standard:PRO:CA:0.03279");
		params.push_back("fa_standard:PRO:CB:-0.29508");
		params.push_back("fa_standard:PRO:CG:-0.29508");
		params.push_back("fa_standard:PRO:HA:0.14754");
		params.push_back("fa_standard:PRO:N:-0.47541");
		params.push_back("fa_standard:PRO:O:-0.76511");
		params.push_back("fa_standard:SER:1HB:0.07337");
		params.push_back("fa_standard:SER:2HB:0.07337");
		params.push_back("fa_standard:SER:C:0.76511");
		params.push_back("fa_standard:SER:CA:0.11475");
		params.push_back("fa_standard:SER:CB:0.04076");
		params.push_back("fa_standard:SER:H:0.50819");
		params.push_back("fa_standard:SER:HA:0.14754");
		params.push_back("fa_standard:SER:HG:0.35053");
		params.push_back("fa_standard:SER:N:-0.77049");
		params.push_back("fa_standard:SER:O:-0.76511");
		params.push_back("fa_standard:SER:OG:-0.53802");
		params.push_back("fa_standard:THR:1HG2:0.08381");
		params.push_back("fa_standard:THR:2HG2:0.08381");
		params.push_back("fa_standard:THR:3HG2:0.08381");
		params.push_back("fa_standard:THR:C:0.76511");
		params.push_back("fa_standard:THR:CA:0.11475");
		params.push_back("fa_standard:THR:CB:0.13037");
		params.push_back("fa_standard:THR:CG2:-0.25143");
		params.push_back("fa_standard:THR:H:0.50819");
		params.push_back("fa_standard:THR:HA:0.14754");
		params.push_back("fa_standard:THR:HB:0.08381");
		params.push_back("fa_standard:THR:HG1:0.40042");
		params.push_back("fa_standard:THR:N:-0.77049");
		params.push_back("fa_standard:THR:O:-0.76511");
		params.push_back("fa_standard:THR:OG1:-0.61459");
		params.push_back("fa_standard:TRP:1HB:0.10831");
		params.push_back("fa_standard:TRP:2HB:0.10831");
		params.push_back("fa_standard:TRP:C:0.76511");
		params.push_back("fa_standard:TRP:CA:0.11475");
		params.push_back("fa_standard:TRP:CB:-0.21662");
		params.push_back("fa_standard:TRP:CD1:0.04212");
		params.push_back("fa_standard:TRP:CD2:-0.02407");
		params.push_back("fa_standard:TRP:CE2:0.15644");
		params.push_back("fa_standard:TRP:CE3:-0.13839");
		params.push_back("fa_standard:TRP:CG:-0.03610");
		params.push_back("fa_standard:TRP:CH2:-0.13839");
		params.push_back("fa_standard:TRP:CZ2:-0.13839");
		params.push_back("fa_standard:TRP:CZ3:-0.13839");
		params.push_back("fa_standard:TRP:H:0.50819");
		params.push_back("fa_standard:TRP:HA:0.14754");
		params.push_back("fa_standard:TRP:HD1:0.13839");
		params.push_back("fa_standard:TRP:HE1:0.45730");
		params.push_back("fa_standard:TRP:HE3:0.13839");
		params.push_back("fa_standard:TRP:HH2:0.13839");
		params.push_back("fa_standard:TRP:HZ2:0.13839");
		params.push_back("fa_standard:TRP:HZ3:0.13839");
		params.push_back("fa_standard:TRP:N:-0.77049");
		params.push_back("fa_standard:TRP:NE1:-0.73409");
		params.push_back("fa_standard:TRP:O:-0.76511");
		params.push_back("fa_standard:TYR:1HB:0.07353");
		params.push_back("fa_standard:TYR:2HB:0.07353");
		params.push_back("fa_standard:TYR:C:0.76511");
		params.push_back("fa_standard:TYR:CA:0.11475");
		params.push_back("fa_standard:TYR:CB:-0.14706");
		params.push_back("fa_standard:TYR:CD1:-0.09395");
		params.push_back("fa_standard:TYR:CD2:-0.09395");
		params.push_back("fa_standard:TYR:CE1:-0.09395");
		params.push_back("fa_standard:TYR:CE2:-0.09395");
		params.push_back("fa_standard:TYR:CG:0.00000");
		params.push_back("fa_standard:TYR:CZ:0.08987");
		params.push_back("fa_standard:TYR:H:0.50819");
		params.push_back("fa_standard:TYR:HA:0.14754");
		params.push_back("fa_standard:TYR:HD1:0.09395");
		params.push_back("fa_standard:TYR:HD2:0.09395");
		params.push_back("fa_standard:TYR:HE1:0.09395");
		params.push_back("fa_standard:TYR:HE2:0.09395");
		params.push_back("fa_standard:TYR:HH:0.35130");
		params.push_back("fa_standard:TYR:N:-0.77049");
		params.push_back("fa_standard:TYR:O:-0.76511");
		params.push_back("fa_standard:TYR:OH:-0.44117");
		params.push_back("fa_standard:VAL:1HG1:0.12155");
		params.push_back("fa_standard:VAL:1HG2:0.12155");
		params.push_back("fa_standard:VAL:2HG1:0.12155");
		params.push_back("fa_standard:VAL:2HG2:0.12155");
		params.push_back("fa_standard:VAL:3HG1:0.12155");
		params.push_back("fa_standard:VAL:3HG2:0.12155");
		params.push_back("fa_standard:VAL:C:0.76511");
		params.push_back("fa_standard:VAL:CA:0.11475");
		params.push_back("fa_standard:VAL:CB:-0.12155");
		params.push_back("fa_standard:VAL:CG1:-0.36465");
		params.push_back("fa_standard:VAL:CG2:-0.36465");
		params.push_back("fa_standard:VAL:H:0.50819");
		params.push_back("fa_standard:VAL:HA:0.14754");
		params.push_back("fa_standard:VAL:HB:0.12155");
		params.push_back("fa_standard:VAL:N:-0.77049");
		params.push_back("fa_standard:VAL:O:-0.76511");
		params.push_back("fa_standard:DAL:1HB:0.12155");
		params.push_back("fa_standard:DAL:2HB:0.12155");
		params.push_back("fa_standard:DAL:3HB:0.12155");
		params.push_back("fa_standard:DAL:C:0.76511");
		params.push_back("fa_standard:DAL:CA:0.11475");
		params.push_back("fa_standard:DAL:CB:-0.36465");
		params.push_back("fa_standard:DAL:H:0.50819");
		params.push_back("fa_standard:DAL:HA:0.14754");
		params.push_back("fa_standard:DAL:N:-0.77049");
		params.push_back("fa_standard:DAL:O:-0.76511");
		params.push_back("fa_standard:DAR:1HB:0.07619");
		params.push_back("fa_standard:DAR:1HD:0.07619");
		params.push_back("fa_standard:DAR:1HG:0.07619");
		params.push_back("fa_standard:DAR:1HH1:0.29787");
		params.push_back("fa_standard:DAR:1HH2:0.29787");
		params.push_back("fa_standard:DAR:2HB:0.07619");
		params.push_back("fa_standard:DAR:2HD:0.07619");
		params.push_back("fa_standard:DAR:2HG:0.07619");
		params.push_back("fa_standard:DAR:2HH1:0.29787");
		params.push_back("fa_standard:DAR:2HH2:0.29787");
		params.push_back("fa_standard:DAR:C:0.76511");
		params.push_back("fa_standard:DAR:CA:0.11475");
		params.push_back("fa_standard:DAR:CB:-0.08557");
		params.push_back("fa_standard:DAR:CD:0.14210");
		params.push_back("fa_standard:DAR:CG:-0.08557");
		params.push_back("fa_standard:DAR:CZ:0.40572");
		params.push_back("fa_standard:DAR:H:0.50819");
		params.push_back("fa_standard:DAR:HA:0.14754");
		params.push_back("fa_standard:DAR:HE:0.28589");
		params.push_back("fa_standard:DAR:N:-0.77049");
		params.push_back("fa_standard:DAR:NE:-0.39713");
		params.push_back("fa_standard:DAR:NH1:-0.45704");
		params.push_back("fa_standard:DAR:NH2:-0.45704");
		params.push_back("fa_standard:DAR:O:-0.76511");
		params.push_back("fa_standard:DAN:1HB:0.07322");
		params.push_back("fa_standard:DAN:1HD2:0.26032");
		params.push_back("fa_standard:DAN:2HB:0.07322");
		params.push_back("fa_standard:DAN:2HD2:0.24405");
		params.push_back("fa_standard:DAN:C:0.76511");
		params.push_back("fa_standard:DAN:CA:0.11475");
		params.push_back("fa_standard:DAN:CB:-0.14643");
		params.push_back("fa_standard:DAN:CG:0.44743");
		params.push_back("fa_standard:DAN:H:0.50819");
		params.push_back("fa_standard:DAN:HA:0.14754");
		params.push_back("fa_standard:DAN:N:-0.77049");
		params.push_back("fa_standard:DAN:ND2:-0.50437");
		params.push_back("fa_standard:DAN:O:-0.76511");
		params.push_back("fa_standard:DAN:OD1:-0.44743");
		params.push_back("fa_standard:DAS:1HB:0.13276");
		params.push_back("fa_standard:DAS:2HB:0.13276");
		params.push_back("fa_standard:DAS:C:0.76511");
		params.push_back("fa_standard:DAS:CA:0.11475");
		params.push_back("fa_standard:DAS:CB:-0.29888");
		params.push_back("fa_standard:DAS:CG:0.75106");
		params.push_back("fa_standard:DAS:H:0.50819");
		params.push_back("fa_standard:DAS:HA:0.14754");
		params.push_back("fa_standard:DAS:N:-0.77049");
		params.push_back("fa_standard:DAS:O:-0.76511");
		params.push_back("fa_standard:DAS:OD1:-0.85885");
		params.push_back("fa_standard:DAS:OD2:-0.85885");
		params.push_back("fa_standard:DCS:1HB:0.12155");
		params.push_back("fa_standard:DCS:2HB:0.12155");
		params.push_back("fa_standard:DCS:C:0.76511");
		params.push_back("fa_standard:DCS:CA:0.11475");
		params.push_back("fa_standard:DCS:CB:-0.14856");
		params.push_back("fa_standard:DCS:H:0.50819");
		params.push_back("fa_standard:DCS:HA:0.14754");
		params.push_back("fa_standard:DCS:HG:0.21609");
		params.push_back("fa_standard:DCS:N:-0.77049");
		params.push_back("fa_standard:DCS:O:-0.76511");
		params.push_back("fa_standard:DCS:SG:-0.31063");
		params.push_back("fa_standard:DGN:1HB:0.09669");
		params.push_back("fa_standard:DGN:1HE2:0.34378");
		params.push_back("fa_standard:DGN:1HG:0.09669");
		params.push_back("fa_standard:DGN:2HB:0.09669");
		params.push_back("fa_standard:DGN:2HE2:0.32229");
		params.push_back("fa_standard:DGN:2HG:0.09669");
		params.push_back("fa_standard:DGN:C:0.76511");
		params.push_back("fa_standard:DGN:CA:0.11475");
		params.push_back("fa_standard:DGN:CB:-0.19338");
		params.push_back("fa_standard:DGN:CD:0.59087");
		params.push_back("fa_standard:DGN:CG:-0.19338");
		params.push_back("fa_standard:DGN:H:0.50819");
		params.push_back("fa_standard:DGN:HA:0.14754");
		params.push_back("fa_standard:DGN:N:-0.77049");
		params.push_back("fa_standard:DGN:NE2:-0.66607");
		params.push_back("fa_standard:DGN:O:-0.76511");
		params.push_back("fa_standard:DGN:OE1:-0.59087");
		params.push_back("fa_standard:DGU:1HB:0.05862");
		params.push_back("fa_standard:DGU:1HG:0.05862");
		params.push_back("fa_standard:DGU:2HB:0.05862");
		params.push_back("fa_standard:DGU:2HG:0.05862");
		params.push_back("fa_standard:DGU:C:0.76511");
		params.push_back("fa_standard:DGU:CA:0.11475");
		params.push_back("fa_standard:DGU:CB:-0.16925");
		params.push_back("fa_standard:DGU:CD:0.50590");
		params.push_back("fa_standard:DGU:CG:-0.25364");
		params.push_back("fa_standard:DGU:H:0.50819");
		params.push_back("fa_standard:DGU:HA:0.14754");
		params.push_back("fa_standard:DGU:N:-0.77049");
		params.push_back("fa_standard:DGU:O:-0.76511");
		params.push_back("fa_standard:DGU:OE1:-0.65874");
		params.push_back("fa_standard:DGU:OE2:-0.65874");
		params.push_back("fa_standard:DHI_D:1HB:0.08456");
		params.push_back("fa_standard:DHI_D:2HB:0.08456");
		params.push_back("fa_standard:DHI_D:C:0.76511");
		params.push_back("fa_standard:DHI_D:CA:0.11475");
		params.push_back("fa_standard:DHI_D:CB:-0.08456");
		params.push_back("fa_standard:DHI_D:CD2:0.20671");
		params.push_back("fa_standard:DHI_D:CE1:0.23490");
		params.push_back("fa_standard:DHI_D:CG:-0.04698");
		params.push_back("fa_standard:DHI_D:H:0.50819");
		params.push_back("fa_standard:DHI_D:HA:0.14754");
		params.push_back("fa_standard:DHI_D:HD1:0.30067");
		params.push_back("fa_standard:DHI_D:HD2:0.09396");
		params.push_back("fa_standard:DHI_D:HE1:0.12215");
		params.push_back("fa_standard:DHI_D:N:-0.77049");
		params.push_back("fa_standard:DHI_D:ND1:-0.33825");
		params.push_back("fa_standard:DHI_D:NE2:-0.65772");
		params.push_back("fa_standard:DHI_D:O:-0.76511");
		params.push_back("fa_standard:DHI:1HB:0.08456");
		params.push_back("fa_standard:DHI:2HB:0.08456");
		params.push_back("fa_standard:DHI:C:0.76511");
		params.push_back("fa_standard:DHI:CA:0.11475");
		params.push_back("fa_standard:DHI:CB:-0.07517");
		params.push_back("fa_standard:DHI:CD2:-0.04698");
		params.push_back("fa_standard:DHI:CE1:0.23490");
		params.push_back("fa_standard:DHI:CG:0.20671");
		params.push_back("fa_standard:DHI:H:0.50819");
		params.push_back("fa_standard:DHI:HA:0.14754");
		params.push_back("fa_standard:DHI:HD2:0.08456");
		params.push_back("fa_standard:DHI:HE1:0.12215");
		params.push_back("fa_standard:DHI:HE2:0.30067");
		params.push_back("fa_standard:DHI:N:-0.77049");
		params.push_back("fa_standard:DHI:ND1:-0.65772");
		params.push_back("fa_standard:DHI:NE2:-0.33825");
		params.push_back("fa_standard:DHI:O:-0.76511");
		params.push_back("fa_standard:DIL:1HD1:0.12155");
		params.push_back("fa_standard:DIL:1HG1:0.12155");
		params.push_back("fa_standard:DIL:1HG2:0.12155");
		params.push_back("fa_standard:DIL:2HD1:0.12155");
		params.push_back("fa_standard:DIL:2HG1:0.12155");
		params.push_back("fa_standard:DIL:2HG2:0.12155");
		params.push_back("fa_standard:DIL:3HD1:0.12155");
		params.push_back("fa_standard:DIL:3HG2:0.12155");
		params.push_back("fa_standard:DIL:C:0.76511");
		params.push_back("fa_standard:DIL:CA:0.11475");
		params.push_back("fa_standard:DIL:CB:-0.12155");
		params.push_back("fa_standard:DIL:CD1:-0.36465");
		params.push_back("fa_standard:DIL:CG1:-0.24310");
		params.push_back("fa_standard:DIL:CG2:-0.36465");
		params.push_back("fa_standard:DIL:H:0.50819");
		params.push_back("fa_standard:DIL:HA:0.14754");
		params.push_back("fa_standard:DIL:HB:0.12155");
		params.push_back("fa_standard:DIL:N:-0.77049");
		params.push_back("fa_standard:DIL:O:-0.76511");
		params.push_back("fa_standard:DLE:1HB:0.12155");
		params.push_back("fa_standard:DLE:1HD1:0.12155");
		params.push_back("fa_standard:DLE:1HD2:0.12155");
		params.push_back("fa_standard:DLE:2HB:0.12155");
		params.push_back("fa_standard:DLE:2HD1:0.12155");
		params.push_back("fa_standard:DLE:2HD2:0.12155");
		params.push_back("fa_standard:DLE:3HD1:0.12155");
		params.push_back("fa_standard:DLE:3HD2:0.12155");
		params.push_back("fa_standard:DLE:C:0.76511");
		params.push_back("fa_standard:DLE:CA:0.11475");
		params.push_back("fa_standard:DLE:CB:-0.24310");
		params.push_back("fa_standard:DLE:CD1:-0.36465");
		params.push_back("fa_standard:DLE:CD2:-0.36465");
		params.push_back("fa_standard:DLE:CG:-0.12155");
		params.push_back("fa_standard:DLE:H:0.50819");
		params.push_back("fa_standard:DLE:HA:0.14754");
		params.push_back("fa_standard:DLE:HG:0.12155");
		params.push_back("fa_standard:DLE:N:-0.77049");
		params.push_back("fa_standard:DLE:O:-0.76511");
		params.push_back("fa_standard:DLY:1HB:0.08943");
		params.push_back("fa_standard:DLY:1HD:0.08943");
		params.push_back("fa_standard:DLY:1HE:0.05026");
		params.push_back("fa_standard:DLY:1HG:0.08943");
		params.push_back("fa_standard:DLY:1HZ:0.32446");
		params.push_back("fa_standard:DLY:2HB:0.08943");
		params.push_back("fa_standard:DLY:2HD:0.08943");
		params.push_back("fa_standard:DLY:2HE:0.05026");
		params.push_back("fa_standard:DLY:2HG:0.08943");
		params.push_back("fa_standard:DLY:2HZ:0.32446");
		params.push_back("fa_standard:DLY:3HZ:0.32446");
		params.push_back("fa_standard:DLY:C:0.76511");
		params.push_back("fa_standard:DLY:CA:0.11475");
		params.push_back("fa_standard:DLY:CB:-0.17497");
		params.push_back("fa_standard:DLY:CD:-0.17497");
		params.push_back("fa_standard:DLY:CE:0.20694");
		params.push_back("fa_standard:DLY:CG:-0.17497");
		params.push_back("fa_standard:DLY:H:0.50819");
		params.push_back("fa_standard:DLY:HA:0.14754");
		params.push_back("fa_standard:DLY:N:-0.77049");
		params.push_back("fa_standard:DLY:NZ:-0.29249");
		params.push_back("fa_standard:DLY:O:-0.76511");
		params.push_back("fa_standard:DME:1HB:0.12155");
		params.push_back("fa_standard:DME:1HE:0.12155");
		params.push_back("fa_standard:DME:1HG:0.12155");
		params.push_back("fa_standard:DME:2HB:0.12155");
		params.push_back("fa_standard:DME:2HE:0.12155");
		params.push_back("fa_standard:DME:2HG:0.12155");
		params.push_back("fa_standard:DME:3HE:0.12155");
		params.push_back("fa_standard:DME:C:0.76511");
		params.push_back("fa_standard:DME:CA:0.11475");
		params.push_back("fa_standard:DME:CB:-0.24310");
		params.push_back("fa_standard:DME:CE:-0.29712");
		params.push_back("fa_standard:DME:CG:-0.18908");
		params.push_back("fa_standard:DME:H:0.50819");
		params.push_back("fa_standard:DME:HA:0.14754");
		params.push_back("fa_standard:DME:N:-0.77049");
		params.push_back("fa_standard:DME:O:-0.76511");
		params.push_back("fa_standard:DME:SD:-0.12155");
		params.push_back("fa_standard:DPH:1HB:0.14323");
		params.push_back("fa_standard:DPH:2HB:0.14323");
		params.push_back("fa_standard:DPH:C:0.76511");
		params.push_back("fa_standard:DPH:CA:0.11475");
		params.push_back("fa_standard:DPH:CB:-0.28647");
		params.push_back("fa_standard:DPH:CD1:-0.18302");
		params.push_back("fa_standard:DPH:CD2:-0.18302");
		params.push_back("fa_standard:DPH:CE1:-0.18302");
		params.push_back("fa_standard:DPH:CE2:-0.18302");
		params.push_back("fa_standard:DPH:CG:0.00000");
		params.push_back("fa_standard:DPH:CZ:-0.18302");
		params.push_back("fa_standard:DPH:H:0.50819");
		params.push_back("fa_standard:DPH:HA:0.14754");
		params.push_back("fa_standard:DPH:HD1:0.18302");
		params.push_back("fa_standard:DPH:HD2:0.18302");
		params.push_back("fa_standard:DPH:HE1:0.18302");
		params.push_back("fa_standard:DPH:HE2:0.18302");
		params.push_back("fa_standard:DPH:HZ:0.18302");
		params.push_back("fa_standard:DPH:N:-0.77049");
		params.push_back("fa_standard:DPH:O:-0.76511");
		params.push_back("fa_standard:DPR:1HB:0.14754");
		params.push_back("fa_standard:DPR:1HD:0.14754");
		params.push_back("fa_standard:DPR:1HG:0.14754");
		params.push_back("fa_standard:DPR:2HB:0.14754");
		params.push_back("fa_standard:DPR:2HD:0.14754");
		params.push_back("fa_standard:DPR:2HG:0.14754");
		params.push_back("fa_standard:DPR:C:0.76511");
		params.push_back("fa_standard:DPR:CA:0.03279");
		params.push_back("fa_standard:DPR:CB:-0.29508");
		params.push_back("fa_standard:DPR:CG:-0.29508");
		params.push_back("fa_standard:DPR:HA:0.14754");
		params.push_back("fa_standard:DPR:N:-0.47541");
		params.push_back("fa_standard:DPR:O:-0.76511");
		params.push_back("fa_standard:DSE:1HB:0.07337");
		params.push_back("fa_standard:DSE:2HB:0.07337");
		params.push_back("fa_standard:DSE:C:0.76511");
		params.push_back("fa_standard:DSE:CA:0.11475");
		params.push_back("fa_standard:DSE:CB:0.04076");
		params.push_back("fa_standard:DSE:H:0.50819");
		params.push_back("fa_standard:DSE:HA:0.14754");
		params.push_back("fa_standard:DSE:HG:0.35053");
		params.push_back("fa_standard:DSE:N:-0.77049");
		params.push_back("fa_standard:DSE:O:-0.76511");
		params.push_back("fa_standard:DSE:OG:-0.53802");
		params.push_back("fa_standard:DTH:1HG2:0.08381");
		params.push_back("fa_standard:DTH:2HG2:0.08381");
		params.push_back("fa_standard:DTH:3HG2:0.08381");
		params.push_back("fa_standard:DTH:C:0.76511");
		params.push_back("fa_standard:DTH:CA:0.11475");
		params.push_back("fa_standard:DTH:CB:0.13037");
		params.push_back("fa_standard:DTH:CG2:-0.25143");
		params.push_back("fa_standard:DTH:H:0.50819");
		params.push_back("fa_standard:DTH:HA:0.14754");
		params.push_back("fa_standard:DTH:HB:0.08381");
		params.push_back("fa_standard:DTH:HG1:0.40042");
		params.push_back("fa_standard:DTH:N:-0.77049");
		params.push_back("fa_standard:DTH:O:-0.76511");
		params.push_back("fa_standard:DTH:OG1:-0.61459");
		params.push_back("fa_standard:DTR:1HB:0.10831");
		params.push_back("fa_standard:DTR:2HB:0.10831");
		params.push_back("fa_standard:DTR:C:0.76511");
		params.push_back("fa_standard:DTR:CA:0.11475");
		params.push_back("fa_standard:DTR:CB:-0.21662");
		params.push_back("fa_standard:DTR:CD1:0.04212");
		params.push_back("fa_standard:DTR:CD2:-0.02407");
		params.push_back("fa_standard:DTR:CE2:0.15644");
		params.push_back("fa_standard:DTR:CE3:-0.13839");
		params.push_back("fa_standard:DTR:CG:-0.03610");
		params.push_back("fa_standard:DTR:CH2:-0.13839");
		params.push_back("fa_standard:DTR:CZ2:-0.13839");
		params.push_back("fa_standard:DTR:CZ3:-0.13839");
		params.push_back("fa_standard:DTR:H:0.50819");
		params.push_back("fa_standard:DTR:HA:0.14754");
		params.push_back("fa_standard:DTR:HD1:0.13839");
		params.push_back("fa_standard:DTR:HE1:0.45730");
		params.push_back("fa_standard:DTR:HE3:0.13839");
		params.push_back("fa_standard:DTR:HH2:0.13839");
		params.push_back("fa_standard:DTR:HZ2:0.13839");
		params.push_back("fa_standard:DTR:HZ3:0.13839");
		params.push_back("fa_standard:DTR:N:-0.77049");
		params.push_back("fa_standard:DTR:NE1:-0.73409");
		params.push_back("fa_standard:DTR:O:-0.76511");
		params.push_back("fa_standard:DTY:1HB:0.07353");
		params.push_back("fa_standard:DTY:2HB:0.07353");
		params.push_back("fa_standard:DTY:C:0.76511");
		params.push_back("fa_standard:DTY:CA:0.11475");
		params.push_back("fa_standard:DTY:CB:-0.14706");
		params.push_back("fa_standard:DTY:CD1:-0.09395");
		params.push_back("fa_standard:DTY:CD2:-0.09395");
		params.push_back("fa_standard:DTY:CE1:-0.09395");
		params.push_back("fa_standard:DTY:CE2:-0.09395");
		params.push_back("fa_standard:DTY:CG:0.00000");
		params.push_back("fa_standard:DTY:CZ:0.08987");
		params.push_back("fa_standard:DTY:H:0.50819");
		params.push_back("fa_standard:DTY:HA:0.14754");
		params.push_back("fa_standard:DTY:HD1:0.09395");
		params.push_back("fa_standard:DTY:HD2:0.09395");
		params.push_back("fa_standard:DTY:HE1:0.09395");
		params.push_back("fa_standard:DTY:HE2:0.09395");
		params.push_back("fa_standard:DTY:HH:0.35130");
		params.push_back("fa_standard:DTY:N:-0.77049");
		params.push_back("fa_standard:DTY:O:-0.76511");
		params.push_back("fa_standard:DTY:OH:-0.44117");
		params.push_back("fa_standard:DVA:1HG1:0.12155");
		params.push_back("fa_standard:DVA:1HG2:0.12155");
		params.push_back("fa_standard:DVA:2HG1:0.12155");
		params.push_back("fa_standard:DVA:2HG2:0.12155");
		params.push_back("fa_standard:DVA:3HG1:0.12155");
		params.push_back("fa_standard:DVA:3HG2:0.12155");
		params.push_back("fa_standard:DVA:C:0.76511");
		params.push_back("fa_standard:DVA:CA:0.11475");
		params.push_back("fa_standard:DVA:CB:-0.12155");
		params.push_back("fa_standard:DVA:CG1:-0.36465");
		params.push_back("fa_standard:DVA:CG2:-0.36465");
		params.push_back("fa_standard:DVA:H:0.50819");
		params.push_back("fa_standard:DVA:HA:0.14754");
		params.push_back("fa_standard:DVA:HB:0.12155");
		params.push_back("fa_standard:DVA:N:-0.77049");
		params.push_back("fa_standard:DVA:O:-0.76511");

		// fd & rpav water
		params.push_back("fa_standard:HOH:O:-0.2896");
		params.push_back("fa_standard:HOH:H1:0.1448");
		params.push_back("fa_standard:HOH:H2:0.1448");

		options[ basic::options::OptionKeys::chemical::set_atomic_charge ].value(params);
	} else {
		TR.Warning << "Flag -beta_jan25 is set but -set_atomic_charge is also specified.  Not changing atom properties!" << std::endl;
	}

	// patch charges
	if ( ! options[ basic::options::OptionKeys::chemical::set_patch_atomic_charge ].user() ) {
		utility::vector1< std::string > params;
		params.push_back("fa_standard:ALA:CtermProteinFull:C:0.35448");
		params.push_back("fa_standard:ALA:CtermProteinFull:O:-0.67724");
		params.push_back("fa_standard:ALA:CtermProteinFull:OXT:-0.67724");
		params.push_back("fa_standard:ALA:NtermProteinFull:1H:0.26154");
		params.push_back("fa_standard:ALA:NtermProteinFull:2H:0.26154");
		params.push_back("fa_standard:ALA:NtermProteinFull:3H:0.26154");
		params.push_back("fa_standard:ALA:NtermProteinFull:CA:0.19184");
		params.push_back("fa_standard:ALA:NtermProteinFull:HA:0.12794");
		params.push_back("fa_standard:ALA:NtermProteinFull:N:-0.10440");
		params.push_back("fa_standard:DAL:CtermProteinFull:C:0.35448");
		params.push_back("fa_standard:DAL:CtermProteinFull:O:-0.67724");
		params.push_back("fa_standard:DAL:CtermProteinFull:OXT:-0.67724");
		params.push_back("fa_standard:DAL:NtermProteinFull:1H:0.26154");
		params.push_back("fa_standard:DAL:NtermProteinFull:2H:0.26154");
		params.push_back("fa_standard:DAL:NtermProteinFull:3H:0.26154");
		params.push_back("fa_standard:DAL:NtermProteinFull:CA:0.19184");
		params.push_back("fa_standard:DAL:NtermProteinFull:HA:0.12794");
		params.push_back("fa_standard:DAL:NtermProteinFull:N:-0.10440");
		params.push_back("fa_standard:GLY:NtermProteinFull:1H:0.25156");
		params.push_back("fa_standard:GLY:NtermProteinFull:1HA:0.11215");
		params.push_back("fa_standard:GLY:NtermProteinFull:2H:0.25156");
		params.push_back("fa_standard:GLY:NtermProteinFull:2HA:0.11215");
		params.push_back("fa_standard:GLY:NtermProteinFull:3H:0.25156");
		params.push_back("fa_standard:GLY:NtermProteinFull:CA:0.13539");
		params.push_back("fa_standard:GLY:NtermProteinFull:N:-0.11438");
		params.push_back("fa_standard:PRO:NtermProteinFull:1H:0.17880");
		params.push_back("fa_standard:PRO:NtermProteinFull:1HB:0.09167");
		params.push_back("fa_standard:PRO:NtermProteinFull:1HD:0.09167");
		params.push_back("fa_standard:PRO:NtermProteinFull:1HG:0.09167");
		params.push_back("fa_standard:PRO:NtermProteinFull:2H:0.17880");
		params.push_back("fa_standard:PRO:NtermProteinFull:2HB:0.09167");
		params.push_back("fa_standard:PRO:NtermProteinFull:2HD:0.09167");
		params.push_back("fa_standard:PRO:NtermProteinFull:2HG:0.09167");
		params.push_back("fa_standard:PRO:NtermProteinFull:CA:0.13233");
		params.push_back("fa_standard:PRO:NtermProteinFull:CB:-0.06516");
		params.push_back("fa_standard:PRO:NtermProteinFull:CG:-0.06516");
		params.push_back("fa_standard:PRO:NtermProteinFull:HA:0.09167");
		params.push_back("fa_standard:PRO:NtermProteinFull:N:-0.00127");
		params.push_back("fa_standard:DPR:NtermProteinFull:1H:0.17880");
		params.push_back("fa_standard:DPR:NtermProteinFull:1HB:0.09167");
		params.push_back("fa_standard:DPR:NtermProteinFull:1HD:0.09167");
		params.push_back("fa_standard:DPR:NtermProteinFull:1HG:0.09167");
		params.push_back("fa_standard:DPR:NtermProteinFull:2H:0.17880");
		params.push_back("fa_standard:DPR:NtermProteinFull:2HB:0.09167");
		params.push_back("fa_standard:DPR:NtermProteinFull:2HD:0.09167");
		params.push_back("fa_standard:DPR:NtermProteinFull:2HG:0.09167");
		params.push_back("fa_standard:DPR:NtermProteinFull:CA:0.13233");
		params.push_back("fa_standard:DPR:NtermProteinFull:CB:-0.06516");
		params.push_back("fa_standard:DPR:NtermProteinFull:CG:-0.06516");
		params.push_back("fa_standard:DPR:NtermProteinFull:HA:0.09167");
		params.push_back("fa_standard:DPR:NtermProteinFull:N:-0.00127");
		options[ basic::options::OptionKeys::chemical::set_patch_atomic_charge ].value(params);
	} else {
		TR.Warning << "Flag -beta_jan25 is set but -set_patch_atomic_charge is also specified.  Not changing atom properties!" << std::endl;
	}

	// hbonds
	if ( ! options[ basic::options::OptionKeys::score::hb_acc_strength ].user() ) {
		utility::vector1< std::string > params;
		params.push_back( "hbacc_AHX:1.15");
		params.push_back( "hbacc_CXA:1.21");
		params.push_back( "hbacc_CXL:1.10");
		params.push_back( "hbacc_HXL:1.15");
		params.push_back( "hbacc_IMD:1.13");
		params.push_back( "hbacc_IME:1.17");
		params.push_back( "hbacc_PBA:1.19");
		params.push_back( "hbacc_H2O:1.23");
		params.push_back( "hbacc_GENERIC_SP2SC:1.19");
		params.push_back( "hbacc_GENERIC_SP3SC:1.19");
		params.push_back( "hbacc_GENERIC_RINGSC:1.19");
		options[ basic::options::OptionKeys::score::hb_acc_strength ].value(params);
	} else {
		TR.Warning << "Flag -beta_jan25 is set but -hb_acc_strength are also specified.  Not changing atom properties!" << std::endl;
	}
	if ( ! options[ basic::options::OptionKeys::score::hb_don_strength ].user() ) {
		utility::vector1< std::string > params;
		params.push_back( "hbdon_AHX:1.00");
		params.push_back( "hbdon_AMO:1.17");
		params.push_back( "hbdon_CXA:1.29");
		params.push_back( "hbdon_GDE:1.11");
		params.push_back( "hbdon_GDH:1.11");
		params.push_back( "hbdon_HXL:0.99");
		params.push_back( "hbdon_IMD:1.18");
		params.push_back( "hbdon_IME:1.42");
		params.push_back( "hbdon_IND:1.15");
		params.push_back( "hbdon_PBA:1.45");
		params.push_back( "hbdon_H2O:1.26");
		params.push_back( "hbdon_GENERIC_SC:1.45");
		options[ basic::options::OptionKeys::score::hb_don_strength ].value(params);
	} else {
		TR.Warning << "Flag -beta_jan25 is set but -hb_don_strength are also specified.  Not changing atom properties!" << std::endl;
	}

	// unmodifypot
	if ( ! options[ basic::options::OptionKeys::score::unmodifypot ].user() ) {
		options[ basic::options::OptionKeys::score::unmodifypot ].value(true);
	}

	// attractive hydrogens
	if ( ! options[ basic::options::OptionKeys::score::fa_Hatr ].user() ) {
		options[ basic::options::OptionKeys::score::fa_Hatr ].value(true);
	}

	// new sp3 acceptor
	if ( ! options[ basic::options::OptionKeys::score::hbond_new_sp3_acc ].user() ) {
		options[ basic::options::OptionKeys::score::hbond_new_sp3_acc ].value(true);
	}

	// sigmoidal dielectric
	if ( ! options[ basic::options::OptionKeys::score::elec_sigmoidal_die ].user() ) {
		options[ basic::options::OptionKeys::score::elec_sigmoidal_die ].value(true);
	}
	if ( ! options[ basic::options::OptionKeys::score::elec_sigmoidal_die_D ].user() ) {
		options[ basic::options::OptionKeys::score::elec_sigmoidal_die_D ].value(79.931);
	}
	if ( ! options[ basic::options::OptionKeys::score::elec_sigmoidal_die_D0 ].user() ) {
		options[ basic::options::OptionKeys::score::elec_sigmoidal_die_D0 ].value(6.648);
	}
	if ( ! options[ basic::options::OptionKeys::score::elec_sigmoidal_die_S ].user() ) {
		options[ basic::options::OptionKeys::score::elec_sigmoidal_die_S ].value(0.441546);
	}

	// lkball
	if ( ! options[ basic::options::OptionKeys::dna::specificity::lk_ball_ramp_width_A2 ].user() ) {
		options[ basic::options::OptionKeys::dna::specificity::lk_ball_ramp_width_A2 ].value(3.709);
	}
	if ( ! options[ basic::options::OptionKeys::dna::specificity::lk_ball_overlap_width_A2 ].user() ) {
		options[ basic::options::OptionKeys::dna::specificity::lk_ball_overlap_width_A2 ].value(2.60);
	}
	if ( ! options[ basic::options::OptionKeys::dna::specificity::lk_ball_overlap_gap ].user() ) {
		options[ basic::options::OptionKeys::dna::specificity::lk_ball_overlap_gap ].value(0.50);
	}
	if ( ! options[ basic::options::OptionKeys::dna::specificity::lk_ball_water_fade ].user() ) {
		options[ basic::options::OptionKeys::dna::specificity::lk_ball_water_fade ].value(1.0);
	}
	if ( ! options[ basic::options::OptionKeys::dna::specificity::lk_ball_for_bb ].user() ) {
		options[ basic::options::OptionKeys::dna::specificity::lk_ball_for_bb ].value(true);
	}
	if ( ! options[ basic::options::OptionKeys::dna::specificity::lk_ball_bridge_angle_widthscale ].user() ) {
		options[ basic::options::OptionKeys::dna::specificity::lk_ball_bridge_angle_widthscale ].value(2.8);
	}
	if ( ! options[ basic::options::OptionKeys::dna::specificity::lk_ball_waters_sp2 ].user() ) {
		utility::vector1< core::Real > params;
		params.push_back(2.828);
		params.push_back(134.5);
		params.push_back(0.0);
		params.push_back(2.828);
		params.push_back(134.5);
		params.push_back(180.0);
		options[ basic::options::OptionKeys::dna::specificity::lk_ball_waters_sp2 ].value(params);
	}
	if ( ! options[ basic::options::OptionKeys::dna::specificity::lk_ball_waters_sp3 ].user() ) {
		utility::vector1< core::Real > params;
		params.push_back(2.828);
		params.push_back(109.3);
		params.push_back(120.0);
		params.push_back(2.828);
		params.push_back(109.3);
		params.push_back(240.0);
		options[ basic::options::OptionKeys::dna::specificity::lk_ball_waters_sp3 ].value(params);
	}
	if ( ! options[ basic::options::OptionKeys::dna::specificity::lk_ball_waters_ring ].user() ) {
		utility::vector1< core::Real > params;
		params.push_back(2.828);
		params.push_back(180);
		params.push_back(0);
		options[ basic::options::OptionKeys::dna::specificity::lk_ball_waters_ring ].value(params);
	}
	if ( ! options[ basic::options::OptionKeys::dna::specificity::lk_ball_waters_donor ].user() ) {
		options[ basic::options::OptionKeys::dna::specificity::lk_ball_waters_donor ].value(2.828);
	}


	// roland's new smoothed p_aa_pp and fa_dun libraries
	if ( ! options[ basic::options::OptionKeys::corrections::shapovalov_lib_fixes_enable ].user() ) {
		options[ basic::options::OptionKeys::corrections::shapovalov_lib_fixes_enable ].value(true);
	}
	if ( ! options[ basic::options::OptionKeys::corrections::shapovalov_lib::shap_dun10_enable ].user() ) {
		options[ basic::options::OptionKeys::corrections::shapovalov_lib::shap_dun10_enable ].value(false);
	}
	if ( ! options[ basic::options::OptionKeys::corrections::shapovalov_lib::shap_p_aa_pp_smooth_level ].user() ) {
		options[ basic::options::OptionKeys::corrections::shapovalov_lib::shap_p_aa_pp_smooth_level ].value("1");
	}

	// corrected torsions
	if ( ! options[ basic::options::OptionKeys::corrections::score::rama_pp_map ].user() ) {
		options[ basic::options::OptionKeys::corrections::score::rama_pp_map ].value("scoring/score_functions/rama/fd_beta_nov2016");
	}
	if ( ! options[ basic::options::OptionKeys::corrections::score::dun10_dir ].user() ) {
		options[ basic::options::OptionKeys::corrections::score::dun10_dir ].value("rotamer/beta_nov2016");
	}


	// bbdep omega
	if ( ! options[ basic::options::OptionKeys::corrections::score::bbdep_omega ].user() ) {
		options[ basic::options::OptionKeys::corrections::score::bbdep_omega ].value(true);
	}

	// updated fa_elec countpair logic
	if ( ! options[ basic::options::OptionKeys::score::elec_representative_cp_flip ].user() ) {
		options[ basic::options::OptionKeys::score::elec_representative_cp_flip ].value(true);
	}

	// updated cartbonded parameters
	if ( ! options[ basic::options::OptionKeys::score::bonded_params_dir ].user() ) {
		options[ basic::options::OptionKeys::score::bonded_params_dir ].value("scoring/score_functions/bondlength_bondangle");
	}

	// weights file
	if ( setweight && !options[ score::weights ].user() ) {
		if ( options[corrections::beta_jan25_cart]() || options[optimization::nonideal]() ) {
			options[ score::weights ].value( "beta_jan25_cart.wts" );
		} else {
			options[ score::weights ].value( "beta_jan25.wts" );
		}
	} else {
		TR.Warning << "Flag -beta_jan25 is set but either -weights are also specified or this is being used with the genpot score function.  Not changing input weights file!" << std::endl;
	}

	// protocol option but it should be default ...
	if ( !options[optimization::default_max_cycles]() ) {
		options[ optimization::default_max_cycles ].value( 200 );
	}
}

void
init_nonideal_correction( utility::options::OptionCollection & options ) {
	//PTC set up for flexible bond geometries
	if ( options[optimization::nonideal] ) {
		// add jumps to missing densities so cart_bonded_length behaves correctly
		if ( ! options[ in::missing_density_to_jump ].user() ) {
			options[ in::missing_density_to_jump ].value( true );
		}

		// use dualspace relax protocol - start with internal coordinate, finish with Cartesian
		if ( ! options[ relax::dualspace ].user() ) {
			options[ relax::dualspace ].value( true );
		}

		// use flexible bond angles during internal coordinate stage of dualspace relax
		if ( ! options[ relax::minimize_bond_angles ].user() ) {
			options[ relax::minimize_bond_angles ].value( true );
		}

		// limit default_max_cycles
		if ( ! options[ optimization::default_max_cycles ].user() ) {
			options[ optimization::default_max_cycles ].value( 200 );
		}

		// enable for rtmin/minpack - use nonideal unless nonideal option is defined or cartesian is specified as true
		if ( (! options[ optimization::scmin_nonideal ].user()) || options[ optimization::scmin_cartesian ] ) {
			options[ optimization::scmin_nonideal ].value( true );
		}

		// talaris2013_cart.wts - enables cart_bonded, disables pro_close, new reference weights
		if ( ! options[ score::weights ].user() ) {
			options[ score::weights ].value( "talaris2013_cart.wts" );
		} else {
			scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
			if ( !sfxn->ready_for_nonideal_scoring() ) {
				TR.Warning << "options[ optimization::nonideal ] - scorefunction not set up for nonideal/Cartesian scoring (safe to ignore this warning if setting weights via RosettaScripts)" << std::endl;
			}
		}
	}
}

void
init_restore_score12prime( utility::options::OptionCollection & options ) {
	if ( options[ corrections::score::score12prime ] ) {
		revert_to_pre_talaris_2013_defaults( options );
		options[ score::weights ].value( "score12prime.wts" );
		options[ score::patch ].value( "" );
	}
}

void init_spades_score_function_correction( utility::options::OptionCollection & options )
{
	if ( options[ score::water_hybrid_sf ] ) {
		utility::vector1< std::string > params;
		if ( options[ basic::options::OptionKeys::chemical::set_atom_properties ].user() ) {
			TR.Warning << "Flag -score::water_hybrid_sf is set, but -set_atom_properties (or another flag such as -beta_jan25 that modifies this flag) is also set. Adding additional properties." << std::endl;
			params = options[ basic::options::OptionKeys::chemical::set_atom_properties ]();
			//for ( auto param : params ) {
			// TR << "Before spades adjustment: parameter " << param << std::endl;
			//}
		}
		params.push_back("fa_standard:OW:LK_DGFREE:0" );
		params.push_back("fa_standard:OW:LK_VOLUME:0" );
		options[ basic::options::OptionKeys::chemical::set_atom_properties ].value( params );
		//params = options[ basic::options::OptionKeys::chemical::set_atom_properties ]();
		//for ( auto param : params ) {
		// TR << "After spades adjustment: parameter " << param << std::endl;
		//}
	}
}


void
init_gen_potential_settings( utility::options::OptionCollection & options, bool setweight ) {
	if ( ! options[ corrections::gen_potential]() ) return;

	if ( ! options[ basic::options::OptionKeys::score::fa_Hatr ].user() ) {
		options[ basic::options::OptionKeys::score::fa_Hatr ].value(true);
	}

	// load proper ATS
	if ( ! options[ basic::options::OptionKeys::corrections::chemical::alternate_fullatom_ats ].user() ) {
		options[ basic::options::OptionKeys::corrections::chemical::alternate_fullatom_ats ].value( "fa_standard_genpot" );
	}

	// load new ring conf libs
	if ( ! options[ basic::options::OptionKeys::rings::ring_conformer_dbpath ].user() ) {
		options[ basic::options::OptionKeys::rings::ring_conformer_dbpath ].value( "chemical/ring_conformer_sets/fd_alternate" );
	}

	// gen_torsion behavior
	if ( options[ basic::options::OptionKeys::score::genbonded_score_hybrid ]() &&
			!(options[ basic::options::OptionKeys::score::count_pair_hybrid ].user() ||
			options[ basic::options::OptionKeys::score::count_pair_full ].user() ) ) {
		options[ basic::options::OptionKeys::score::count_pair_hybrid ].value( true );
		TR.Warning << "Flag -score::genbonded_score_hybrid set but count_pair behavior not defined; auto set to -score:count_pair_hybrid" << std::endl;
	}
	if ( options[ basic::options::OptionKeys::score::genbonded_score_full ]() &&
			!(options[ basic::options::OptionKeys::score::count_pair_hybrid ].user() ||
			options[ basic::options::OptionKeys::score::count_pair_full ].user() ) ) {
		options[ basic::options::OptionKeys::score::count_pair_full ].value( true );
		TR.Warning << "Flag -score::genbonded_score_full set but count_pair behavior not defined; auto set to -score:count_pair_full" << std::endl;
	}

	// modified countpair logic for ligands
	if ( ! options[ basic::options::OptionKeys::score::count_pair_hybrid ].user() &&
			! options[ basic::options::OptionKeys::score::count_pair_full ].user() ) {
		options[ basic::options::OptionKeys::score::count_pair_hybrid ].value( true );
	}

	// allow cyclic connection to be considered by cart_bonded
	if ( ! options[ basic::options::OptionKeys::score::cart_bonded_skip_cutpoints ].user() ) {
		options[ basic::options::OptionKeys::score::cart_bonded_skip_cutpoints ].value( false );
	}

	// hpark: recover original behavior as of May 23 2020
	/*
	// hbond params
	if ( !options[ basic::options::OptionKeys::corrections::score::lj_hbond_hdis ].user() ) {
	options[ basic::options::OptionKeys::corrections::score::lj_hbond_hdis ].value(2.3);
	}
	if ( !options[ basic::options::OptionKeys::corrections::score::lj_hbond_OH_donor_dis ].user() ) {
	options[ basic::options::OptionKeys::corrections::score::lj_hbond_OH_donor_dis ].value(3.3);
	}
	*/

	// consistency check
	{
		bool const bond_hybrid( options[ basic::options::OptionKeys::score::genbonded_score_hybrid ]() );
		bool const bond_full( options[ basic::options::OptionKeys::score::genbonded_score_full ]() );
		bool const cp_hybrid( options[ basic::options::OptionKeys::score::count_pair_hybrid ]() );
		bool const cp_full( options[ basic::options::OptionKeys::score::count_pair_full ]() );

		// assert not defined both together
		if ( cp_full && cp_hybrid ) {
			utility_exit_with_message("-count_pair_full/_hybrid cannot be used together!");
		}
		if ( bond_full && bond_hybrid ) {
			utility_exit_with_message("-genbonded_full/_hybrid cannot be used together!");
		}

		// assert not used as a mixture
		if ( (bond_full && cp_hybrid ) || (bond_hybrid && cp_full ) ) {
			utility_exit_with_message("[-count_pair_full/genbonded_pair_full] & [-count_pair_hybrid/genbonded_pair_hybrid] should be paired!");
		}

		// bond_full mode should be used -only if- cp_full is on (not the other way around; cp can work separately)
		if ( (bond_full && !cp_full) ) {
			utility_exit_with_message("-genbonded_score_full & -count_pair_full should be used together!");
		}
		// bond_hybrid mode should be used -only if- cp_hybrid is on (not the other way around; cp can work separately)
		if ( (bond_hybrid && !cp_hybrid) ) {
			utility_exit_with_message("-genbonded_score_hybrid & -count_pair_hybrid should be used together!");
		}
	}

	if ( setweight && !options[ score::weights ].user() ) {
		if ( options[corrections::beta_cart]() || options[optimization::nonideal]() ) {
			options[ score::weights ].value( "beta_genpot_cart.wts" );
		} else {
			options[ score::weights ].value( "beta_genpot.wts" );
		}
	} else {
		TR.Warning << "Flag -gen_potential is set but -weights are also specified.  Not changing input weights file!" << std::endl;
	}
}


void
init_dna_correction( utility::options::OptionCollection & options ) {

	utility::vector1< std::string > dna_residue_types;
	dna_residue_types.push_back("ADE");
	dna_residue_types.push_back("CYT");
	dna_residue_types.push_back("GUA");
	dna_residue_types.push_back("THY");

	if ( options[ corrections::newdna ] ) {

		{ // cloning atomtypes
			utility::vector1< std::string > mods;

			/// reassign the N9/N1 atoms from Ntrp atomtype to Npro atomtype
			// mods.push_back( "fa_standard:OOC:OPdd" ); // phosphate oxygens
			// mods.push_back( "fa_standard:OH:OEdd" ); // phosphate ester oxygens, was ONH2, doesnt make much difference
			// mods.push_back( "fa_standard:OH:ORdd" ); // ribose oxygen
			mods.push_back( "fa_standard:Phos:Pdna" ); // Phos in phosphate backbone
			mods.push_back( "fa_standard:OCbb:OBdd" ); // carbonyl oxygens in the bases

			if ( options[ basic::options::OptionKeys::chemical::clone_atom_types ].user() ) {
				/// add these at the end so they take precedence
				utility::vector1< std::string > const user_mods
					( options[ basic::options::OptionKeys::chemical::clone_atom_types ]() );
				for ( Size i=1; i<= user_mods.size(); ++i ) {
					if ( std::find( mods.begin(), mods.end(), user_mods[i] ) == mods.end() ) {
						mods.push_back( user_mods[i] );
					}
				}
			}
			options[ basic::options::OptionKeys::chemical::clone_atom_types ].value( mods );

		}

		{ // reassigning atomtype properties
			utility::vector1< std::string > mods;
			mods.push_back( "fa_standard:Pdna:LK_DGFREE:-4.1" ); // from way back, Jim H., based on S
			mods.push_back( "fa_standard:Pdna:LK_VOLUME:14.7" ); // -- ditto --

			if ( options[ basic::options::OptionKeys::chemical::set_atom_properties ].user() ) {
				/// add these at the end so they take precedence
				utility::vector1< std::string > const user_mods
					( options[ basic::options::OptionKeys::chemical::set_atom_properties ]() );
				for ( Size i=1; i<= user_mods.size(); ++i ) {
					if ( std::find( mods.begin(), mods.end(), user_mods[i] ) == mods.end() ) {
						mods.push_back( user_mods[i] );
					}
				}
			}
			options[ basic::options::OptionKeys::chemical::set_atom_properties ].value( mods );
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

			if ( options[ basic::options::OptionKeys::chemical::reassign_atom_types ].user() ) {
				/// add these at the end so they take precedence
				utility::vector1< std::string > const user_mods
					( options[ basic::options::OptionKeys::chemical::reassign_atom_types ]() );
				for ( Size i=1; i<= user_mods.size(); ++i ) {
					if ( std::find( mods.begin(), mods.end(), user_mods[i] ) == mods.end() ) {
						mods.push_back( user_mods[i] );
					}
				}
			}
			options[ basic::options::OptionKeys::chemical::reassign_atom_types ].value( mods );
		}

		{ // reassigning icoor values

			utility::vector1< std::string > mods;

			for ( Size ii=1; ii<= dna_residue_types.size(); ++ii ) {
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


			if ( options[ basic::options::OptionKeys::chemical::reassign_icoor ].user() ) {
				/// add these at the end so they take precedence
				utility::vector1< std::string > const user_mods
					( options[ basic::options::OptionKeys::chemical::reassign_icoor ]() );
				for ( Size i=1; i<= user_mods.size(); ++i ) {
					if ( std::find( mods.begin(), mods.end(), user_mods[i] ) == mods.end() ) {
						mods.push_back( user_mods[i] );
					}
				}
			}
			options[ basic::options::OptionKeys::chemical::reassign_icoor ].value( mods );
		}


	}


}

} // namespace
} // namespace
