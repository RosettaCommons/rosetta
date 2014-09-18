// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/init/score_function_corrections.cc
/// @brief  Initialize Score Funciton Corrections
/// @author Matthew O'Meara
///

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

void
init_crystal_refinement_correction() {
	//fpd  crystal-refinement specific changes
	if( option[cryst::crystal_refine]) {
		// use -correct icoor
		if ( ! option[ corrections::chemical::icoor_05_2009 ].user() ) {
			option[corrections::chemical::icoor_05_2009].value( true );
		}

		// use -correct rama fixes
		if ( ! option[ corrections::score::rama_not_squared ].user() ) {
			option[corrections::score::rama_not_squared].value( true );
		}
		if ( ! option[ corrections::score::rama_map ].user() ) {
			option[ corrections::score::rama_map ].value("scoring/score_functions/rama/Rama09_noEH_kernel25_it08.dat");
		}

		// use dun10
		if ( ! option[ corrections::score::dun10 ].user() ) {
			option[corrections::score::dun10].value( true );
		}

		// use bicubic interpolation
		if ( ! option[ corrections::score::use_bicubic_interpolation ].user() ) {
			option[corrections::score::use_bicubic_interpolation].value( true );
		}

		// read pdbs properly
		if ( ! option[ in::missing_density_to_jump ].user() ) {
			option[in::missing_density_to_jump].value( true );
		}
		if ( ! option[ in::preserve_crystinfo ].user() ) {
			option[in::preserve_crystinfo].value( true );
		}

		// special case for centroid scoring
		if ( ! option[ out::file::no_output_cen ].user() ) {
			option[out::file::no_output_cen].value( true );
		}

		// internal coordinate minimization: behave properly
		if ( ! option[ optimization::scale_rb ].user() ) {
			option[optimization::scale_rb].value( 1.0 );
		}

	}
}


void
init_beta_correction() {
	// beta score function - experimental energy function changes are accumulated here until made default
	// please share issues/feedback with Patrick Conway - ptconway@gmail.com
  if( option[corrections::beta]) {
		//nonideal
		if ( ! option[ optimization::nonideal ].user() ) {
      option[ optimization::nonideal ].value( true );
    }
		//optimized cart_bonded parameters - currently same as default
		if ( option[ optimization::nonideal ] ) {
			option[ score::bonded_params_dir ].value( "scoring/score_functions/bondlength_bondangle_beta" );
		}
		/*	--coming soon--
		//length dependent hbond_sr_bb to model helix cooperativity
		if ( ! option[ corrections::nonideal ].user() ) {
      option[ corrections::nonideal ].value( true );
    }
		//arovirt - aromatic cation-pi interaction modeling via virtual atoms
		if ( ! option[ corrections::nonideal ].user() ) {
      option[ corrections::nonideal ].value( true );
    }
		//dunbrack hbond exclusion - downweight rotamers that hbond with self/neighbor (don't double count hbond terms)
		if ( ! option[ corrections::nonideal ].user() ) {
      option[ corrections::nonideal ].value( true );
    }
		//upweight serine/threonine H-HG/H-HH intra_rep so they don't clash
		if ( ! option[ corrections::nonideal ].user() ) {
      option[ corrections::nonideal ].value( true );
    }
		*/

		// talaris2013_beta.wts - currently same as talaris2013_cart.wts
		// overrides weight files set in the other options - make sure all of these previous options are happy!
    if ( ! option[ score::weights ].user() ) {
      option[ score::weights ].value( "talaris2013_beta.wts" );
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

void
init_score_function_corrections(){
	init_revert_to_pre_talaris_2013_mistake();
	init_hbond_sp2_correction();
	init_facts_correction();
	init_correct_correction();
	init_crystal_refinement_correction();
	init_beta_correction();
	init_nonideal_correction();
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

} // namespace
} // namespace
