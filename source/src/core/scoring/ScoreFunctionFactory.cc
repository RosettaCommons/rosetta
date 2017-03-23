// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/conformation/ScoreFunctionFactory.cc
/// @brief Manages generation of ScoreFunction objects from weights and patch files, RosettaScripts tags, etc
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/DockingScoreFunction.hh>
#include <core/scoring/MinScoreScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/chemical/AA.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/keys/OptionKeyList.hh>

// option key includes
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/mistakes.OptionKeys.gen.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <sstream>

using basic::T;

static THREAD_LOCAL basic::Tracer TR( "core.scoring.ScoreFunctionFactory" );

namespace core {
namespace scoring {

ScoreFunctionOP
ScoreFunctionFactory::create_score_function( std::string const & weights_tag )
{
	return create_score_function( basic::options::option, weights_tag );
}

ScoreFunctionOP
ScoreFunctionFactory::create_score_function( utility::options::OptionCollection const & options, std::string const & weights_tag )
{
	utility::vector1< std::string > patch_tags;
	return create_score_function( options, weights_tag, patch_tags );
}


ScoreFunctionOP
ScoreFunctionFactory::create_score_function( std::string const & weights_tag, utility::vector1< std::string > const & patch_tags ) {
	return create_score_function( basic::options::option, weights_tag, patch_tags );
}

ScoreFunctionOP
ScoreFunctionFactory::create_score_function(
	utility::options::OptionCollection const & options,
	std::string const & weights_tag_in,
	utility::vector1< std::string > const & patch_tags_in
) {
	using namespace basic::options::OptionKeys;
	std::string weights_tag( weights_tag_in );
	utility::vector1< std::string > patch_tags( patch_tags_in ); // copy the input patches; we're going to modify them

	// create a new scorefunction
	ScoreFunctionOP scorefxn;
	if ( options[ score::min_score_score ].user() ) {
		scorefxn = ScoreFunctionOP( new MinScoreScoreFunction( options[ score::min_score_score ]() ) );
	} else if ( options[ score::docking_interface_score ]() ) {
		scorefxn = ScoreFunctionOP( new DockingScoreFunction );
	} else if ( options[ basic::options::OptionKeys::symmetry::symmetry_definition ].user() ) {
		scorefxn = ScoreFunctionOP( new SymmetricScoreFunction( options ) );
	} else {
		scorefxn = ScoreFunctionOP( new ScoreFunction( options ) );
	}

	// Avoid loading the score12 patch if we're using
	// 1) standard weights,
	// 2) the score12 patch, and
	// 3) the flag "score12prime"
	if ( weights_tag == PRE_TALARIS_2013_STANDARD_WTS &&
			options[ corrections::score::score12prime ] ) {
		bool sc12patch = false;
		for ( Size ii = 1; ii <= patch_tags.size(); ++ii ) {
			if ( patch_tags[ ii ] == SCORE12_PATCH ) {
				patch_tags[ ii ] = "NOPATCH";
				sc12patch = true;
			}
		}
		if ( sc12patch ) {
			weights_tag = "score12prime";
		}
	}

	runtime_assert(validate_beta(weights_tag, options));

	load_weights_file( weights_tag, scorefxn );

	for ( utility::vector1< std::string >::const_iterator it = patch_tags.begin(); it != patch_tags.end(); ++it ) {
		std::string const& patch_tag( *it );
		if ( patch_tag.size() && patch_tag != "NOPATCH" ) {
			//   TR.Debug << "SCOREFUNCTION: apply patch "  << patch_tag << std::endl;
			scorefxn->apply_patch_from_file( patch_tag );
		}
	}

	// allow user to change weights via options system
	apply_user_defined_reweighting_( options, scorefxn );

	scorefxn->name( weights_tag );
	return scorefxn;
}


/// @details If requested tag is beta_nov15 or beta_july15, but the user did not
/// pass the relevant options-system option, ERROR!  This is because this
/// scorefunction family has overrides to parameters (LK solvation params, etc).
/// Those are loaded from the command-line flag, not the weights file. Using
/// only the weights file will give you mismatched weights/params and much
/// sadness.
bool
ScoreFunctionFactory::validate_beta(
	std::string const & weights_tag,
	utility::options::OptionCollection const & options
)
{

	bool const sf_maybe_beta(weights_tag.find("beta") != std::string::npos);
	core::Size const weights_length(weights_tag.length());

	//if the scorefunction does't have beta in it, we're OK.
	//If the length is less than 4, also abort, because later substr() operations
	//will go out of range.
	//(As the code stands this latter check is irrelevant, but I'm leaving it in
	//in case we check for a string that isn't literally "beta", which happens to
	//be the same char length as ".wts"
	if ( !sf_maybe_beta || weights_length < 4 ) return true;

	//determine if weights_tag has .wts or not
	//This Size is protected from overflow by the above if < 4
	core::Size const fname_length(weights_tag.length() - 4); //4 represents ".wts"
	std::string const weights_tag_extension(weights_tag.substr(fname_length)); // might or might not actually be an extension
	bool const weights_tag_has_extension(weights_tag_extension == ".wts");
	//this ternary creates the de-extended weights tag if it was extended
	std::string const weights_tag_no_extension(weights_tag_has_extension ? weights_tag.substr(0, fname_length): weights_tag);

	//Determine which user options are active and concerning
	//Note this is checking a function argument options, not the global options
	//not sure why create_score_function is like that.
	using namespace basic::options::OptionKeys;
	bool const betanov16_active(options[corrections::beta_nov16].value()
		|| options[corrections::beta_nov16_cart].value() );
	bool const betanov15_active(options[corrections::beta_nov15].value()
		|| options[corrections::beta_nov15_cart].value() );
	bool const betajuly15_active(options[corrections::beta_july15].value()
		|| options[corrections::beta_july15_cart].value() );

	if ( (weights_tag_no_extension == (BETA_NOV16)) && !betanov16_active ) {
		utility_exit_with_message(BETA_NOV16 + "(.wts) requested, but -corrections::beta_nov16 not set to true. This leads to a garbage scorefunction.  Exiting.");
		return false; //can't get here
	} else if ( (weights_tag_no_extension == (BETA_NOV15)) && !betanov15_active ) {
		utility_exit_with_message(BETA_NOV15 + "(.wts) requested, but -corrections::beta_nov15 not set to true. This leads to a garbage scorefunction.  Exiting.");
		return false; //can't get here
	} else if ( (weights_tag_no_extension == (BETA_JULY15)) && !betajuly15_active ) {
		utility_exit_with_message(BETA_JULY15 + "(.wts) requested, but -corrections::beta_july15 not set to true. This leads to a garbage scorefunction.  Exiting.");
		return false; //can't get here
	} else if ( sf_maybe_beta && !betanov16_active && !betanov15_active && !betajuly15_active ) {
		TR.Warning << "**************************************************************************\n"
			<< "*****************************************************\n"
			<< "****************************************************\n"
			<< weights_tag << " may be a 'beta' scorefunction, but ScoreFunctionFactory thinks the beta flags weren't set.  "
			<< "Your scorefunction may be garbage!\n"
			<< "**************************************************************************\n"
			<< "*****************************************************\n"
			<< "****************************************************" << std::endl;
		//return something between true and false, if there was such a thing
	}

	return true;
} //ScoreFunctionFactory::validate_beta

ScoreFunctionOP
ScoreFunctionFactory::create_score_function( std::string const & weights_tag, std::string const & patch_tag )
{
	return create_score_function( basic::options::option, weights_tag, patch_tag );
}

ScoreFunctionOP
ScoreFunctionFactory::create_score_function(
	utility::options::OptionCollection const & options,
	std::string const & weights_tag,
	std::string const & patch_tag )
{
	utility::vector1< std::string > patch_tags;
	patch_tags.push_back( patch_tag );
	return create_score_function( options, weights_tag, patch_tags );
}

/// @brief A documentation function which reports the set of options read by the create_score_function variants
void
ScoreFunctionFactory::list_read_options( utility::options::OptionKeyList & opts )
{
	using namespace basic::options::OptionKeys;

	MinScoreScoreFunction::list_options_read( opts );
	DockingScoreFunction::list_options_read( opts );
	ScoreFunction::list_options_read( opts );
	symmetry::SymmetricScoreFunction::list_options_read( opts );

	opts
		+ score::min_score_score
		+ score::docking_interface_score
		+ basic::options::OptionKeys::symmetry::symmetry_definition
		+ corrections::score::score12prime
		+ score::set_weights
		+ abinitio::rg_reweight
		+ score::ref_offset
		+ score::ref_offsets
		+ basic::options::OptionKeys::corrections::beta_nov15
		+ basic::options::OptionKeys::corrections::beta_nov15_cart
		+ basic::options::OptionKeys::corrections::beta_july15
		+ basic::options::OptionKeys::corrections::beta_july15_cart;
}

void ScoreFunctionFactory::apply_user_defined_reweighting_(
	utility::options::OptionCollection const & options,
	core::scoring::ScoreFunctionOP scorefxn
) {
	// do some reweighting here. This code could be much more simple if the options system could
	// produce a std::pair< std::string, core::Real >. For now these are separate options for
	// the different reweights.
	using namespace basic::options::OptionKeys;

	/// new mechanism: set multiple weights using string vector option
	if ( options[ score::set_weights ].user() ) apply_set_weights( scorefxn, options[ score::set_weights ]() );

	if ( options[ abinitio::rg_reweight ].user() ) {
		scorefxn->set_weight( rg, scorefxn->get_weight( rg ) * options[ abinitio::rg_reweight ]() );
	}
	// offset reference energies using user options, for example: -score:ref_offsets TRP 0.9 HIS 0.3
	if ( !options[ score::ref_offsets ].user() && !options[ score::ref_offset ].user() ) return;

	// get the ref weights from the EnergyMethodOptions object
	methods::EnergyMethodOptions energy_method_options(scorefxn->energy_method_options());
	if ( !energy_method_options.has_method_weights(ref) ) {
		// utility_exit_with_message("option -score:ref_offsets requires preexisting reference energies");
	} else {
		utility::vector1<core::Real> ref_weights(energy_method_options.method_weights(ref));

		if ( options[ score::ref_offset ].user() ) {
			Real const offset = options[ score::ref_offset ]();
			for ( Size n = 1; n <= ref_weights.size(); n++ ) {
				ref_weights[ n ] += offset;
			}
		} else {
			runtime_assert(  options[ score::ref_offsets ].user() );
			// get the offsets vector and make sure it contains pairs
			utility::vector1<std::string> const & ref_offsets( options[ score::ref_offsets ]() );
			if ( ref_offsets.size() % 2 != 0 ) {
				utility_exit_with_message("option -score:ref_offsets requires pairs of 3 character residue types and offsets");
			}

			// iterate over all pairs
			for ( auto iter(ref_offsets.begin()), iter_end(ref_offsets.end());
					iter != iter_end; ++iter ) {
				// get the aa type from the pair
				std::istringstream aa_iss(*iter);
				core::chemical::AA aa;
				if ( !(aa_iss >> aa) ) {
					utility_exit_with_message(aa_iss.str()+" is not a valid 3 character residue type for -score:ref_offsets");
				}
				// get the offset from the pair
				std::istringstream offset_iss(*(++iter));
				core::Real offset;
				if ( !(offset_iss >> offset) ) {
					utility_exit_with_message(offset_iss.str()+" is not a valid offset for -score:ref_offsets");
				}
				// offset the weight
				ref_weights[aa] += offset;
			}
		}
		// load the ref weights back into the EnergyMethodOptions object
		energy_method_options.set_method_weights(ref, ref_weights);
		scorefxn->set_energy_method_options(energy_method_options);
	}
}


void ScoreFunctionFactory::load_weights_file( std::string weights_tag, ScoreFunctionOP scorefxn )
{
	scorefxn->initialize_from_file(weights_tag);
}

std::string const TALARIS_2014( "talaris2014" );
std::string const TALARIS_2013( "talaris2013" );
std::string const TALARIS_2013_CART( "talaris2013_cart" );
std::string const PRE_TALARIS_2013_STANDARD_WTS( "pre_talaris_2013_standard" );
std::string const CENTROID_WTS( "cen_std" );
std::string const SOFT_REP_WTS( "soft_rep" );
std::string const SOFT_REP_DESIGN_WTS( "soft_rep_design" );
std::string const DNA_INT_WTS( "dna_no_gb" );
std::string const DNA_INT_WTS_GB( "dna" );
std::string const MM_STD_WTS( "mm_std" );
std::string const RNA_LORES_WTS( "rna/denovo/rna_lores" );
std::string const RNA_HIRES_WTS( "rna/denovo/rna_hires" );
std::string const RNA_LORES_PLUS_HIRES_WTS( "rna/denovo/rna_lores_plus_hires" );
std::string const MEMB_HIGHRES_WTS( "membrane_highres" ); //pba

std::string const SCORE12_PATCH( "score12" );
std::string const SCORE13( "score13" );
std::string const DOCK_PATCH( "docking" );
std::string const DOCK_LOW_PATCH( "docking_cen" );

std::string const SCORE4_SMOOTH_CART( "score4_smooth_cart" );

std::string const BETA_NOV16( "beta_nov16" );
std::string const BETA_NOV15( "beta_nov15" );
std::string const BETA_JULY15( "beta_july15" );

core::scoring::ScoreFunctionOP
get_score_function( bool const is_fullatom /* default true */ )
{
	return get_score_function( basic::options::option, is_fullatom );
}


core::scoring::ScoreFunctionOP
get_score_function(
	utility::options::OptionCollection const & options,
	bool const is_fullatom
)
{
	using namespace basic::options::OptionKeys;

	if ( options[ score::empty ]() ) return core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction( options ) );

	std::string weight_set = options[ score::weights ];
	utility::vector1< std::string > patch_tags = options[ score::patch ]();

	if ( !options[ score::weights ].user() && !is_fullatom ) {
		// Defalt score of centroid is cen_wts when is_fullatom is false and user has not specified a score weights
		weight_set = CENTROID_WTS;
	} else {

		// Default score is talaris2014 if the user has not specified a score weights file or a patch file
		// on the command line.

		//TR << "get_score_function1: weight set " << weight_set << " same ? " << ( weight_set == "pre_talaris_2013_standard.wts") << std::endl;
		//TR << "options[ score::weights ].user() ? " << options[ score::weights ].user() << std::endl;
		//TR << "options[ score::patch ].user() ? " << options[ score::patch ].user() << std::endl;

		if ( ( weight_set == "pre_talaris_2013_standard.wts" && !options[ score::weights ].user() ) &&
				( !options[ score::patch ].user() ) ) {
			patch_tags.push_back( "score12" );
			//TR << "pushing back score12 patch" << std::endl;
			if ( options[corrections::correct] ) {
				//TR << "setting weight set to score12_w_corrections" << std::endl;
				weight_set = "score12_w_corrections";
				patch_tags.clear();
			} else if ( options[ corrections::hbond_sp2_correction ] ) {
				patch_tags.clear();
				if ( options[ corrections::score::dun10 ] ) {
					weight_set = "sp2_correction";
				} else {
					weight_set = "sp2_correction_dun02";
				}
			} else if ( options[ corrections::score::score12prime ] ) {
				weight_set = "score12prime";
				patch_tags.clear();
			}
		}
		//TR << "get_score_function2: weight set " << weight_set << std::endl;
	}

	T("core.scoring.ScoreFunctionFactory") << "SCOREFUNCTION: " << utility::CSI_Green() << weight_set << utility::CSI_Reset() << std::endl;

	core::scoring::ScoreFunctionOP scorefxn;

	//create_score_function() handles symmetry, etc.
	if ( patch_tags.size() == 0 ) {
		scorefxn = scoring::ScoreFunctionFactory::create_score_function( weight_set );
	} else {
		for ( Size ii = 1; ii <= patch_tags.size(); ++ii ) {
			if ( patch_tags[ii]!="" ) {
				T("core.scoring.ScoreFunctionFactory") << "SCOREFUNCTION PATCH: " << patch_tags[ii] << std::endl;
			}
		}
		scorefxn = scoring::ScoreFunctionFactory::create_score_function( weight_set, patch_tags );
	}

	// add in constraint weights if specified by the user.
	// mtyka: No No No we dont want this here. Add constraints and weights ouside this function.
	//  just setting the weights alone isnt gonna get oyu far anyway.
	// VKM: There are special cases in which constraint weights need to be turned on, I think (e.g.
	// for metalloproteins), and this is the best place to ensure that it happens consistently.

	//Turn on constraints if the user has used the auto_setup_metals flag.  Constraints are added automatically on PDB import.
	if ( options[in::auto_setup_metals].user() ) {
		if ( scorefxn->get_weight(metalbinding_constraint) < 1.0e-10 ) {
			T("core.scoring.ScoreFunctionFactory") << "The -auto_setup_metals flag was used with no metalbinding_constraint weight set in the weights file.  Setting to 1.0." << std::endl ;
			scorefxn->set_weight(metalbinding_constraint, 1.0); // Turn on the atom_pair_constraint weight if and only if it isn't already turned on.
			// If it is already turned on, then the automatic constraint adder will adjust constraint strengths appropriately, which means that we
			// don't need to set this to 1.0 -- any nonzero value is fine.
		}
	}

	// Turn on carbohydrate energy method weights if the user has supplied the -include_sugars flag.
	if ( options[ in::include_sugars ].user() ) {
		if ( TR.Info.visible() ) {
			TR.Info << "The -include_sugars flag was used with no sugar_bb weight set in the weights file.  " <<
				"Setting sugar_bb weight to 1.0 by default." << std::endl;
		}
		scorefxn->set_weight( sugar_bb, 1.0);
	}

	return scorefxn;
}

void
list_read_options_in_get_score_function( utility::options::OptionKeyList & opts )
{
	using namespace basic::options::OptionKeys;

	// the call to create_score_function reads a significant number of options
	ScoreFunctionFactory::list_read_options( opts );

	opts
		+ score::empty
		+ score::weights
		+ score::patch
		+ corrections::correct
		+ corrections::hbond_sp2_correction
		+ corrections::score::dun10
		+ corrections::score::score12prime
		+ in::auto_setup_metals
		//+ in::metals_distance_constraint_multiplier
		//+ in::metals_angle_constraint_multiplier
		+ in::include_sugars;
}

core::scoring::ScoreFunctionOP get_score_function_legacy(
	std::string pre_talaris_2013_weight_set,
	std::string pre_talaris_2013_patch_file
)
{
	return get_score_function_legacy( basic::options::option, pre_talaris_2013_weight_set, pre_talaris_2013_patch_file );
}
core::scoring::ScoreFunctionOP get_score_function_legacy(
	utility::options::OptionCollection const & options,
	std::string pre_talaris_2013_weight_set,
	std::string pre_talaris_2013_patch_file
)
{
	if ( options[ basic::options::OptionKeys::mistakes::restore_pre_talaris_2013_behavior ] ) {
		return ScoreFunctionFactory::create_score_function( options, pre_talaris_2013_weight_set, pre_talaris_2013_patch_file );
	}
	return get_score_function( options );
}

/// @brief A documentation function which reports the set of options read by get_score_function_legacy.
void
list_read_options_in_get_score_function_legacy( utility::options::OptionKeyList & opts )
{
	ScoreFunctionFactory::list_read_options( opts );
	opts + basic::options::OptionKeys::mistakes::restore_pre_talaris_2013_behavior;
}

std::string
get_score_functionName(
	bool const is_fullatom /* default true */
) {
	using basic::options::option;
	using namespace basic::options::OptionKeys;

	if ( option[ score::empty ]() ) return "empty";

	std::string weight_set = option[ score::weights ];
	utility::vector1< std::string > patch_tags = option[ score::patch ]();


	if ( !option[ score::weights ].user() && !is_fullatom ) {
		// Default score of centroid is cen_wts when is_fullatom is false and user has not specified a score weights
		weight_set = CENTROID_WTS;
	} else {
		/// Default score is score12 if the user has not specified a score weights file or a patch file
		/// on the command line.  If the user has specified that they would like the standard weight set,
		/// and has not also asked for the score12 patch, then do not apply the score12 patch to it.
		if ( ( weight_set == "pre_talaris_2013_standard" && !option[ score::weights ].user() ) &&
				( !option[ score::patch ].user() ) ) {
			patch_tags.push_back( "score12" );
			if ( option[ corrections::correct ] ) {
				weight_set = "score12_w_corrections";
				patch_tags.clear();
			} else if ( option[ corrections::score::score12prime ] ) {
				weight_set = "score12prime";
				patch_tags.clear();
			}
		}
	}

	if ( patch_tags.size() != 0 &&
			weight_set == PRE_TALARIS_2013_STANDARD_WTS &&
			option[ corrections::score::score12prime ] ) {

		bool sc12patch = false;
		for ( Size ii = 1; ii <= patch_tags.size(); ++ii ) {
			if ( patch_tags[ ii ] == SCORE12_PATCH ) {
				patch_tags[ ii ] = "NOPATCH";
				sc12patch = true;
			}
		}
		if ( sc12patch ) {
			weight_set = "score12prime";
		}
	}

	std::stringstream patch_string;
	for ( Size ii=1; ii <= patch_tags.size(); ++ii ) {
		if ( patch_tags[ii] == "NOPATCH" ) continue;
		patch_string << "_" << patch_tags[ii];
	}

	return weight_set + patch_string.str();
}

void
apply_set_weights( ScoreFunctionOP scorefxn, utility::vector1< std::string > const & settings )
{
	std::string const errmsg("proper format for -set_weights is a list of paired strings, e.g: '-set_weights fa_atr 0.6 -fa_rep 0.55 -fa_sol 0.9' ");
	if ( settings.size()%2 != 0 ) utility_exit_with_message( errmsg );
	for ( Size i=0; i< settings.size()/2; ++i ) {
		if ( !ObjexxFCL::is_float( settings[ 2*i+2 ] ) ) utility_exit_with_message( errmsg );
		ScoreType const t( score_type_from_name( settings[ 2*i + 1] ) );
		Real const value( ObjexxFCL::float_of( settings[ 2*i + 2 ] ) );
		TR << "Setting/modifying scorefxn weight from command line: " << t << ' ' << value << std::endl;
		scorefxn->set_weight( t, value );
	}
}


} // namespace scoring
} // namespace core
