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
#include <core/scoring/DockingScoreFunction.hh>
#include <core/scoring/MinScoreScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/symmetry/util.hh>
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
#include <basic/citation_manager/CitationManager.hh>
#include <sstream>

#include <utility/thread/threadsafe_creation.hh>
#include <utility/string_util.hh>

static basic::Tracer TR( "core.scoring.ScoreFunctionFactory" );

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
	return ScoreFunctionFactory::get_instance()->create_score_function_nonstatic( options, weights_tag_in, patch_tags_in );
}

/// @details If requested tag is talaris2013/talaris2014, but the user did not
/// pass the relevant options-system option, ERROR! This is because this
/// scorefunction family has overrides to parameters (LK solvation params, etc).
/// Those are loaded from the command-line flag, not the weights file. Using
/// only the weights file will give you mismatched weights/params and much
/// sadness.
bool
ScoreFunctionFactory::validate_talaris(
	std::string const & weights_tag,
	utility::options::OptionCollection const & options
)
{
	bool sf_maybe_talaris(weights_tag.find("talaris") != std::string::npos);
	if ( weights_tag.find("pre_talaris") != std::string::npos ) {
		sf_maybe_talaris = false;
	}
	core::Size const weights_length(weights_tag.length());

	if ( !sf_maybe_talaris || weights_length < 11 ) return true;

	core::Size const fname_length(weights_tag.length() - 4); //4 represents ".wts"
	std::string const weights_tag_extension(weights_tag.substr(fname_length)); // might or might not actually be an extension
	bool const weights_tag_has_extension(weights_tag_extension == ".wts");
	//this ternary creates the de-extended weights tag if it was extended
	std::string const weights_tag_no_extension(weights_tag_has_extension ? weights_tag.substr(0, fname_length): weights_tag);

	using namespace basic::options::OptionKeys;
	bool const talaris_active( options[corrections::restore_talaris_behavior].value());

	if ( (weights_tag_no_extension == (TALARIS_2014)) && !talaris_active ) {
		utility_exit_with_message(TALARIS_2014 + "(.wts) requested, but -corrections::restore_talaris_behavior not set to true. This leads to a garbage scorefunction.  Exiting.");
		return false; //can't get here
	} else if ( (weights_tag_no_extension == (TALARIS_2013)) && !talaris_active ) {
		utility_exit_with_message(TALARIS_2013 + "(.wts) requested, but -corrections::restore_talaris_behavior not set to true. This leads to a garbage scorefunction.  Exiting.");
		return false; //can't get here
	} else if ( sf_maybe_talaris && !talaris_active ) {
		TR.Warning << "**************************************************************************\n"
			<< "*****************************************************\n"
			<< "****************************************************\n"
			<< weights_tag << " may be a 'talaris' scorefunction, but ScoreFunctionFactory thinks the -restore_talaris_behavior flags weren't set.  "
			<< "Your scorefunction may be garbage!\n"
			<< "**************************************************************************\n"
			<< "*****************************************************\n"
			<< "****************************************************" << std::endl;
	}

	return true;
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
	bool const genpot_active(options[corrections::gen_potential].value());
	bool const betanov16_active(options[corrections::beta_nov16].value()
		|| options[corrections::beta_nov16_cart].value() );
	bool const betajuly15_active(options[corrections::beta_july15].value()
		|| options[corrections::beta_july15_cart].value() );
	if ( (weights_tag_no_extension == (BETA_GENPOT)) && !genpot_active ) {
		utility_exit_with_message(BETA_GENPOT + "(.wts) requested, but -corrections::gen_potential not set to true. This leads to a garbage scorefunction.  Exiting.");
		return false; //can't get here
	} else if ( (weights_tag_no_extension == (BETA_NOV16)) && !betanov16_active ) {
		utility_exit_with_message(BETA_NOV16 + "(.wts) requested, but -corrections::beta_nov16 not set to true. This leads to a garbage scorefunction.  Exiting.");
		return false; //can't get here
	} else if ( (weights_tag_no_extension == (BETA_JULY15)) && !betajuly15_active ) {
		utility_exit_with_message(BETA_JULY15 + "(.wts) requested, but -corrections::beta_july15 not set to true. This leads to a garbage scorefunction.  Exiting.");
		return false; //can't get here
	} else if ( sf_maybe_beta && !betanov16_active && !genpot_active ) {
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

/********************** PRIVATE MEMBER FUNCTIONS *******************************/


/// @brief Nonstatic version for storing the loaded scorefunction in this object.
/// @details Loads scorefunction from disk if it has not yet been loaded (and caches it), or retrieves cached copy
/// if it has been loaded already.  Returns clone of cached copy.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
ScoreFunctionOP
ScoreFunctionFactory::create_score_function_nonstatic(
	utility::options::OptionCollection const & options,
	std::string const & weights_tag,
	utility::vector1< std::string > const & patch_tags
) {
	// The following is for threadsafe insertion of an entry into a map.  We create a std::function object, which allows
	// a function to be packaged with a set of parameters passed to it, for calling later.  We then pass it to a function
	// located in utility::thread which (a) locks the map mutex (to ensure that only one thread does the subsequent steps),
	// (b) checks whether the map already contains the key (e.g. if anothe thread has already inserted it into the map),
	// (c) calls the function if and only if the key is not in the map (so that only one thread ever loads the data from
	// disk, and this is not repeated by a second thread that might be waiting on the mutex, (d) adds the key-value pair
	// to the map, and (e) unlocks the mutex.  Any following threads waiting on the mutex will discover that the key is
	// already in the map, and will not try to load the data again or add it again:
	std::function< core::scoring::ScoreFunctionOP () > creator(
		std::bind( &ScoreFunctionFactory::load_score_function_from_disk, std::cref( options ), std::cref( weights_tag ), std::cref( patch_tags ) )
	);

	std::ostringstream concatenated_patch_tags;
	for ( core::Size i(1), imax(patch_tags.size()); i<=imax; ++i ) {
		concatenated_patch_tags << patch_tags[i];
		if ( i<imax ) {
			concatenated_patch_tags << ",";
		}
	}

	return ( utility::thread::safely_check_map_for_key_and_insert_if_absent(
		creator,
		SAFELY_PASS_MUTEX( loaded_scorefxns_mutex_ ),
		ScoreFunctionKey( weights_tag, concatenated_patch_tags.str(), &options ),
		loaded_scorefxns_
		)
		)->clone();
}

/// @brief Load a scorefunction from disk.
/// @details TRIGGERS READ FROM DISK!  This function is needed for threadsafe creation.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
/*static*/
core::scoring::ScoreFunctionOP
ScoreFunctionFactory::load_score_function_from_disk(
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
		scorefxn = utility::pointer::make_shared< MinScoreScoreFunction >( options[ score::min_score_score ]() );
	} else if ( options[ score::docking_interface_score ]() ) {
		scorefxn = utility::pointer::make_shared< DockingScoreFunction >();
	} else {
		scorefxn = utility::pointer::make_shared< ScoreFunction >( options );
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

	runtime_assert(validate_talaris(weights_tag, options) );
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

	// Register this scorefunction with the citation manager:
	basic::citation_manager::CitationCollectionList citations;
	scorefxn->provide_citation_info( citations );
	basic::citation_manager::CitationManager::get_instance()->add_citations( citations );

	return scorefxn;
}

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
	utility::vector1< std::string > const patch_tags {patch_tag};
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

	opts
		+ score::min_score_score
		+ score::docking_interface_score
		+ basic::options::OptionKeys::symmetry::symmetry_definition
		+ corrections::score::score12prime
		+ score::set_weights
		+ abinitio::rg_reweight
		+ score::ref_offset
		+ score::ref_offsets
		//+ basic::options::OptionKeys::corrections::beta_nov15
		//+ basic::options::OptionKeys::corrections::beta_nov15_cart
		+ basic::options::OptionKeys::corrections::beta_july15
		+ basic::options::OptionKeys::corrections::beta_july15_cart
		+ basic::options::OptionKeys::corrections::gen_potential
		+ basic::options::OptionKeys::corrections::beta_nov16
		+ basic::options::OptionKeys::corrections::beta_nov16_cart;

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

std::string const REF_2015( "ref2015" );
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

std::string const BETA_GENPOT( "beta_genpot" );
std::string const BETA_NOV16( "beta_nov16" );
std::string const BETA_NOV15( "beta_nov15" );
std::string const BETA_JULY15( "beta_july15" );

core::scoring::ScoreFunctionOP
get_score_function( bool const is_fullatom /* default true */ )
{
	return get_score_function( basic::options::option, is_fullatom );
}

core::scoring::ScoreFunctionOP
get_score_function( pose::Pose const & pose, bool const is_fullatom )
{
	return get_score_function(pose, basic::options::option, is_fullatom );

}

core::scoring::ScoreFunctionOP
get_score_function( pose::Pose const & pose, utility::options::OptionCollection const & options, bool const is_fullatom )
// Looks like this only exists to check if pose is symmetric.  Since scorefxn is symmetry-agnostic now, this could be removed
// Would require changing all instances of get_score_function to the get_score_function( options, is_fullatom ) signature.
{
	core::scoring::ScoreFunctionOP scorefxn =  get_score_function( options, is_fullatom );
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		return scorefxn;
	}
	return scorefxn;
}


core::scoring::ScoreFunctionOP
get_score_function(
	utility::options::OptionCollection const & options,
	bool const is_fullatom
)
{
	using namespace basic::options::OptionKeys;

	if ( options[ score::empty ]() ) return utility::pointer::make_shared< core::scoring::ScoreFunction >( options );

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

	TR << "SCOREFUNCTION: " << utility::CSI_Green() << weight_set << utility::CSI_Reset() << std::endl;

	core::scoring::ScoreFunctionOP scorefxn;

	//create_score_function() handles symmetry, etc.
	if ( patch_tags.size() == 0 ) {
		scorefxn = scoring::ScoreFunctionFactory::create_score_function( weight_set );
	} else {
		for ( Size ii = 1; ii <= patch_tags.size(); ++ii ) {
			if ( patch_tags[ii]!="" ) {
				TR << "SCOREFUNCTION PATCH: " << patch_tags[ii] << std::endl;
			}
		}
		scorefxn = scoring::ScoreFunctionFactory::create_score_function(options, weight_set, patch_tags );
	}

	// add in constraint weights if specified by the user.
	// mtyka: No No No we dont want this here. Add constraints and weights ouside this function.
	//  just setting the weights alone isnt gonna get oyu far anyway.
	// VKM: There are special cases in which constraint weights need to be turned on, I think (e.g.
	// for metalloproteins), and this is the best place to ensure that it happens consistently.

	//Turn on constraints if the user has used the auto_setup_metals flag.  Constraints are added automatically on PDB import.
	if ( options[in::auto_setup_metals].value() ) {
		if ( scorefxn->get_weight(metalbinding_constraint) < 1.0e-10 ) {
			TR << "The -auto_setup_metals flag was used with no metalbinding_constraint weight set in the weights file.  Setting to 1.0." << std::endl ;
			scorefxn->set_weight(metalbinding_constraint, 1.0); // Turn on the atom_pair_constraint weight if and only if it isn't already turned on.
			// If it is already turned on, then the automatic constraint adder will adjust constraint strengths appropriately, which means that we
			// don't need to set this to 1.0 -- any nonzero value is fine.
		}
	}

	// Turn on carbohydrate energy method weights if the user has supplied the -include_sugars flag.
	if ( options[ in::include_sugars ].value() && ! options[ score::force_sugar_bb_zero].value() ) {
		if ( TR.Info.visible() && scorefxn->get_weight( sugar_bb ) == 0 ) {

			TR.Info << "The -include_sugars flag was used with no sugar_bb weight set in the weights file.  " <<
				"Setting sugar_bb weight to 0.5 by default." << std::endl;
		}
		scorefxn->set_weight_if_zero( sugar_bb, 0.5); //JAB - only set the weight if we have not added or changed it already.
	}
	// JAB - turn on intra-rep to get less bad structures and energies.
	//  Also, recommend beta - as the LKBridge term helps sugars significantly.
	if ( options[ in::include_sugars].value() && scorefxn->get_weight( fa_intra_rep_xover4 ) == 0 ) {
		TR.Info << " The -include_sugars flag was used without fa_intra_rep_xover4 term in the scorefunction."<<
			" Setting this term's weight to 0.55. It is generally recommended to use the -beta scorefunction (Rosetta-ICO) with sugars,"<<
			" which includes this and other desired terms such as those bridging waters" << std::endl;
		scorefxn->set_weight_if_zero(fa_intra_rep_xover4, 0.55);
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
		+ in::include_sugars
		+ score::force_sugar_bb_zero;
}

core::scoring::ScoreFunctionOP get_score_function_legacy(
	std::string const & pre_talaris_2013_weight_set,
	std::string const & pre_talaris_2013_patch_file
)
{
	return get_score_function_legacy( basic::options::option, pre_talaris_2013_weight_set, pre_talaris_2013_patch_file );
}
core::scoring::ScoreFunctionOP get_score_function_legacy(
	utility::options::OptionCollection const & options,
	std::string const & pre_talaris_2013_weight_set,
	std::string const &  pre_talaris_2013_patch_file
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
	std::string result = weight_set + patch_string.str();

	//result sometimes includes extensions. For consistency, let's remove them
	auto index_of_first_period = result.find_first_of(".");
	if ( index_of_first_period != std::string::npos ) {
		return result.substr( 0, index_of_first_period );
	}

	//result also sometimes ends with a _. Let's remove that too
	if ( result[ result.size() - 1 ] == '_' ) {
		result.pop_back();
	}

	return result;
}

std::string
basename_for_score_function( std::string const & sfxn_name ){
	using namespace utility;
	if ( startswith( sfxn_name, REF_2015 ) ) return REF_2015;
	if ( startswith( sfxn_name, BETA_GENPOT ) ) return BETA_GENPOT;
	if ( startswith( sfxn_name, BETA_JULY15 ) ) return BETA_JULY15;
	if ( startswith( sfxn_name, BETA_NOV15 ) ) return BETA_NOV15;
	if ( startswith( sfxn_name, BETA_NOV16 ) ) return BETA_NOV16;
	if ( startswith( sfxn_name, CENTROID_WTS ) ) return CENTROID_WTS;
	if ( startswith( sfxn_name, DNA_INT_WTS ) ) return DNA_INT_WTS;
	if ( startswith( sfxn_name, DNA_INT_WTS_GB ) ) return DNA_INT_WTS_GB;
	if ( startswith( sfxn_name, DOCK_LOW_PATCH ) ) return DOCK_LOW_PATCH;
	if ( startswith( sfxn_name, DOCK_PATCH ) ) return DOCK_PATCH;
	if ( startswith( sfxn_name, MEMB_HIGHRES_WTS ) ) return MEMB_HIGHRES_WTS;
	if ( startswith( sfxn_name, MM_STD_WTS ) ) return MM_STD_WTS;
	if ( startswith( sfxn_name, PRE_TALARIS_2013_STANDARD_WTS ) ) return PRE_TALARIS_2013_STANDARD_WTS;
	if ( startswith( sfxn_name, RNA_HIRES_WTS ) ) return RNA_HIRES_WTS;
	if ( startswith( sfxn_name, RNA_LORES_PLUS_HIRES_WTS ) ) return RNA_LORES_PLUS_HIRES_WTS;
	if ( startswith( sfxn_name, RNA_LORES_WTS ) ) return RNA_LORES_WTS;
	if ( startswith( sfxn_name, SCORE12_PATCH ) ) return SCORE12_PATCH;
	if ( startswith( sfxn_name, SCORE13 ) ) return SCORE13;
	if ( startswith( sfxn_name, SCORE4_SMOOTH_CART ) ) return SCORE4_SMOOTH_CART;
	if ( startswith( sfxn_name, SOFT_REP_DESIGN_WTS ) ) return SOFT_REP_DESIGN_WTS;
	if ( startswith( sfxn_name, SOFT_REP_WTS ) ) return SOFT_REP_WTS;
	if ( startswith( sfxn_name, TALARIS_2013 ) ) return TALARIS_2013;
	if ( startswith( sfxn_name, TALARIS_2014 ) ) return TALARIS_2014;

	TR << "core::scoring::basename_for_score_function found no match for " << sfxn_name << std::endl;
	return "";
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
