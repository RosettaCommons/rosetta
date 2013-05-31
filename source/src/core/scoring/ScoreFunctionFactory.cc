// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/conformation/ResidueFactory.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/DockingScoreFunction.hh>
#include <core/scoring/MinScoreScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/chemical/AA.hh>

// AUTO-REMOVED #include <basic/database/open.hh>

#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/after_opts.hh>

#include <basic/Tracer.hh>

// AUTO-REMOVED #include <utility/io/izstream.hh>
#include <utility/vector1.hh>
// option key includes

#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/mistakes.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <sstream>

using basic::T;

static basic::Tracer tr("core.scoring.ScoreFunctionFactory");

namespace core {
namespace scoring {

ScoreFunctionOP
ScoreFunctionFactory::create_score_function( std::string weights_tag )
{
	utility::vector1< std::string > patch_tags;
	return create_score_function( weights_tag, patch_tags );
	/* removing duplicated code OL 5/24/2012 */
	/*
	ScoreFunctionOP scorefxn( new ScoreFunction );

	load_weights_file( weights_tag, scorefxn );

	// allow user to change weights via options system
	apply_user_defined_reweighting_( scorefxn );

	if ( basic::options::option[ basic::options::OptionKeys::symmetry::symmetry_definition ].user() )	{
		scorefxn = new SymmetricScoreFunction( scorefxn );
	}
	if ( basic::options::option[ basic::options::OptionKeys::score::docking_interface_score ]() ) {
		scorefxn = new DockingScoreFunction( scorefxn );
	}

	return scorefxn;
	*/
}


ScoreFunctionOP
ScoreFunctionFactory::create_score_function( std::string weights_tag, utility::vector1< std::string > patch_tags ) {
	using basic::options::option;
	using namespace basic::options::OptionKeys;

	// create a new scorefunction
	ScoreFunctionOP scorefxn( new ScoreFunction );

	/// Avoid loading the score12 patch if we're using
	/// 1) standard weights,
	/// 2) the score12 patch, and
	/// 3) the flag "score12prime"
	if ( weights_tag == PRE_TALARIS_2013_STANDARD_WTS &&
			basic::options::option[ basic::options::OptionKeys::corrections::score::score12prime ] ) {
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

	load_weights_file( weights_tag, scorefxn );

	for ( utility::vector1< std::string >::const_iterator it = patch_tags.begin(); it != patch_tags.end(); ++it ) {
		std::string const& patch_tag( *it );
		if ( patch_tag.size() && patch_tag != "NOPATCH" ) {
			//			tr.Debug << "SCOREFUNCTION: apply patch "  << patch_tag << std::endl;
			scorefxn->apply_patch_from_file( patch_tag );
		}
	}

	// allow user to change weights via options system
	apply_user_defined_reweighting_( scorefxn );

	if ( basic::options::option[ basic::options::OptionKeys::symmetry::symmetry_definition ].user() )	{
		scorefxn = new SymmetricScoreFunction( scorefxn );
	}
	if ( basic::options::option[ basic::options::OptionKeys::score::docking_interface_score ]() ) {
		scorefxn = new DockingScoreFunction( scorefxn );
	}
	if ( basic::options::option[ basic::options::OptionKeys::score::min_score_score ].user() ) {
		scorefxn = new MinScoreScoreFunction( scorefxn, basic::options::option[ basic::options::OptionKeys::score::min_score_score ]() );
	}
	scorefxn->name( weights_tag );
	return scorefxn;
}


ScoreFunctionOP
ScoreFunctionFactory::create_score_function( std::string weights_tag, std::string const & patch_tag )
{
	using basic::options::option;
	using namespace basic::options::OptionKeys;

	// create a new scorefunction
	ScoreFunctionOP scorefxn( new ScoreFunction );
	utility::vector1< std::string > patch_tags;
	patch_tags.push_back( patch_tag );
	return create_score_function( weights_tag, patch_tags );
	/* REMOVING DUPLICATED CODE OL 5/24/2012 */
	/*	std::string patch_tag_local( patch_tag );
	if ( weights_tag == STANDARD_WTS && patch_tag == SCORE12_PATCH &&
			basic::options::option[ basic::options::OptionKeys::corrections::score::score12prime ] ) {
		weights_tag = "score12prime";
		patch_tag_local = "";
	}

	load_weights_file( weights_tag, scorefxn );

	if ( patch_tag_local.size() && patch_tag_local != "NOPATCH" ) {
		scorefxn->apply_patch_from_file( patch_tag_local );
	}

	// allow user to change weights via options system
	apply_user_defined_reweighting_( scorefxn );


	if ( basic::options::option[ basic::options::OptionKeys::symmetry::symmetry_definition ].user() )	{
		scorefxn = new SymmetricScoreFunction( scorefxn );
	}
	if ( basic::options::option[ basic::options::OptionKeys::score::docking_interface_score ]() ) {
		scorefxn = new DockingScoreFunction( scorefxn );
	}

	return scorefxn;
	*/
}

void ScoreFunctionFactory::apply_user_defined_reweighting_( core::scoring::ScoreFunctionOP scorefxn ) {
	// do some reweighting here. This code could be much more simple if the options system could
	// produce a std::pair< std::string, core::Real >. For now these are separate options for
	// the different reweights.
	using basic::options::option;
	using namespace basic::options::OptionKeys;


	/// new mechanism: set multiple weights using string vector option
	if ( option[ score::set_weights ].user() ) {
		std::string const errmsg("proper format for -set_weights is a list of paired strings, e.g: '-set_weights fa_atr 0.6 -fa_rep 0.55 -fa_sol 0.9' ");
		utility::vector1< std::string > const settings( option[ score::set_weights ]() );
		if ( settings.size()%2 != 0 ) utility_exit_with_message( errmsg );
		for ( Size i=0; i< settings.size()/2; ++i ) {
			if ( !ObjexxFCL::is_float( settings[ 2*i+2 ] ) ) utility_exit_with_message( errmsg );
			ScoreType const t( score_type_from_name( settings[ 2*i + 1] ) );
			Real const value( ObjexxFCL::float_of( settings[ 2*i + 2 ] ) );
			tr << "Setting/modifying scorefxn weight from command line: " << t << ' ' << value << std::endl;
			scorefxn->set_weight( t, value );
		}
	}


	if ( option[ abinitio::rg_reweight ].user() ) {
		scorefxn->set_weight( rg, scorefxn->get_weight( rg ) * option[ abinitio::rg_reweight ]() );
	}
	// offset reference energies using user options, for example: -score:ref_offsets TRP 0.9 HIS 0.3
	if ( option[ score::ref_offsets ].user() ) {

		// get the ref weights from the EnergyMethodOptions object
		methods::EnergyMethodOptions energy_method_options(scorefxn->energy_method_options());
		if (!energy_method_options.has_method_weights(ref)) {
			utility_exit_with_message("option -score:ref_offsets requires preexisting reference energies");
		}
		utility::vector1<core::Real> ref_weights(energy_method_options.method_weights(ref));

		// get the offsets vector and make sure it contains pairs
		utility::vector1<std::string> const & ref_offsets( option[ score::ref_offsets ]() );
		if (ref_offsets.size() % 2 != 0) {
			utility_exit_with_message("option -score:ref_offsets requires pairs of 3 character residue types and offsets");
		}

		// iterate over all pairs
		for (utility::vector1<std::string>::const_iterator iter(ref_offsets.begin()), iter_end(ref_offsets.end());
		     iter != iter_end; ++iter) {
			// get the aa type from the pair
			std::istringstream aa_iss(*iter);
			core::chemical::AA aa;
			if (!(aa_iss >> aa)) {
				utility_exit_with_message(aa_iss.str()+" is not a valid 3 character residue type for -score:ref_offsets");
			}
			// get the offset from the pair
			std::istringstream offset_iss(*(++iter));
			core::Real offset;
			if (!(offset_iss >> offset)) {
				utility_exit_with_message(offset_iss.str()+" is not a valid offset for -score:ref_offsets");
			}
			// offset the weight
			ref_weights[aa] += offset;
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

std::string const TALARIS_2013( "talaris2013" );
std::string const PRE_TALARIS_2013_STANDARD_WTS( "pre_talaris_2013_standard" );
std::string const CENTROID_WTS( "cen_std" );
std::string const SOFT_REP_WTS( "soft_rep" );
std::string const SOFT_REP_DESIGN_WTS( "soft_rep_design" );
std::string const DNA_INT_WTS( "dna_no_gb" );
std::string const DNA_INT_WTS_GB( "dna" );
std::string const MM_STD_WTS( "mm_std" );
std::string const RNA_LORES_WTS( "rna_lores" );
std::string const RNA_HIRES_WTS( "rna_hires" );
std::string const RNA_LORES_PLUS_HIRES_WTS( "rna_lores_plus_hires" );
std::string const MEMB_HIGHRES_WTS( "membrane_highres" ); //pba

std::string const SCORE12_PATCH( "score12" );
std::string const SCORE13( "score13" );
std::string const DOCK_PATCH( "docking" );
std::string const DOCK_LOW_PATCH( "docking_cen" );


core::scoring::ScoreFunctionOP getScoreFunction( bool const is_fullatom /* default true */ ) {
	using basic::options::option;
	using namespace basic::options::OptionKeys;

	core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction() );

	//if( option[ score::empty ]() || option[ abinitio::membrane ]() /*fullatom not implemented for membrane yet */) return scorefxn;

	if( option[ score::empty ]() ) return scorefxn;

	std::string weight_set = option[ score::weights ];
	utility::vector1< std::string > patch_tags = option[ score::patch ]();

	if ( !option[ score::weights ].user() && !is_fullatom ){

		// Defalt score of centroid is cen_wts when is_fullatom is false and user has not specified a score weights
		weight_set = CENTROID_WTS;

	}else{

		/// Default score is score12 if the user has not specified a score weights file or a patch file
		/// on the command line.  If the user has specified that they would like the standard weight set,
		/// and has not also asked for the score12 patch, then do not apply the score12 patch to it.

		//tr << "getScoreFunction1: weight set " << weight_set << " same ? " << ( weight_set == "pre_talaris_2013_standard.wts") << std::endl;
		//tr << "option[ score::weights ].user() ? " << option[ score::weights ].user() << std::endl;
		//tr << "option[ score::patch ].user() ? " << option[ score::patch ].user() << std::endl;

		if ( ( weight_set == "pre_talaris_2013_standard.wts" && !option[ score::weights ].user() ) &&
				 ( !option[ score::patch ].user() ) ) {
			patch_tags.push_back( "score12" );
			//tr << "pushing back score12 patch" << std::endl;
			if( basic::options::option[basic::options::OptionKeys::corrections::correct]) {
				//tr << "setting weight set to score12_w_corrections" << std::endl;
				weight_set = "score12_w_corrections";
				patch_tags.clear();
			} else if ( basic::options::option[ basic::options::OptionKeys::corrections::hbond_sp2_correction ] ){
				patch_tags.clear();
				if( basic::options::option[ basic::options::OptionKeys::corrections::score::dun10 ] ){
					weight_set = "sp2_correction";
				} else {
					weight_set = "sp2_correction_dun02";
				}
			} else if ( basic::options::option[ basic::options::OptionKeys::corrections::score::score12prime ] ) {
				weight_set = "score12prime";
				patch_tags.clear();
			}
		}
		//tr << "getScoreFunction2: weight set " << weight_set << std::endl;

	}

	T("core.scoring.ScoreFunctionFactory") << "SCOREFUNCTION: " << weight_set << std::endl;
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

	// add in constraint weights if specified by the user. maybe we want other constraint
	// types to be weighted as well ...
	//
	// mtyka: No No No we dont want this here. Add constraints and weights ouside this function.
	//  just setting the weights alone isnt gonna get oyu far anyway.
	//if ( option[ constraints::cst_weight ].user() ) {
	//	scorefxn->set_weight( atom_pair_constraint, option[ constraints::cst_weight ]() );
	//	}

	// this is already done in create_score_function
	//	if ( basic::options::option[ basic::options::OptionKeys::symmetry::symmetry_definition ].user() )
	//  {
	//		return core::scoring::ScoreFunctionOP( new core::scoring::symmetry::SymmetricScoreFunction( scorefxn ) );
	//	}

	return scorefxn;
}

core::scoring::ScoreFunctionOP getScoreFunctionLegacy(
	std::string pre_talaris_2013_weight_set,
	std::string pre_talaris_2013_patch_file
)
{
	if ( basic::options::option[ basic::options::OptionKeys::mistakes::restore_pre_talaris_2013_behavior ] ) {
		return ScoreFunctionFactory::create_score_function( pre_talaris_2013_weight_set, pre_talaris_2013_patch_file );
	}
	return getScoreFunction();
}

std::string
getScoreFunctionName(
	bool const is_fullatom /* default true */
) {
	using basic::options::option;
	using namespace basic::options::OptionKeys;

	if( option[ score::empty ]() ) return "empty";

	std::string weight_set = option[ score::weights ];
	utility::vector1< std::string > patch_tags = option[ score::patch ]();


	if ( !option[ score::weights ].user() && !is_fullatom ){

		// Defalt score of centroid is cen_wts when is_fullatom is false and user has not specified a score weights
		weight_set = CENTROID_WTS;

	}else{

		/// Default score is score12 if the user has not specified a score weights file or a patch file
		/// on the command line.  If the user has specified that they would like the standard weight set,
		/// and has not also asked for the score12 patch, then do not apply the score12 patch to it.
		if ( ( weight_set == "pre_talaris_2013_standard" && !option[ score::weights ].user() ) &&
				 ( !option[ score::patch ].user() ) ) {
			patch_tags.push_back( "score12" );
			if( option[ corrections::correct ]) {
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
	for(Size ii=1; ii <= patch_tags.size(); ++ii){
		if( patch_tags[ii] == "NOPATCH" ) continue;
		patch_string << "_" << patch_tags[ii];
	}

	return weight_set + patch_string.str();
}



} // namespace scoring
} // namespace core


