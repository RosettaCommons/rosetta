// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/parser/DataLoader.cc
/// @brief  Implementation of the XML parser's DataLoader base class (ctor & dstor)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <protocols/jd2/parser/ScoreFunctionLoader.hh>
#include <protocols/jd2/parser/StandardLoaderCreators.hh>

// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/etable/EtableOptions.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/tag/Tag.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/mistakes.OptionKeys.gen.hh>

// Boost Headers
#include <boost/foreach.hpp>

#include <basic/datacache/DataMap.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>

namespace protocols {
namespace jd2 {
namespace parser {

static thread_local basic::Tracer TR( "protocols.jd2.parser.ScoreFunctionLoader" );

ScoreFunctionLoader::ScoreFunctionLoader() {}
ScoreFunctionLoader::~ScoreFunctionLoader() {}

void ScoreFunctionLoader::load_data(
	core::pose::Pose const &,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) const
{
	using namespace utility::tag;
	typedef utility::vector0< TagCOP > TagCOPs;

	TagCOPs const & scorefxn_tags( tag->getTags() );

	BOOST_FOREACH(TagCOP scorefxn_tag, scorefxn_tags){
		using namespace core::scoring;
		using namespace core::scoring::symmetry;

		ScoreFunctionOP in_scorefxn;
		std::string const scorefxn_name( scorefxn_tag->getName() );
		std::string const default_sfxn_weights(
			basic::options::option[ basic::options::OptionKeys::mistakes::restore_pre_talaris_2013_behavior ] ?
			"pre_talaris_2013_standard.wts" : "talaris2013" );
		std::string const scorefxn_weights( scorefxn_tag->getOption<std::string>( "weights", default_sfxn_weights ) );
		if(  scorefxn_tag->hasOption( "weights" ) && scorefxn_tag->hasOption( "patch" ) ) {
			std::string const scorefxn_patch( scorefxn_tag->getOption<std::string>( "patch" ) );
			in_scorefxn = ScoreFunctionFactory::create_score_function( scorefxn_weights, scorefxn_patch);
			TR << "defined score function \"" << scorefxn_name << "\" with weights \""
				<< scorefxn_weights << "\" and patch \"" << scorefxn_patch << "\"\n";
		} else if ( scorefxn_tag->hasOption( "weights" ) ) {
			in_scorefxn = ScoreFunctionFactory::create_score_function( scorefxn_weights );
			TR << "defined score function \"" << scorefxn_name << "\" with weights \""
				<< scorefxn_weights << "\"\n";
		} else {
			in_scorefxn = new ScoreFunction;
			in_scorefxn->reset();
			TR << "***WARNING***: No weights/patch defined. Defining " << scorefxn_name << " with all-zero weights.\n";
		}
		BOOST_FOREACH(TagCOP mod_tag, scorefxn_tag->getTags()){
			if( mod_tag->getName() == "Reweight" ) {
				std::string const scoretype_name( mod_tag->getOption<std::string>( "scoretype" ) );
				core::Real const weight( mod_tag->getOption<core::Real>( "weight" ) );
				TR<<"setting "<<scorefxn_name<<" weight " << scoretype_name << " to " << weight<<'\n';
				core::scoring::ScoreType const type = score_type_from_name( scoretype_name );
				in_scorefxn->set_weight( type, weight );
			}

			// Set energy method options:
			if( mod_tag->getName() == "Set" ){
				core::scoring::methods::EnergyMethodOptions emoptions( in_scorefxn->energy_method_options() );
				emoptions.hbond_options().parse_my_tag(mod_tag);
				emoptions.etable_options().parse_my_tag(mod_tag);

				if( mod_tag->hasOption( "softrep_etable" )) {
					if ( mod_tag->getOption<bool>( "softrep_etable" )) {
						emoptions.etable_type( core::scoring::FA_STANDARD_SOFT );

					}
				}

				if( mod_tag->hasOption( "fa_elec_min_dis" )) {
					emoptions.elec_min_dis( mod_tag->getOption<core::Real>( "fa_elec_min_dis" ) );
				}
				if( mod_tag->hasOption( "fa_elec_max_dis" )) {
					emoptions.elec_max_dis( mod_tag->getOption<core::Real>( "fa_elec_max_dis" ) );
				}
				if( mod_tag->hasOption( "fa_elec_dielectric" )) {
					emoptions.elec_die( mod_tag->getOption<core::Real>( "fa_elec_dielectric" ) );
				}
				if( mod_tag->hasOption( "fa_elec_no_dis_dep_die" )) {
					emoptions.elec_no_dis_dep_die( mod_tag->getOption<bool>( "fa_elec_no_dis_dep_die" ) );
				}
				if( mod_tag->hasOption( "exclude_protein_protein_fa_elec" )) {
					emoptions.exclude_protein_protein_fa_elec( mod_tag->getOption<bool>( "exclude_protein_protein_fa_elec" ) );
				}
				if( mod_tag->hasOption( "exclude_DNA_DNA" )) {
					emoptions.exclude_DNA_DNA( mod_tag->getOption<bool>( "exclude_DNA_DNA" ) );
				}

				if( mod_tag->hasOption( "pb_bound_tag" )) {
					emoptions.pb_bound_tag( mod_tag->getOption<std::string>("pb_bound_tag" ) );
					TR << "User defined bound tag: " << emoptions.pb_bound_tag() << std::endl;
				}
				if( mod_tag->hasOption( "pb_unbound_tag" )) {
					emoptions.pb_unbound_tag( mod_tag->getOption<std::string>("pb_unbound_tag" ) );
					TR << "User defined unbound tag: " << emoptions.pb_unbound_tag() << std::endl;
				}
				in_scorefxn->set_energy_method_options( emoptions );
			}
		} // Mod tags

		// weights for arbitrary ScoreFunctions should be tampered with only as a consequence of user input--NEVER by default

		// hotspot hash constraint
		if ( scorefxn_tag->hasOption("hs_hash") ) {
			core::Real hotspot_hash( 0.0 ); // APL FIX THIS!  This used to be initialized when the HotspotHashingConstraints were read in.
			core::Real const hs_hash( scorefxn_tag->getOption<core::Real>( "hs_hash", hotspot_hash ) );
			TR<<"setting "<<scorefxn_name<<" backbone_stub_constraint to "<<hs_hash<<'\n';
			in_scorefxn->set_weight( backbone_stub_constraint, hs_hash );
		}

		//fpd should we symmetrize scorefunction?
		bool const scorefxn_symm( scorefxn_tag->getOption<bool>( "symmetric", 0 ) );
		if (scorefxn_symm) {
			in_scorefxn = core::scoring::symmetry::symmetrize_scorefunction( *in_scorefxn );
			TR<<"symmetrizing "<<scorefxn_name<<'\n';
		}

		// auto-generate and set cache-tags for bound and unbound energy states, if PB term is used.
		if( !in_scorefxn->has_zero_weight(PB_elec) ) {

			core::scoring::methods::EnergyMethodOptions emoptions( in_scorefxn->energy_method_options() );
			// Don't overwrite if it's already set, by "Set" modifier.
			if( emoptions.pb_bound_tag() == "" ){
				//std::string bound_tag = scorefxn_name + "_" + "bound";
				emoptions.pb_bound_tag( "bound" );
			}
			if( emoptions.pb_unbound_tag() == "" ){
				//std::string unbound_tag = scorefxn_name + "_" + "unbound";
				emoptions.pb_unbound_tag( "unbound" );
			}
			in_scorefxn->set_energy_method_options( emoptions );
		}

		bool const data_add_status = data.add( "scorefxns" , scorefxn_name, in_scorefxn );

		if( !data_add_status )
			utility_exit_with_message( "scorefxn " + scorefxn_name + " already exists in the basic::datacache::DataMap, possibly as a default scorefxn. Please rename." );

	}//end user-defined scorefxns
	TR.flush();
}

DataLoaderOP
ScoreFunctionLoaderCreator::create_loader() const { return new ScoreFunctionLoader; }

std::string
ScoreFunctionLoaderCreator::keyname() const { return "SCOREFXNS"; }


} //namespace parser
} //namespace jd2
} //namespace protocols
