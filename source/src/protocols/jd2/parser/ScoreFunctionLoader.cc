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

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/mistakes.OptionKeys.gen.hh>
#include <basic/datacache/DataMap.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>

// Boost Headers
#include <boost/foreach.hpp>

namespace protocols {
namespace jd2 {
namespace parser {

static THREAD_LOCAL basic::Tracer TR( "protocols.jd2.parser.ScoreFunctionLoader" );

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

	BOOST_FOREACH ( TagCOP scorefxn_tag, scorefxn_tags ) {
		using namespace core::scoring;
		using namespace core::scoring::symmetry;

		ScoreFunctionOP in_scorefxn;
		std::string scorefxn_name( scorefxn_tag->getName() );
		if ( scorefxn_name == "ScoreFunction" ) {
			scorefxn_name = scorefxn_tag->getOption< std::string >( "name" );
		}
		std::string const default_sfxn_weights(
			basic::options::option[ basic::options::OptionKeys::mistakes::restore_pre_talaris_2013_behavior ] ?
			"pre_talaris_2013_standard.wts" : "talaris2013" );
		std::string const scorefxn_weights( scorefxn_tag->getOption<std::string>( "weights", default_sfxn_weights ) );
		if (  scorefxn_tag->hasOption( "weights" ) && scorefxn_tag->hasOption( "patch" ) ) {
			std::string const scorefxn_patch( scorefxn_tag->getOption<std::string>( "patch" ) );
			in_scorefxn = ScoreFunctionFactory::create_score_function( scorefxn_weights, scorefxn_patch);
			TR << "defined score function \"" << scorefxn_name << "\" with weights \""
				<< scorefxn_weights << "\" and patch \"" << scorefxn_patch << "\"\n";
		} else if ( scorefxn_tag->hasOption( "weights" ) ) {
			in_scorefxn = ScoreFunctionFactory::create_score_function( scorefxn_weights );
			TR << "defined score function \"" << scorefxn_name << "\" with weights \""
				<< scorefxn_weights << "\"\n";
		} else {
			in_scorefxn = ScoreFunctionOP( new ScoreFunction );
			in_scorefxn->reset();
			TR << "***WARNING***: No weights/patch defined. Defining " << scorefxn_name << " with all-zero weights.\n";
		}
		BOOST_FOREACH ( TagCOP mod_tag, scorefxn_tag->getTags() ) {
			std::string const tagname( mod_tag->getName() ); //Get the name of the tag
			if ( tagname == "Reweight" ) {
				std::string const scoretype_name( mod_tag->getOption<std::string>( "scoretype" ) );
				core::Real const weight( mod_tag->getOption<core::Real>( "weight" ) );
				TR<<"setting "<<scorefxn_name<<" weight " << scoretype_name << " to " << weight<<'\n';
				core::scoring::ScoreType const type = score_type_from_name( scoretype_name );
				in_scorefxn->set_weight( type, weight );
			} else if ( tagname == "Set" ) { // Set energy method options:
				core::scoring::methods::EnergyMethodOptions emoptions( in_scorefxn->energy_method_options() );
				emoptions.hbond_options().parse_my_tag(mod_tag);
				emoptions.etable_options().parse_my_tag(mod_tag);

				// Set up the list of aa_composition score term setup files:
				if ( mod_tag->hasOption("aa_composition_setup_file") ) {
					std::stringstream filelist( mod_tag->getOption<std::string>("aa_composition_setup_file") );
					utility::vector1 < std::string > filevect;
					while ( !filelist.eof() ) {
						std::string tempstring("");
						filelist >> tempstring;
						if ( !filelist.fail() ) filevect.push_back( tempstring );
					}
					emoptions.append_aa_composition_setup_files( filevect );
				}

				if ( mod_tag->hasOption( "softrep_etable" ) ) {
					if ( mod_tag->getOption<bool>( "softrep_etable" ) ) {
						emoptions.etable_type( core::scoring::FA_STANDARD_SOFT );

					}
				}

				if ( mod_tag->hasOption( "fa_elec_min_dis" ) ) {
					emoptions.elec_min_dis( mod_tag->getOption<core::Real>( "fa_elec_min_dis" ) );
				}
				if ( mod_tag->hasOption( "fa_elec_max_dis" ) ) {
					emoptions.elec_max_dis( mod_tag->getOption<core::Real>( "fa_elec_max_dis" ) );
				}
				if ( mod_tag->hasOption( "fa_elec_dielectric" ) ) {
					emoptions.elec_die( mod_tag->getOption<core::Real>( "fa_elec_dielectric" ) );
				}
				if ( mod_tag->hasOption( "fa_elec_no_dis_dep_die" ) ) {
					emoptions.elec_no_dis_dep_die( mod_tag->getOption<bool>( "fa_elec_no_dis_dep_die" ) );
				}
				if ( mod_tag->hasOption( "exclude_protein_protein_fa_elec" ) ) {
					emoptions.exclude_protein_protein_fa_elec( mod_tag->getOption<bool>( "exclude_protein_protein_fa_elec" ) );
				}
				if ( mod_tag->hasOption( "exclude_DNA_DNA" ) ) {
					emoptions.exclude_DNA_DNA( mod_tag->getOption<bool>( "exclude_DNA_DNA" ) );
				}

				if ( mod_tag->hasOption( "pb_bound_tag" ) ) {
					emoptions.pb_bound_tag( mod_tag->getOption<std::string>("pb_bound_tag" ) );
					TR << "User defined bound tag: " << emoptions.pb_bound_tag() << std::endl;
				}
				if ( mod_tag->hasOption( "pb_unbound_tag" ) ) {
					emoptions.pb_unbound_tag( mod_tag->getOption<std::string>("pb_unbound_tag" ) );
					TR << "User defined unbound tag: " << emoptions.pb_unbound_tag() << std::endl;
				}
				if ( mod_tag->hasOption( "scale_sc_dens" ) ) {
					core::Real scale_sc_dens = mod_tag->getOption<core::Real>("scale_sc_dens" );
					emoptions.set_density_sc_scale_byres( scale_sc_dens );
					TR << "User defined sidechain density reweighing: " << scale_sc_dens << std::endl;
				}
				if ( mod_tag->hasOption( "scale_sc_dens_byres" ) ) {
					utility::vector1< std::string > scale_sc_dens_byres
						= utility::string_split_multi_delim( mod_tag->getOption<std::string>("scale_sc_dens_byres" ), " ,");

					for ( int i=1; i<=(int)scale_sc_dens_byres.size(); ++i ) {
						utility::vector1< std::string > tag = utility::string_split( scale_sc_dens_byres[i] , ':');
						core::chemical::AA aa_i = core::chemical::aa_from_oneletter_code( tag[1][0] );
						core::Real value = std::atof( tag[2].c_str() );
						if ( aa_i<=core::chemical::num_canonical_aas ) {
							emoptions.set_density_sc_scale_byres( aa_i, value );
						}
					}
					TR << "User defined per-residue sidechain density reweighing: " << std::endl;
					for ( int i=1; i<=(int)core::chemical::num_canonical_aas; ++i ) {
						TR << "   " << (core::chemical::AA)i << " " << emoptions.get_density_sc_scale_byres()[i] << std::endl;
					}
				}
				in_scorefxn->set_energy_method_options( emoptions );
			} // tagname == ?
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
		if ( scorefxn_symm ) {
			in_scorefxn = core::scoring::symmetry::symmetrize_scorefunction( *in_scorefxn );
			TR<<"symmetrizing "<<scorefxn_name<<'\n';
		}

		// auto-generate and set cache-tags for bound and unbound energy states, if PB term is used.
		if ( !in_scorefxn->has_zero_weight(PB_elec) ) {

			core::scoring::methods::EnergyMethodOptions emoptions( in_scorefxn->energy_method_options() );
			// Don't overwrite if it's already set, by "Set" modifier.
			if ( emoptions.pb_bound_tag() == "" ) {
				//std::string bound_tag = scorefxn_name + "_" + "bound";
				emoptions.pb_bound_tag( "bound" );
			}
			if ( emoptions.pb_unbound_tag() == "" ) {
				//std::string unbound_tag = scorefxn_name + "_" + "unbound";
				emoptions.pb_unbound_tag( "unbound" );
			}
			in_scorefxn->set_energy_method_options( emoptions );
		}

		bool const data_add_status = data.add( "scorefxns" , scorefxn_name, in_scorefxn );

		if ( !data_add_status ) {
			utility_exit_with_message( "scorefxn " + scorefxn_name + " already exists in the basic::datacache::DataMap, possibly as a default scorefxn. Please rename." );
		}

	}//end user-defined scorefxns
	TR.flush();
}

std::string score_function_subtag_complex_type_namer( std::string const & element_name ) { return "sfxn_" + element_name + "_type"; }
std::string score_function_subtag_group() { return "sfxn_subtag"; }
std::string score_function_tag_group() { return "sfxn_tag"; }

std::string
ScoreFunctionLoader::loader_name() { return "SCOREFXNS"; }

std::string
ScoreFunctionLoader::score_function_loader_ct_namer( std::string const & element_name )
{
	return "sfxn_loader_" + element_name + "_type";
}

void
ScoreFunctionLoader::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	// "Reweight" tag
	AttributeList reweight_attributes;
	reweight_attributes
		+ XMLSchemaAttribute::required_attribute( "scoretype", xs_string )
		+ XMLSchemaAttribute( "weight", xs_decimal );

	// "Set" tag
	AttributeList set_attributes;
	core::scoring::hbonds::HBondOptions::append_schema_attributes( set_attributes );
	core::scoring::etable::EtableOptions::append_schema_attributes( set_attributes );
	set_attributes
		+ XMLSchemaAttribute( "aa_composition_setup_file", xs_string )
		+ XMLSchemaAttribute( "softrep_etable", xsct_rosetta_bool )
		+ XMLSchemaAttribute( "fa_elec_min_dis", xs_decimal )
		+ XMLSchemaAttribute( "fa_elec_max_dis", xs_decimal )
		+ XMLSchemaAttribute( "fa_elec_dielectric", xs_decimal )
		+ XMLSchemaAttribute( "fa_elec_no_dis_dep_die", xsct_rosetta_bool )
		+ XMLSchemaAttribute( "exclude_protein_protein_fa_elec", xsct_rosetta_bool )
		+ XMLSchemaAttribute( "exclude_DNA_DNA", xsct_rosetta_bool )
		+ XMLSchemaAttribute( "pb_bound_tag", xs_string )
		+ XMLSchemaAttribute( "pb_unbound_tag", xs_string )
		+ XMLSchemaAttribute( "scale_sc_dens", xs_decimal )
		+ XMLSchemaAttribute( "scale_sc_dens_byres", xs_decimal );

	// ScoreFunction complex type
	AttributeList scorefunction_ct_attributes;
	scorefunction_ct_attributes
		+ XMLSchemaAttribute::required_attribute( "name", xs_string )
		+ XMLSchemaAttribute( "weights", xs_string )
		+ XMLSchemaAttribute( "patch", xs_string )
		+ XMLSchemaAttribute( "symmetric", xsct_rosetta_bool );
	XMLSchemaSimpleSubelementList subelements;
	subelements
		.add_simple_subelement( "Reweight", reweight_attributes )
		.add_simple_subelement( "Set", set_attributes );
	XMLSchemaComplexTypeGenerator scorefunction_ct_gen;
	scorefunction_ct_gen
		.element_name( "ScoreFunction" )
		.complex_type_naming_func( & score_function_loader_ct_namer )
		.add_attributes( scorefunction_ct_attributes )
		.set_subelements_repeatable( subelements )
		.write_complex_type_to_schema( xsd );

	// SCOREFXNS complex type
	XMLSchemaSimpleSubelementList loader_subelements;
	loader_subelements.add_already_defined_subelement( "ScoreFunction", & score_function_loader_ct_namer );
	XMLSchemaComplexTypeGenerator loader_ct_gen;
	loader_ct_gen
		.element_name( loader_name() )
		.complex_type_naming_func( & score_function_loader_ct_namer )
		.set_subelements_repeatable( loader_subelements )
		.write_complex_type_to_schema( xsd );
}

DataLoaderOP
ScoreFunctionLoaderCreator::create_loader() const { return DataLoaderOP( new ScoreFunctionLoader ); }

std::string
ScoreFunctionLoaderCreator::keyname() const { return ScoreFunctionLoader::loader_name(); }

ScoreFunctionLoaderCreator::DerivedNameFunction
ScoreFunctionLoaderCreator::schema_ct_naming_function() const
{
	return & ScoreFunctionLoader::score_function_loader_ct_namer;
}

} //namespace parser
} //namespace jd2
} //namespace protocols
