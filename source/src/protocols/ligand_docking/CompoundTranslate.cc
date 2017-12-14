// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Gordon Lemmon (glemmon@gmail.com)

// Unit Headers
#include <protocols/ligand_docking/CompoundTranslate.hh>
#include <protocols/ligand_docking/CompoundTranslateCreator.hh>
#include <protocols/ligand_docking/Translate.hh>

#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <numeric/random/random_permutation.hh>

#include <core/kinematics/Jump.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>

#include <utility/excn/Exceptions.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


//Auto Headers
using basic::Error;
using basic::Warning;

namespace protocols {
namespace ligand_docking {

static basic::Tracer translate_tracer( "protocols.ligand_docking.CompoundTranslate", basic::t_debug );


// XRW TEMP std::string
// XRW TEMP CompoundTranslateCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return CompoundTranslate::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP CompoundTranslateCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new CompoundTranslate );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP CompoundTranslate::mover_name()
// XRW TEMP {
// XRW TEMP  return "CompoundTranslate";
// XRW TEMP }

/// @brief
CompoundTranslate::CompoundTranslate():
	//utility::pointer::ReferenceCount(),
	Mover("CompoundTranslate")
{}

CompoundTranslate::CompoundTranslate(CompoundTranslate const & /*that*/) = default;

CompoundTranslate::~CompoundTranslate() = default;

protocols::moves::MoverOP CompoundTranslate::clone() const {
	return protocols::moves::MoverOP( new CompoundTranslate( *this ) );
}

protocols::moves::MoverOP CompoundTranslate::fresh_instance() const {
	return protocols::moves::MoverOP( new CompoundTranslate );
}

// XRW TEMP std::string CompoundTranslate::get_name() const{
// XRW TEMP  return "CompoundTranslate";
// XRW TEMP }

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
CompoundTranslate::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
)
{
	if ( tag->getName() != "CompoundTranslate" ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "This should be impossible");
	}
	if ( ! tag->hasOption("randomize_order") ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "CompoundTranslate needs a 'randomize_order' option");
	}
	if ( ! tag->hasOption("allow_overlap") ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "CompoundTranslate needs an 'allow_overlap' option");
	}

	// parsing randomize_order tag
	bool randomize = tag->getOption< bool >("randomize_order");
	if ( randomize ) {
		randomize_order_= true;
	} else {
		randomize_order_= false;
	}
	// parsing allow_overlap tag
	bool allow_overlap = tag->getOption< bool >("allow_overlap");
	if ( allow_overlap ) {
		allow_overlap_= true;
	} else {
		allow_overlap_= false;
	}

	for ( utility::tag::TagCOP subtag : tag->getTags() ) {
		std::string const & name = subtag->getName();
		if ( name == "Translate" ) {
			TranslateOP translate( new Translate() );
			translate->parse_my_tag( subtag, datamap, filters, movers, pose);
			translates_.push_back(translate);
		} else if ( name == "Translates" ) {
			if ( ! subtag->hasOption("chain") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'Translates' mover requires chain tag");
			if ( ! subtag->hasOption("distribution") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'Translates' mover requires distribution tag");
			if ( ! subtag->hasOption("angstroms") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'Translates' mover requires angstroms tag");
			if ( ! subtag->hasOption("cycles") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'Translates' mover requires cycles tag");

			std::string const & chain = subtag->getOption<std::string>("chain");
			utility::vector1<core::Size> chain_ids = core::pose::get_chain_ids_from_chain(chain, pose);

			for ( core::Size const chain_id : chain_ids ) {
				Translate_info translate_info;
				translate_info.chain_id= chain_id;
				translate_info.jump_id = core::pose::get_jump_id_from_chain_id(chain_id, pose);
				std::string distribution_str= subtag->getOption<std::string>("distribution");
				translate_info.distribution= get_distribution(distribution_str);
				translate_info.angstroms = subtag->getOption<core::Real>("angstroms");
				translate_info.cycles = subtag->getOption<core::Size>("cycles");
				if ( subtag->hasOption("force") ) {
					if ( subtag->getOption< bool >("force") ) {
						translate_info.force= true;
					}
				}
				translates_.push_back(protocols::ligand_docking::TranslateOP( new Translate(translate_info) ));
			}
		} else {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "CompoundTranslate only takes Translate or Translates child tags");
		}
	}
}

void CompoundTranslate::apply(core::pose::Pose & pose) {
	if ( randomize_order_ ) {
		numeric::random::random_permutation(translates_, numeric::random::rg());
	}

	std::set<core::Size> chains_to_translate;

	// TranslateOPs::iterator begin= translates_.begin(); // Unused variable causes warning.
	// TranslateOPs::iterator const end= translates_.end(); // Unused variable causes warning.

	for ( TranslateOP translate : translates_ ) {
		core::Size chain_id= translate->get_chain_id(pose);
		chains_to_translate.insert(chain_id);
	}

	if ( allow_overlap_ ) {
		for ( TranslateOP translate : translates_ ) {
			translate->add_excluded_chains(chains_to_translate.begin(), chains_to_translate.end());
			translate->apply(pose);
		}
	} else { // remove each chain from the exclusion list so that placed chains are in the grid
		for ( TranslateOP translate : translates_ ) {
			translate->add_excluded_chains(chains_to_translate.begin(), chains_to_translate.end());
			translate->apply(pose);
			core::Size chain_id= translate->get_chain_id(pose);
			chains_to_translate.erase(chain_id);
		}
	}
}

std::string CompoundTranslate::get_name() const {
	return mover_name();
}

std::string CompoundTranslate::mover_name() {
	return "CompoundTranslate";
}

void CompoundTranslate::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	// Main attributes

	XMLSchemaRestriction restriction_type;
	restriction_type.name( "distribution_string" );
	restriction_type.base_type( xs_string );
	restriction_type.add_restriction( xsr_pattern, "uniform|gaussian" );
	xsd.add_top_level_element( restriction_type );


	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute("randomize_order", xs_string, "Randomize order of perturbations.")
		+ XMLSchemaAttribute::required_attribute("allow_overlap", xs_string, "XRW TO DO");

	// Subelements
	AttributeList subelement_attributes;
	subelement_attributes
		+ XMLSchemaAttribute::required_attribute("chain", xs_string, "Chain to be translated.")
		+ XMLSchemaAttribute::required_attribute("distribution", "distribution_string", "Distribution from which to sample distances. Valid options are uniform or gaussian")
		+ XMLSchemaAttribute::required_attribute("angstroms", xsct_real, "Uniform translations of up to X angstroms")
		+ XMLSchemaAttribute::required_attribute("cycles", xsct_non_negative_integer, "Number of cycles to run.")
		+ XMLSchemaAttribute("force", xsct_rosetta_bool, "Force neighbor atoms to move if clashing, default=false");

	XMLSchemaSimpleSubelementList subelement_list;

	subelement_list.add_simple_subelement("Translate", subelement_attributes, "Configure a Translation element that can be used together with others.")
		.add_simple_subelement( "Translates", subelement_attributes, "Configure a Translation element that can be used together with others.");

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "Combines multiple individually configurable Translate movers", attlist, subelement_list );
}

std::string CompoundTranslateCreator::keyname() const {
	return CompoundTranslate::mover_name();
}

protocols::moves::MoverOP
CompoundTranslateCreator::create_mover() const {
	return protocols::moves::MoverOP( new CompoundTranslate );
}

void CompoundTranslateCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	CompoundTranslate::provide_xml_schema( xsd );
}


} //namespace ligand_docking
} //namespace protocols
