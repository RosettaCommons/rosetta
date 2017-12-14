// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/ResfileReader.cc
/// @brief  implementation of resfile reader and its command classes
/// @author Gordon Lemmon (glemmon@gmail.com)

// Unit Headers
#include <protocols/ligand_docking/Rotates.hh>
#include <protocols/ligand_docking/RotatesCreator.hh>
#include <protocols/ligand_docking/Rotate.hh>
#include <protocols/ligand_docking/DistributionMap.hh>

#include <basic/datacache/DataMap.hh>
#include <core/pose/util.hh> // includes Pose.hh
#include <core/pose/chains_util.hh>

// Utility Headers
#include <numeric/random/random.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/tag/Tag.hh>

#include <algorithm>

#include <utility/excn/Exceptions.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

using basic::Error;
using basic::Warning;

namespace protocols {
namespace ligand_docking {

static basic::Tracer rotates_tracer( "protocols.ligand_docking.ligand_options.rotates", basic::t_debug );

// XRW TEMP std::string
// XRW TEMP RotatesCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return Rotates::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP RotatesCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new Rotates );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP Rotates::mover_name()
// XRW TEMP {
// XRW TEMP  return "Rotates";
// XRW TEMP }

/// @brief
Rotates::Rotates(): Mover("Rotates")
{}

Rotates::Rotates(Rotates const & that):
	//utility::pointer::ReferenceCount(),
	protocols::moves::Mover( that )
{}

Rotates::~Rotates() = default;

protocols::moves::MoverOP Rotates::clone() const {
	return protocols::moves::MoverOP( new Rotates( *this ) );
}

protocols::moves::MoverOP Rotates::fresh_instance() const {
	return protocols::moves::MoverOP( new Rotates );
}

// XRW TEMP std::string Rotates::get_name() const{
// XRW TEMP  return "Rotates";
// XRW TEMP }

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
Rotates::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data_map*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & pose
)
{
	if ( tag->getName() != "Rotates" ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "This should be impossible");
	}
	if ( ! tag->hasOption("distribution") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'Rotates' mover requires 'distribution' tag");
	if ( ! tag->hasOption("degrees") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'Rotates' mover requires 'degrees' tag");
	if ( ! tag->hasOption("cycles") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'Rotates' mover requires 'cycles' tag");
	if ( tag->hasOption("chain") && tag->hasOption("chains") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'Rotates' mover cannot have both a 'chain' and a 'chains' tag");
	if ( ! (tag->hasOption("chain") || tag->hasOption("chains") ) ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'Rotates' mover requires either a 'chain' or a 'chains' tag");

	utility::vector1<std::string> chain_strs;
	if ( tag->hasOption("chain") ) {
		chain_strs.push_back( tag->getOption<std::string>("chain") );
	} else if ( tag->hasOption("chains") ) {
		std::string const chains_str = tag->getOption<std::string>("chains");
		chain_strs= utility::string_split(chains_str, ',');
	}

	std::string const distribution_str= tag->getOption<std::string>("distribution");
	Distribution distribution= get_distribution(distribution_str);
	auto const degrees = tag->getOption<core::Size>("degrees");
	auto const cycles = tag->getOption<core::Size>("cycles");

	for ( std::string const & chain : chain_strs ) {
		utility::vector1<core::Size> chain_ids = core::pose::get_chain_ids_from_chain(chain, pose);
		for ( core::Size const chain_id : chain_ids ) {
			Rotate_info rotate_info;
			rotate_info.chain_id = chain_id;
			rotate_info.jump_id = core::pose::get_jump_id_from_chain_id(chain_id, pose);
			rotate_info.distribution= distribution;
			rotate_info.degrees = degrees;
			rotate_info.cycles = cycles;
			rotates_.push_back( protocols::ligand_docking::RotateOP( new Rotate(rotate_info) ) );
		}
	}
}

void Rotates::apply(core::pose::Pose & pose){
	for ( RotateOP rotate : rotates_ ) {
		rotate->apply(pose);
	}
}

std::string Rotates::get_name() const {
	return mover_name();
}

std::string Rotates::mover_name() {
	return "Rotates";
}

void Rotates::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	XMLSchemaRestriction restriction_type;
	restriction_type.name( "distribution_string" );
	restriction_type.base_type( xs_string );
	restriction_type.add_restriction( xsr_pattern, "uniform|gaussian" );
	xsd.add_top_level_element( restriction_type );

	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute("distribution", "distribution_string", "Sampling distribution; Either \"uniform\" or \"gaussian\"")
		+ XMLSchemaAttribute::required_attribute("degrees", xsct_non_negative_integer, "How degrees should be rotated around. Recommended=360")
		+ XMLSchemaAttribute::required_attribute("cycles", xsct_non_negative_integer, "Number of cycles. Recommended: 1000")
		+ XMLSchemaAttribute("chain", xs_string, "Chain to be rotated. Not compatible with \"chains\" option.")
		+ XMLSchemaAttribute("chains", xs_string, "Comma-separated list of chain IDs. Not compatible with \"chain\" option.");
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Perform a course random rotation "
		" throughout all rotational degrees of freedom.", attlist );
}

std::string RotatesCreator::keyname() const {
	return Rotates::mover_name();
}

protocols::moves::MoverOP
RotatesCreator::create_mover() const {
	return protocols::moves::MoverOP( new Rotates );
}

void RotatesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	Rotates::provide_xml_schema( xsd );
}



} //namespace ligand_docking
} //namespace protocols
