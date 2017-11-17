// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/qsar/scoring_grid/LigandPropertyScore.cc
/// @author Sam DeLuca

#include <protocols/qsar/scoring_grid/LigandPropertyScore.hh>
#include <protocols/qsar/scoring_grid/LigandPropertyScoreCreator.hh>

#include <protocols/qsar/scoring_grid/schema_util.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/UltraLightResidue.hh>

#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tools/make_vector.hh>


namespace protocols {
namespace qsar {
namespace scoring_grid {

GridBaseOP LigandPropertyScoreCreator::create_grid(utility::tag::TagCOP tag) const
{
	GridBaseOP ligand_property_score( new LigandPropertyScore() );
	ligand_property_score->parse_my_tag(tag);
	return ligand_property_score;
}

GridBaseOP LigandPropertyScoreCreator::create_grid() const
{
	return GridBaseOP( new LigandPropertyScore() );
}

std::string LigandPropertyScoreCreator::keyname() const
{
	return LigandPropertyScore::grid_name();
}

void LigandPropertyScoreCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LigandPropertyScore::provide_xml_schema( xsd );
}

//std::string LigandPropertyScoreCreator::grid_name()
//{
// return "LigandPropertyScore";
//}

LigandPropertyScore::LigandPropertyScore()
{

}

/// @brief Make a copy of the grid, respecting the subclassing.
GridBaseOP LigandPropertyScore::clone() const {
	return GridBaseOP( new LigandPropertyScore( *this ) );
}

void LigandPropertyScore::parse_my_tag(utility::tag::TagCOP tag)
{
	if ( !tag->hasOption("parameter") ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Grid Score term LigandPropertyScore needs tag 'parameter'");
	}

	parameter_tag_ = tag->getOption<std::string>("parameter");

}

core::Real LigandPropertyScore::score(core::conformation::UltraLightResidue const & residue, core::Real const , qsarMapCOP ) const
{
	core::Real property_score = residue.residue()->type().get_numeric_property(parameter_tag_);
	return property_score;
}

core::Real LigandPropertyScore::score(core::conformation::Residue const & residue, core::Real const , qsarMapCOP ) const
{
	core::Real property_score = residue.type().get_numeric_property(parameter_tag_);
	return property_score;
}

std::string LigandPropertyScore::get_type() const
{
	//return "LigandPropertyScore";
	return grid_name();
}

utility::json_spirit::Value LigandPropertyScore::serialize() const
{
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;

	Pair type_record("type",Value(get_type()));
	Pair parameter_record("parameter_tag",Value(parameter_tag_));

	return Value(utility::tools::make_vector(type_record,parameter_record));
}

void LigandPropertyScore::deserialize(utility::json_spirit::mObject data)
{
	parameter_tag_ = data["parameter_tag"].get_str();
}

std::string LigandPropertyScore::grid_name()
{
	return "LigandPropertyScore";
}

void LigandPropertyScore::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::required_attribute( "parameter", xs_string, "The numeric property that the"
		" ResidueType of the incoming Residue holds that the LigandPropertyScore is going to query for" );

	xsd_type_definition_w_attributes( xsd, grid_name(), "A scoring grid that computes a score simply based"
		" on a static property of the ligand's chemical composition. Not exactly a grid -- but still capable"
		" of delivering a bonus or a penalty for a given ResidueType and fits smoothly within the scoring"
		" grid machinery", attributes );

}

/// @brief Print a brief summary about this grid to the provided output stream
void LigandPropertyScore::show( std::ostream & out ) const {
	out << "LigandPropertyScore grid using parameter " << parameter_tag_ << std::endl;
}

std::string
LigandPropertyScore::hash_fingerprint() const {
	std::stringstream ss;
	const char sep('\t');
	ss << grid_name();
	ss << sep << parameter_tag_;
	return ss.str();
}
}
}
}
