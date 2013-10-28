// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/qsar/scoring_grid/LigandPropertyScore.cc
/// @author Sam DeLuca

#include <protocols/qsar/scoring_grid/LigandPropertyScore.hh>
#include <protocols/qsar/scoring_grid/LigandPropertyScoreCreator.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/UltraLightResidue.hh>

#include <utility/tag/Tag.hh>
#include <utility/tools/make_vector.hh>
namespace protocols {
namespace qsar {
namespace scoring_grid {

GridBaseOP LigandPropertyScoreCreator::create_grid(utility::tag::TagCOP const tag) const
{
	GridBaseOP ligand_property_score = new LigandPropertyScore();
	ligand_property_score->parse_my_tag(tag);
	return ligand_property_score;
}

GridBaseOP LigandPropertyScoreCreator::create_grid() const
{
	return new LigandPropertyScore();
}

std::string LigandPropertyScoreCreator::keyname() const
{
	return grid_name();
}

std::string LigandPropertyScoreCreator::grid_name()
{
	return "LigandPropertyScore";
}

LigandPropertyScore::LigandPropertyScore()
{

}

void LigandPropertyScore::parse_my_tag(utility::tag::TagCOP const tag)
{
	if(!tag->hasOption("parameter"))
	{
		throw utility::excn::EXCN_RosettaScriptsOption("Grid Score term LigandPropertyScore needs tag 'parameter'");
	}

	parameter_tag_ = tag->getOption<std::string>("parameter");

}

core::Real LigandPropertyScore::score(core::conformation::UltraLightResidue const & residue, core::Real const , qsarMapOP )
{
	core::Real property_score = residue.residue()->type().get_numeric_property(parameter_tag_);
	return property_score;
}

core::Real LigandPropertyScore::score(core::conformation::Residue const & residue, core::Real const , qsarMapOP )
{
	core::Real property_score = residue.type().get_numeric_property(parameter_tag_);
	return property_score;
}

std::string LigandPropertyScore::get_type()
{
	return "LigandPropertyScore";
}

utility::json_spirit::Value LigandPropertyScore::serialize()
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

}
}
}
