// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/qsar/scoring_grid/SiteGrid.cc
/// @author Sam DeLuca

#include <protocols/qsar/scoring_grid/SiteGrid.hh>
#include <protocols/qsar/scoring_grid/SiteGridCreator.hh>

#include <protocols/qsar/scoring_grid/schema_util.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/UltraLightResidue.hh>
#include <core/chemical/AtomType.hh>
#include <core/id/AtomID.hh>

#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tools/make_vector.hh>
#include <utility/json_spirit/json_spirit_value.h>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/io/mpistream.hh>

#include <utility/string_util.hh>

namespace protocols {
namespace qsar {
namespace scoring_grid {


std::string SiteGridCreator::keyname() const
{
	return SiteGrid::grid_name();
}

GridBaseOP SiteGridCreator::create_grid(utility::tag::TagCOP tag) const
{
	GridBaseOP site_grid( new SiteGrid() );

	site_grid->parse_my_tag(tag);

	return site_grid;
}

GridBaseOP SiteGridCreator::create_grid() const
{
	return GridBaseOP( new SiteGrid() );
}

void SiteGridCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SiteGrid::provide_xml_schema( xsd );
}

//std::string SiteGridCreator::grid_name()
//{
// return "SiteGrid";
//}

GridBaseOP SiteGrid::clone() const {
	return GridBaseOP( new SiteGrid( *this ) );
}

SiteGrid::SiteGrid() : SingleGrid("SiteGrid")
{
	atr_weight_ = -1;
	rep_weight_ = 1;
	atr_radius_ = 4.75;
	rep_radius_ = 2.25;
	sites_.clear();
}


SiteGrid::~SiteGrid() = default;

utility::json_spirit::Value SiteGrid::serialize() const
{
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;

	Pair atr("atr",Value(atr_radius_));
	Pair rep("rep",Value(rep_radius_));
	Pair base_data("base_data",SingleGrid::serialize());

	return Value(utility::tools::make_vector(atr,rep,base_data));

}

void SiteGrid::deserialize(utility::json_spirit::mObject data)
{
	atr_radius_ = data["atr"].get_real();
	rep_radius_ = data["rep"].get_real();
	SingleGrid::deserialize(data["base_data"].get_obj());
}

void SiteGrid::parse_my_tag(utility::tag::TagCOP const tag)
{
	//Set residues to emphasize
	std::string const all_residues_str = tag->getOption<std::string>("residues");
	for ( auto residue : utility::string_split(all_residues_str, ',') ) {
		sites_.push_back(utility::string2Size(residue));
	}

	if ( tag->hasOption("atr") ) {
		atr_weight_= tag->getOption<core::Real>("atr");
	}

	if ( tag->hasOption("rep") ) {
		rep_weight_= tag->getOption<core::Real>("rep");
	}

}

void SiteGrid::refresh(core::pose::Pose const & pose, core::Vector const & )
{

	for ( auto residue_number : sites_ ) {
		core::conformation::Residue const & residue(pose.residue(residue_number));

		if ( !residue.is_protein() ) {
			continue;
		}

		for ( core::Size atom_index = 1; atom_index <= residue.nheavyatoms(); ++atom_index ) {

			//set attraction
			set_sphere(residue.xyz(atom_index),atr_radius_,atr_weight_);

			//set neutral
			for ( core::Size atom_index = 1; atom_index <= residue.nheavyatoms(); ++atom_index ) {
				set_sphere(residue.xyz(atom_index),rep_radius_,0);
			}

			//set repulsive
			if ( residue.has("CB") ) set_sphere(residue.xyz("CB"), rep_radius_, rep_weight_);
			if ( residue.has("N") ) set_sphere(residue.xyz("N"), rep_radius_, rep_weight_);
			if ( residue.has("CA") ) set_sphere(residue.xyz("CA"), rep_radius_, rep_weight_);
			if ( residue.has("C") ) set_sphere(residue.xyz("C"), rep_radius_, rep_weight_);
			if ( residue.has("O") ) set_sphere(residue.xyz("O"), rep_radius_, rep_weight_);
		}

	}

}

void SiteGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, core::Size const & )
{
	refresh(pose,center);
}

void SiteGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, utility::vector1<core::Size> )
{
	refresh(pose,center);
}

std::string SiteGrid::grid_name()
{
	return "SiteGrid";
}

void SiteGrid::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute( "grid_name", xs_string, "The name used to insert the scoring grid into the GridManager" )
		+ XMLSchemaAttribute::attribute_w_default("residues", xs_string, "Residue numbers to add additional attractive or repulsive term", "")
		+ XMLSchemaAttribute::attribute_w_default("atr", xsct_real, "Score for attractive term of grid (negative scoring is better)", "-1.0")
		+ XMLSchemaAttribute::attribute_w_default("rep", xsct_real, "Score for repulsive term of grid (negative scoring is better)", "1.0");


	xsd_type_definition_w_attributes( xsd, grid_name(), "A scoring grid that treats selected residue atoms as both attractive within 4.75A (getting a default score of -1) and repulsive within 2.25A (getting a default score of +1) atr and rep changes scoring weights; ", attributes );

}

std::string
SiteGrid::hash_fingerprint() const {
	std::stringstream ss;
	const char sep('\t');
	ss << grid_name();
	ss << sep << get_type(); // Only thing of interest from parent class
	ss << sep << atr_radius_ << sep << rep_radius_;
	ss << sep << atr_weight_ << sep << rep_weight_;

	for ( auto const & site : sites_ ) {
		ss << site;
	}
	return ss.str();
}

}
}
}
