// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/qsar/scoring_grid/ClassicGrid.cc
/// @author Sam DeLuca

#include <protocols/qsar/scoring_grid/ClassicGrid.hh>
#include <protocols/qsar/scoring_grid/ClassicGridCreator.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

#include <utility/tools/make_vector.hh>
#include <utility/json_spirit/json_spirit_value.h>
#include <utility/tag/Tag.hh>

namespace protocols {
namespace qsar {
namespace scoring_grid {

std::string ClassicGridCreator::keyname() const
{
	return ClassicGridCreator::grid_name();
}

GridBaseOP ClassicGridCreator::create_grid(utility::tag::TagCOP tag) const
{
	GridBaseOP classic_grid( new ClassicGrid() );
	classic_grid->parse_my_tag(tag);
	return classic_grid;
}

GridBaseOP ClassicGridCreator::create_grid() const
{
	return GridBaseOP( new ClassicGrid() );
}

std::string ClassicGridCreator::grid_name()
{
	return "ClassicGrid";
}

ClassicGrid::ClassicGrid(): SingleGrid("ClassicGrid"), atr_radius_(4.75), rep_radius_(2.25)
{
	//
}

utility::json_spirit::Value ClassicGrid::serialize()
{
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;

	Pair atr("atr",Value(atr_radius_));
	Pair rep("rep",Value(rep_radius_));
	Pair base_data("base_data",SingleGrid::serialize());

	return Value(utility::tools::make_vector(atr,rep,base_data));

}

void ClassicGrid::deserialize(utility::json_spirit::mObject data)
{
	atr_radius_ = data["atr"].get_real();
	rep_radius_ = data["rep"].get_real();
	SingleGrid::deserialize(data["base_data"].get_obj());
}

void ClassicGrid::parse_my_tag(utility::tag::TagCOP const /*tag*/)
{

}

void ClassicGrid::refresh(core::pose::Pose const & pose, core::Vector const & )
{
	//set attractive zones
	for ( core::Size residue_index = 1; residue_index <=pose.total_residue(); ++residue_index ) {
		core::conformation::Residue const & residue(pose.residue(residue_index));
		if ( !residue.is_protein() ) continue;
		for ( core::Size atom_index = 1; atom_index <= residue.nheavyatoms(); ++atom_index ) {
			set_sphere(residue.xyz(atom_index),atr_radius_,-1);
		}
	}

	//set neutral zones
	for ( core::Size residue_index = 1; residue_index <= pose.total_residue(); ++residue_index ) {
		core::conformation::Residue const & residue(pose.residue(residue_index));
		if ( !residue.is_protein() ) continue;
		for ( core::Size atom_index = 1; atom_index <= residue.nheavyatoms(); ++atom_index ) {
			set_sphere(residue.xyz(atom_index),rep_radius_,0);
		}
	}

	//set repulsive zones
	for ( core::Size residue_index = 1; residue_index <= pose.total_residue(); ++residue_index ) {
		core::conformation::Residue const & residue = pose.residue(residue_index);
		if ( !residue.is_protein() ) continue;
		if ( residue.has("CB") ) set_sphere(residue.xyz("CB"), rep_radius_, 1);
		if ( residue.has("N") ) set_sphere(residue.xyz("N"), rep_radius_, 1);
		if ( residue.has("CA") ) set_sphere(residue.xyz("CA"), rep_radius_, 1);
		if ( residue.has("C") ) set_sphere(residue.xyz("C"), rep_radius_, 1);
		if ( residue.has("O") ) set_sphere(residue.xyz("O"), rep_radius_, 1);
	}
}

void ClassicGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, core::Size const & )
{
	refresh(pose,center);
}

void ClassicGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, utility::vector1<core::Size> )
{
	refresh(pose,center);
}


}
}
}


