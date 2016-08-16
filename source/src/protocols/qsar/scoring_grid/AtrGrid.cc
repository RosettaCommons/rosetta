// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/qsar/scoring_grid/AtrGrid.cc
/// @author Sam DeLuca


#include <protocols/qsar/scoring_grid/AtrGrid.hh>
#include <protocols/qsar/scoring_grid/AtrGridCreator.hh>

#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/tools/make_vector.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/json_spirit/json_spirit_value.h>

namespace protocols {
namespace qsar {
namespace scoring_grid {

std::string AtrGridCreator::keyname()const
{
	return AtrGridCreator::grid_name();
}

GridBaseOP AtrGridCreator::create_grid(utility::tag::TagCOP tag) const
{
	GridBaseOP atr_grid( new AtrGrid() );

	atr_grid->parse_my_tag(tag);

	return atr_grid;
}

GridBaseOP AtrGridCreator::create_grid() const
{

	return GridBaseOP( new AtrGrid() );
}


std::string AtrGridCreator::grid_name()
{
	return "AtrGrid";
}


AtrGrid::AtrGrid() :
	SingleGrid("AtrGrid"),
	inner_radius_(2.25),
	outer_radius_(4.75),
	bb_(-1),
	sc_(-1),
	ligand_(-1)
{
	//
}

void
AtrGrid::parse_my_tag(utility::tag::TagCOP tag){

	if ( tag->hasOption("bb") || tag->hasOption("sc") || tag->hasOption("ligand") ) {
		// the user MUST provide all 3 if he/she is providing any of these 3 options
		if ( !(tag->hasOption("bb") && tag->hasOption("sc") && tag->hasOption("ligand") ) ) {
			throw utility::excn::EXCN_RosettaScriptsOption("'AtrGrid' requires bb, sc, and ligand if any one of these are used");
		}
		bb_= tag->getOption<core::Real>("bb");
		sc_= tag->getOption<core::Real>("sc");
		ligand_= tag->getOption<core::Real>("ligand");
	}

	if ( tag->hasOption("inner_radius") || tag->hasOption("outer_radius") ) {
		// the user MUST provide both if he/she is providing either of these options
		if ( !(tag->hasOption("inner_radius") && tag->hasOption("outer_radius")) ) {
			throw utility::excn::EXCN_RosettaScriptsOption("'AtrGrid' requires outer_radius and inner_radius if either of these options are used");
		}
		inner_radius_= tag->getOption<core::Real>("inner_radius");
		outer_radius_= tag->getOption<core::Real>("outer_radius");
	}
}

utility::json_spirit::Value AtrGrid::serialize()
{
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;

	// [inner_radius, outer_radius]
	Pair radius("radius",utility::tools::make_vector(Value(inner_radius_),Value(outer_radius_)));
	Pair bb("bb",Value(bb_));
	Pair sc("sc",Value(sc_));
	Pair ligand("ligand",Value(ligand_));
	Pair base_data("base_data",SingleGrid::serialize());

	return Value(utility::tools::make_vector(radius,bb,sc,ligand,base_data));

}

void AtrGrid::deserialize(utility::json_spirit::mObject data)
{
	utility::json_spirit::mArray radius(data["radius"].get_array());
	inner_radius_ = radius[0].get_real();
	outer_radius_ = radius[1].get_real();
	bb_ = data["bb"].get_real();
	sc_ = data["sc"].get_real();
	ligand_ = data["ligand"].get_real();
	SingleGrid::deserialize(data["base_data"].get_obj());
}

void AtrGrid::refresh(core::pose::Pose const & pose, core::Vector const & center)
{
	utility::vector1<core::Size> ligand_chain_ids_to_exclude;
	this->refresh(pose, center, ligand_chain_ids_to_exclude);
}

void AtrGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, core::Size const & ligand_chain_id_to_exclude)
{
	utility::vector1<core::Size> ligand_chain_ids_to_exclude;
	ligand_chain_ids_to_exclude.push_back(ligand_chain_id_to_exclude);
	this->refresh(pose, center, ligand_chain_ids_to_exclude);
}

void AtrGrid::set_protein_rings( core::conformation::Residue const & rsd)
{
	for ( core::Size a=1, a_end = rsd.last_backbone_atom(); a <= a_end; ++a ) {
		set_ring(rsd.xyz(a), inner_radius_, outer_radius_, bb_);
	}
	for ( core::Size a = rsd.first_sidechain_atom(), a_end = rsd.nheavyatoms(); a <= a_end; ++a ) {
		set_ring(rsd.xyz(a), inner_radius_, outer_radius_, sc_);
	}

}

void AtrGrid::set_ligand_rings(
	core::conformation::Residue const & rsd,
	utility::vector1<core::Size> ligand_chain_ids_to_exclude

){
	if ( find(
			ligand_chain_ids_to_exclude.begin(),
			ligand_chain_ids_to_exclude.end(),
			rsd.chain()
			) ==  ligand_chain_ids_to_exclude.end()
			) {
		return;
	}
	for ( core::Size a = 1, a_end = rsd.nheavyatoms(); a <= a_end; ++a ) {
		set_ring(rsd.xyz(a), inner_radius_, outer_radius_, ligand_);
	}
}

void AtrGrid::refresh(
	core::pose::Pose const & pose,
	core::Vector const &,
	utility::vector1<core::Size> ligand_chain_ids_to_exclude
){
	// Set neutral core around each sidechain heavy atom, as MOST of these stay put.
	for ( Size r = 1, r_end = pose.total_residue(); r <= r_end; ++r ) {
		core::conformation::Residue const & rsd = pose.residue(r);
		if ( rsd.is_protein() ) set_protein_rings(rsd);
		else {
			set_ligand_rings(rsd, ligand_chain_ids_to_exclude);

		}
		// else don't add this ligand to the grid
	}

}


}
}
}
