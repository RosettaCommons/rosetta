// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/qsar/scoring_grid/RepGrid.cc
/// @author Sam DeLuca

#include <protocols/qsar/scoring_grid/RepGrid.hh>
#include <protocols/qsar/scoring_grid/RepGridCreator.hh>
#include <protocols/qsar/qsarMap.fwd.hh>
#include <core/conformation/Residue.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>

#include <core/pose/Pose.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector0.hh>


namespace protocols {
namespace qsar {
namespace scoring_grid {

static basic::Tracer RepGridTracer("protocols.ligand_docking.scoring_grid.RepGrid");

std::string RepGridCreator::keyname() const
{
	return RepGridCreator::grid_name();
}

GridBaseOP RepGridCreator::create_grid(utility::tag::TagCOP const tag) const
{
	GridBaseOP rep_grid= new RepGrid();

	rep_grid->parse_my_tag(tag);

	return rep_grid;
}

GridBaseOP RepGridCreator::create_grid() const
{
	return new RepGrid();
}

std::string RepGridCreator::grid_name()
{
	return "RepGrid";
}

RepGrid::RepGrid() : SingleGrid("RepGrid"), radius_(2.25), bb_(1), sc_(0), ligand_(1)
{

}

utility::json_spirit::Value RepGrid::serialize()
{
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;

	// [inner_radius, outer_radius]
	Pair radius("radius",radius_);
	Pair bb("bb",Value(bb_));
	Pair sc("sc",Value(sc_));
	Pair ligand("ligand",Value(ligand_));
	Pair base_data("base_data",SingleGrid::serialize());

	return Value(utility::tools::make_vector(radius,bb,sc,ligand,base_data));

}

void RepGrid::deserialize(utility::json_spirit::mObject data)
{
	radius_ = data["radius"].get_real();
	bb_ = data["bb"].get_real();
	sc_ = data["sc"].get_real();
	ligand_ = data["ligand"].get_real();
	SingleGrid::deserialize(data["base_data"].get_obj());
}

void
RepGrid::parse_my_tag(utility::tag::TagCOP const tag){

	if (tag->hasOption("bb") || tag->hasOption("sc") || tag->hasOption("ligand") ){
		// the user MUST provide all 3 if he/she is providing any of these 3 options
		if (!(tag->hasOption("bb") && tag->hasOption("sc") && tag->hasOption("ligand") ) ){
			throw utility::excn::EXCN_RosettaScriptsOption("'RepGrid' requires bb, sc, and ligand if any one of these are used");
		}else{
			bb_= tag->getOption<core::Real>("bb");
			sc_= tag->getOption<core::Real>("sc");
			ligand_= tag->getOption<core::Real>("ligand");
		}
	}
}

void RepGrid::refresh( core::pose::Pose const & pose,  core::Vector const & )
{
	for(core::Size residue_index = 1; residue_index <= pose.total_residue(); ++residue_index)
	{
		//RepGridTracer <<"refreshing residue " <<residue_index << " of " << pose.total_residue() <<std::endl;
		core::conformation::Residue const & residue = pose.residue(residue_index);
		if(!residue.is_protein())
			continue;

		if(residue.has("CB"))
			this->set_sphere(residue.xyz("CB"),this->radius_,1.0);
		if(residue.has("N"))
			this->set_sphere(residue.xyz("N"),this->radius_,1.0);
		if(residue.has("CA"))
			this->set_sphere(residue.xyz("CA"),this->radius_,1.0);
		if(residue.has("C"))
			this->set_sphere(residue.xyz("C"),this->radius_,1.0);
		if(residue.has("O"))
			this->set_sphere(residue.xyz("O"),this->radius_,1.0);
	}
}

void RepGrid::refresh( core::pose::Pose const & pose,  core::Vector const & center, const core::Size & )
{
	//for the repulsive force, the case of no ligands and all ligands are identical
	refresh(pose, center);
}

void RepGrid::refresh(core::pose::Pose const & pose, core::Vector const & center,utility::vector1<core::Size> )
{
	//for the repulsive force, the case of no ligands and all ligands are identical
	refresh(pose,center);
}


}
}
}
