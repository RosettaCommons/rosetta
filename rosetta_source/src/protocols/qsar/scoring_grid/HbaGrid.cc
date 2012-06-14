// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/qsar/scoring_grid/HbaGrid.cc
/// @author Sam DeLuca

#include <protocols/qsar/scoring_grid/HbaGrid.hh>
#include <protocols/qsar/scoring_grid/HbaGridCreator.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
#include <core/id/AtomID.hh>

#include <utility/tag/Tag.hh>
#include <utility/tools/make_vector.hh>
#include <utility/json_spirit/json_spirit_value.h>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/io/mpistream.hh>

#include <basic/database/open.hh>

#include <numeric/interpolation/util.hh>

namespace protocols {
namespace qsar {
namespace scoring_grid {


std::string HbaGridCreator::keyname() const
{
	return HbaGridCreator::grid_name();
}

GridBaseOP HbaGridCreator::create_grid(utility::tag::TagPtr const tag) const
{
	GridBaseOP hba_grid= new HbaGrid();

	hba_grid->parse_my_tag(tag);

	return hba_grid;
}

GridBaseOP HbaGridCreator::create_grid() const
{
	return new HbaGrid();
}

std::string HbaGridCreator::grid_name()
{
	return "HbaGrid";
}

HbaGrid::HbaGrid() : SingleGrid("HbaGrid",1.0)
{
	std::string lj_file(basic::database::full_name("scoring/qsar/hb_table.txt"));
	lj_spline_ = numeric::interpolation::spline_from_file(lj_file,0.05).get_interpolator();
}

HbaGrid::HbaGrid(core::Real weight) : SingleGrid ("HbaGrid",weight)
{
	std::string lj_file(basic::database::full_name("scoring/qsar/hb_table.txt"));
	lj_spline_ = numeric::interpolation::spline_from_file(lj_file,0.05).get_interpolator();
}

HbaGrid::~HbaGrid()
{

}

utility::json_spirit::Value HbaGrid::serialize()
{
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;

	Pair spline_data("spline",lj_spline_->serialize());
	Pair base_data("base_data",SingleGrid::serialize());

	return Value(utility::tools::make_vector(spline_data,base_data));

}

void HbaGrid::deserialize(utility::json_spirit::mObject data)
{
	lj_spline_->deserialize(data["spline"].get_obj());
	SingleGrid::deserialize(data["base_data"].get_obj());
}

void
HbaGrid::parse_my_tag(utility::tag::TagPtr const tag){
	if (!tag->hasOption("weight")){
		utility_exit_with_message("Could not make HbaGrid: you must specify a weight when making a new grid");
	}
	set_weight( tag->getOption<core::Real>("weight") );
}

void HbaGrid::refresh(core::pose::Pose const & pose, core::Vector const & )
{

	this->fill_with_value(0.0);

	for(core::Size residue_index = 1; residue_index <= pose.total_residue();++residue_index)
	{
		core::conformation::Residue const residue = pose.residue(residue_index);
		if(!residue.is_protein())
		{
			continue;
		}
		for(core::Size atom_index=1; atom_index <= residue.natoms();++atom_index)
		{
			core::chemical::AtomType atom_type(residue.atom_type(atom_index));
			if(atom_type.is_acceptor())
			{
				core::id::AtomID atom_id(atom_index,residue_index);
				core::Vector xyz(pose.xyz(atom_id));
				this->set_score_sphere_for_atom(lj_spline_,xyz,5.0);
			}
		}
	}
}

void HbaGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, core::Size const & )
{
	refresh(pose,center);
}

void HbaGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, utility::vector1<core::Size> )
{
	refresh(pose,center);
}

core::Real HbaGrid::score(core::conformation::Residue const & residue, core::Real const max_score, qsarMapOP qsar_map)
{
	core::Real score = 0.0;
	//GridBaseTracer << "map size is: " << qsar_map->size() <<std::endl;
	for(core::Size atom_index = 1; atom_index <= residue.nheavyatoms() && score < max_score;++atom_index)
	{
		core::Vector const & atom_coord(residue.xyz(atom_index));
		if(this->get_grid().is_in_grid(atom_coord.x(),atom_coord.y(),atom_coord.z()))
		{
			core::chemical::AtomType atom_type(residue.atom_type(atom_index));
			if(atom_type.is_donor())
			{
				core::Real grid_value = this->get_point(atom_coord.x(),atom_coord.y(),atom_coord.z());
				score += grid_value;
			}

		}
	}

	return score*this->get_weight();
}

}
}
}
