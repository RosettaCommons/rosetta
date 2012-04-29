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

#include <utility/tag/Tag.hh>
#include <utility/tools/make_vector.hh>
#include <utility/json_spirit/json_spirit_writer.h>
#include <utility/json_spirit/json_spirit_reader.h>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/io/mpistream.hh>


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

HbaGrid::HbaGrid() : SingleGrid("HbaGrid",1.0), radius_(2.4),width_(1.0),magnitude_(-1.0)
{

}

HbaGrid::HbaGrid(core::Real weight) : SingleGrid ("HbaGrid",weight), radius_(2.4),width_(1.0), magnitude_(-1.0)
{

}


utility::json_spirit::Value HbaGrid::serialize()
{
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;

	Pair radius_record("radius",Value(radius_));
	Pair width_record("width",Value(width_));
	Pair magnitude_record("mag",Value(magnitude_));
	Pair base_data("base_data",SingleGrid::serialize());

	return Value(utility::tools::make_vector(radius_record,width_record,magnitude_record,base_data));

}

void HbaGrid::deserialize(utility::json_spirit::mObject data)
{
	radius_ = data["radius"].get_real();
	width_ = data["width"].get_real();
	magnitude_ = data["mag"].get_real();

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
	for(core::Size residue_index = 1; residue_index <=pose.total_residue(); ++residue_index)
	{
		core::conformation::Residue const & residue = pose.residue(residue_index);

		if(!residue.is_protein())
			continue;
		//if(residue.has("O"))
		//	this->diffuse_ring(residue.xyz("O"),radius_,width_,magnitude_);
		if(residue.has("N"))
			this->diffuse_ring(residue.xyz("N"),radius_,width_,magnitude_);


		for(core::Size atom_index=1; atom_index <= residue.natoms(); ++atom_index)
		{
			core::chemical::AtomType atom_type(residue.atom_type(atom_index));

			if(atom_type.is_acceptor())
			{
				this->diffuse_ring(residue.xyz(atom_index),radius_,width_,magnitude_);
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

}
}
}
