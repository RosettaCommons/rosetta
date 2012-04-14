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

#include <protocols/qsar/scoring_grid/HbdGrid.hh>
#include <protocols/qsar/scoring_grid/HbdGridCreator.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
#include <utility/tag/Tag.hh>

#include <core/pose/Pose.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/io/mpistream.hh>


namespace protocols {
namespace qsar {
namespace scoring_grid {


std::string HbdGridCreator::keyname() const
{
	return HbdGridCreator::grid_name();
}

GridBaseOP HbdGridCreator::create_grid(utility::tag::TagPtr const tag) const
{
	GridBaseOP hbd_grid= new HbdGrid();

	hbd_grid->parse_my_tag(tag);

	return hbd_grid;
}

std::string HbdGridCreator::grid_name()
{
	return "HbdGrid";
}

HbdGrid::HbdGrid(): SingleGrid ("HbdGrid",1.0), radius_(2.4),width_(1.0), magnitude_(-1.0)
{

}

HbdGrid::HbdGrid(core::Real weight) : SingleGrid ("HbdGrid",weight), radius_(2.4),width_(1.0), magnitude_(-1.0)
{

}

void
HbdGrid::parse_my_tag(utility::tag::TagPtr const tag){
	if (!tag->hasOption("weight")){
		utility_exit_with_message("Could not make HbdGrid: you must specify a weight when making a new grid");
	}
	set_weight( tag->getOption<core::Real>("weight") );
}

void HbdGrid::refresh(core::pose::Pose const & pose, core::Vector const & )
{
	for(core::Size residue_index = 1; residue_index <=pose.total_residue(); ++residue_index)
	{
		core::conformation::Residue const & residue = pose.residue(residue_index);

		if(!residue.is_protein())
			continue;
		if(residue.has("O"))
			this->diffuse_ring(residue.xyz("O"),radius_,width_,magnitude_);
		//if(residue.has("N"))
		//	this->diffuse_ring(residue.xyz("N"),radius_,width_,magnitude_);


		for(core::Size atom_index=1; atom_index <= residue.natoms(); ++atom_index)
		{
			core::chemical::AtomType atom_type(residue.atom_type(atom_index));

			if(atom_type.is_donor())
			{
				this->diffuse_ring(residue.xyz(atom_index),radius_,width_,magnitude_);
			}
		}
	}
}

void HbdGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, core::Size const & )
{
	refresh(pose,center);
}

void HbdGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, utility::vector1<core::Size> )
{
	refresh(pose,center);
}

}
}
}
