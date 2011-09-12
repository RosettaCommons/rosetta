// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   /git/src/protocols/qsar/scoring_grid/HbaGrid.cc
/// @author Sam DeLuca

#include <protocols/qsar/scoring_grid/HbaGrid.hh>
#include <protocols/qsar/scoring_grid/HbaGridCreator.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
#include <utility/tag/Tag.hh>

namespace protocols {
namespace qsar {
namespace scoring_grid {


std::string HbaGridCreator::keyname() const
{
	return HbaGridCreator::grid_name();
}

GridBaseOP HbaGridCreator::create_grid(utility::tag::TagPtr const tag) const
{
	if (tag->hasOption("weight")){
		return new HbaGrid( tag->getOption<core::Real>("weight") );
	}else{
		return new HbaGrid();
	}
	// This is impossible
	return NULL;
}

std::string HbaGridCreator::grid_name()
{
	return "HbaGrid";
}

HbaGrid::HbaGrid() : GridBase("HbaGrid",0.0), radius_(2.4),width_(0.2),magnitude_(1.0)
{

}

HbaGrid::HbaGrid(core::Real weight) : GridBase ("HbaGrid",weight), radius_(2.4),width_(0.2), magnitude_(1.0)
{

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
