// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   /git/src/protocols/qsar/scoring_grid/AtrGrid.cc
/// @author Sam DeLuca

#include <utility/tag/Tag.hh>
#include <protocols/qsar/scoring_grid/AtrGrid.hh>
#include <protocols/qsar/scoring_grid/AtrGridCreator.hh>
#include <core/conformation/Residue.hh>

namespace protocols {
namespace qsar {
namespace scoring_grid {

std::string AtrGridCreator::keyname() const
{
	return AtrGridCreator::grid_name();
}

GridBaseOP AtrGridCreator::create_grid(utility::tag::TagPtr const tag) const
{
	if (!tag->hasOption("weight")){
		utility_exit_with_message("Could not make AtrGrid: you must specify a weight when making a new grid");

	}else{
		return new AtrGrid( tag->getOption<core::Real>("weight") );
	}
	// This is impossible
	return NULL;
}

std::string AtrGridCreator::grid_name()
{
	return "AtrGrid";
}

AtrGrid::AtrGrid() : GridBase("AtrGrid",1.0),radius_(4.75)
{
	//
}

AtrGrid::AtrGrid(core::Real weight) : GridBase ("AtrGrid",weight), radius_(4.75)
{
 //
}

void AtrGrid::refresh(core::pose::Pose const & pose, core::Vector const & )
{

	for(core::Size residue_id=1 ; residue_id <= pose.total_residue(); ++residue_id)
	{
		core::conformation::Residue const & residue = pose.residue(residue_id);
		if(!residue.is_protein())
			continue;
		for(core::Size atom_index = 1; atom_index <= residue.nheavyatoms();++atom_index)
		{
			this->set_sphere(residue.xyz(atom_index),radius_, -1.0);
		}
	}

	for(core::Size residue_id=1; residue_id <= pose.total_residue(); ++residue_id)
	{
		core::conformation::Residue const & residue = pose.residue(residue_id);
		if(!residue.is_protein())
			continue;
		for(core::Size atom_index = 1; atom_index <= residue.nheavyatoms();++atom_index)
		{
			this->set_sphere(residue.xyz(atom_index),2.25, -1.0);
		}
	}
}

void AtrGrid::refresh(core::pose::Pose const & pose, core::Vector const & , core::Size const & ligand_chain_id_to_exclude)
{
	for(core::Size residue_id=1; residue_id <= pose.total_residue(); ++residue_id)
	{
		core::conformation::Residue const & residue = pose.residue(residue_id);
		if(residue.chain() == ligand_chain_id_to_exclude)
			continue;
		for(core::Size atom_index = 1; atom_index <= residue.nheavyatoms();++atom_index)
		{
			this->set_sphere(residue.xyz(atom_index),radius_, -1.0);
		}
	}

	for(core::Size residue_id=1; residue_id <= pose.total_residue(); ++residue_id)
	{
		core::conformation::Residue const & residue = pose.residue(residue_id);
		if(residue.chain() == ligand_chain_id_to_exclude)
			continue;
		for(core::Size atom_index = 1; atom_index <= residue.nheavyatoms();++atom_index)
		{
			this->set_sphere(residue.xyz(atom_index),2.25, -1.0);
		}
	}
}

void AtrGrid::refresh(core::pose::Pose const & pose, core::Vector const & , utility::vector1<core::Size> ligand_chain_ids_to_exclude)
{
	for(core::Size residue_id=1; residue_id <= pose.total_residue(); ++residue_id)
	{
		core::conformation::Residue const & residue = pose.residue(residue_id);
		if(find(
				ligand_chain_ids_to_exclude.begin(),
				ligand_chain_ids_to_exclude.end(),
				residue.chain()) == ligand_chain_ids_to_exclude.end())
		{
			continue;
		}

		for(core::Size atom_index = 1; atom_index <= residue.nheavyatoms();++atom_index)
		{
			this->set_sphere(residue.xyz(atom_index),radius_, -1.0);
		}
	}

	for(core::Size residue_id=1; residue_id <= pose.total_residue(); ++residue_id)
	{
		core::conformation::Residue const & residue = pose.residue(residue_id);
		if(find(
				ligand_chain_ids_to_exclude.begin(),
				ligand_chain_ids_to_exclude.end(),
				residue.chain()) == ligand_chain_ids_to_exclude.end())
		{
			continue;
		}
		for(core::Size atom_index = 1; atom_index <= residue.nheavyatoms();++atom_index)
		{
			this->set_sphere(residue.xyz(atom_index),2.25, -1.0);
		}
	}

}


}
}
}
