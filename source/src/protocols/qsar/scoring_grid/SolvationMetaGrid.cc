// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/ligand_docking/scoring_grid/SolvationMetaGrid.cc
/// @author Sam DeLuca

#include <protocols/qsar/scoring_grid/SolvationMetaGrid.hh>
#include <protocols/qsar/scoring_grid/SolvationMetaGridCreator.hh>

#include <protocols/qsar/scoring_grid/SingleGrid.hh>
#include <protocols/qsar/scoring_grid/SolvationGrid.fwd.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/UltraLightResidue.hh>

#include <utility/tag/Tag.hh>

#include <map>

namespace protocols {
namespace qsar {
namespace scoring_grid {

std::string SolvationMetaGridCreator::keyname() const
{
	return SolvationMetaGridCreator::grid_name();
}

GridBaseOP SolvationMetaGridCreator::create_grid(utility::tag::TagCOP tag) const
{
	GridBaseOP solvation_meta_grid = new SolvationMetaGrid();
	solvation_meta_grid->parse_my_tag(tag);
	return solvation_meta_grid;
}

GridBaseOP SolvationMetaGridCreator::create_grid() const
{
	return new SolvationMetaGrid();
}

std::string SolvationMetaGridCreator::grid_name()
{
	return "SolvationMetaGrid";
}

SolvationMetaGrid::SolvationMetaGrid() :type_("SolvationMetaGrid")
{

}

SolvationMetaGrid::~SolvationMetaGrid()
{

}

void SolvationMetaGrid::initialize(core::Vector const & center, core::Real width, core::Real resolution)
{
    core::chemical::AtomTypeSetCOP atom_type_set(
        core::chemical::ChemicalManager::get_instance()->atom_type_set("fa_standard"));

    core::ShortSize max_atom_type = atom_type_set->n_atomtypes();
    for(core::ShortSize atom_type = 1; atom_type <= max_atom_type;++atom_type)
    {

    	SolvationGridOP new_grid = new SolvationGrid();
    	new_grid->set_probe_atom_type(atom_type);
    	new_grid->initialize(center, width, resolution);
    	grid_map_[atom_type] = SingleGridOP(new_grid);
    }
}

void SolvationMetaGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, core::Size const & )
{
	refresh(pose,center);
}

void SolvationMetaGrid::refresh(core::pose::Pose const & pose, core::Vector const & center,utility::vector1<core::Size> )
{
	refresh(pose,center);
}

void SolvationMetaGrid::refresh(core::pose::Pose const & pose, core::Vector const & center)
{
	std::map<core::ShortSize,SingleGridOP>::iterator it = grid_map_.begin();
	for(;it != grid_map_.end();++it)
	{
		std::cout << "initializing solvation grid for atom type " <<it->first <<std::endl;
		it->second->refresh(pose,center);
	}
}

void SolvationMetaGrid::parse_my_tag(utility::tag::TagCOP const /*tag*/)
{
}

core::Real SolvationMetaGrid::score(
		core::conformation::UltraLightResidue const & residue,
		core::Real const /*max_score*/,
		qsarMapOP /*qsar_map*/)
{
	core::Real total_score = 0.0;
	for(core::Size atom_index = 1; atom_index <= residue.natoms();++atom_index)
	{
		//core::Vector const & coords = residue.xyz(atom_index);
		core::conformation::Atom current_atom(residue.residue()->atom(atom_index));
		std::map<core::ShortSize,SingleGridOP>::iterator grid_iterator(grid_map_.find(current_atom.type()));
		if(grid_iterator == grid_map_.end())
		{
			utility_exit_with_message("Ligands must be parameterized with the FA_STANDARD atom type set for use in the SolvationMetaGrid");
		}

		SingleGridOP current_grid = grid_iterator->second;
		total_score += current_grid->get_point(residue[atom_index]);
	}

	return total_score;
}

core::Real SolvationMetaGrid::atom_score(
		core::conformation::UltraLightResidue const & residue,
		core::Size atomno,
		qsarMapOP /*qsar_map*/)
{
	core::conformation::Atom current_atom(residue.residue()->atom(atomno));
	std::map<core::ShortSize,SingleGridOP>::iterator grid_iterator(grid_map_.find(current_atom.type()));
	if(grid_iterator == grid_map_.end())
	{
		utility_exit_with_message("Ligands must be parameterized with the FA_STANDARD atom type set for use in the SolvationMetaGrid");
	}

	SingleGridOP current_grid = grid_iterator->second;
	return  current_grid->get_point(residue[atomno]);
}

core::Real SolvationMetaGrid::score(core::conformation::Residue const & residue, core::Real const /*max_score*/, qsarMapOP)
{
	core::Real total_score = 0.0;
	for(core::Size atom_index = 1; atom_index <= residue.natoms();++atom_index)
	{
		//core::Vector const & coords = residue.xyz(atom_index);
		core::conformation::Atom current_atom(residue.atom(atom_index));
		std::map<core::ShortSize,SingleGridOP>::iterator grid_iterator(grid_map_.find(current_atom.type()));
		if(grid_iterator == grid_map_.end())
		{
			utility_exit_with_message("Ligands must be parameterized with the FA_STANDARD atom type set for use in the SolvationMetaGrid");
		}

		SingleGridOP current_grid = grid_iterator->second;
		total_score += current_grid->get_point(current_atom.xyz());
	}

	return total_score;
}

core::Real SolvationMetaGrid::atom_score(core::conformation::Residue const & residue, core::Size atomno, qsarMapOP /*qsar_map*/)
{
	core::conformation::Atom current_atom(residue.atom(atomno));
	std::map<core::ShortSize,SingleGridOP>::iterator grid_iterator(grid_map_.find(current_atom.type()));
	if(grid_iterator == grid_map_.end())
	{
		utility_exit_with_message("Ligands must be parameterized with the FA_STANDARD atom type set for use in the SolvationMetaGrid");
	}

	SingleGridOP current_grid = grid_iterator->second;
	return  current_grid->get_point(current_atom.xyz());
}

std::string SolvationMetaGrid::get_type()
{
	return "SolvationMetaGrid";
}

void SolvationMetaGrid::set_chain(char chain)
{
	std::map<core::ShortSize,SingleGridOP>::iterator it = grid_map_.begin();
	for(;it != grid_map_.end();++it)
	{
		it->second->set_chain(chain);
	}
}

void SolvationMetaGrid::dump_BRIX(std::string const & /*prefix*/)
{
	utility_exit_with_message("SolvationMetaGrid is currently unable to output a BRIX grid, sorry :(");
}

utility::json_spirit::Value SolvationMetaGrid::serialize()
{
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;

	Pair type_record("type",Value(type_));
	std::vector<Value> grid_map_vector;

	for(std::map<core::ShortSize,SingleGridOP>::iterator it = grid_map_.begin();it != grid_map_.end();++it)
	{
		std::vector<Value> current_map_data(utility::tools::make_vector(Value(core::ShortSize(it->first)),Value(it->second->serialize())));
		grid_map_vector.push_back(Value(current_map_data));
	}

	Pair grid_map_data("grids",grid_map_vector);

	#ifdef PYROSETTA
		Value _;  return _;
	#endif

	return Value(utility::tools::make_vector(type_record,grid_map_data));

}

void SolvationMetaGrid::deserialize(utility::json_spirit::mObject data)
{
	type_ = data["type"].get_str();
	utility::json_spirit::mArray grid_map_data(data["grids"].get_array());
	for(utility::json_spirit::mArray::iterator it = grid_map_data.begin();it != grid_map_data.end();++it)
	{
		utility::json_spirit::mArray grid_pair(it->get_array());
		core::ShortSize atom_type = grid_pair[0].get_int();
		SingleGridOP grid = new SolvationGrid();
		grid->deserialize(grid_pair[1].get_obj());
		grid_map_[atom_type] = grid;
	}
}

bool SolvationMetaGrid::is_in_grid(core::conformation::UltraLightResidue const & residue)
{
	for(std::map<core::ShortSize,SingleGridOP>::iterator it = grid_map_.begin();it != grid_map_.end();++it)
	{
		if(!it->second->is_in_grid(residue))
		{
			return false;
		}
	}
	return true;
}

bool SolvationMetaGrid::is_in_grid(core::conformation::Residue const & residue)
{
	for(std::map<core::ShortSize,SingleGridOP>::iterator it = grid_map_.begin();it != grid_map_.end();++it)
	{
		if(!it->second->is_in_grid(residue))
		{
			return false;
		}
	}
	return true;
}

}
}
}
