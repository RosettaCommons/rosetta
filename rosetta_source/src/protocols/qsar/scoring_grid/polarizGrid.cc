// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/qsar/scoring_grid/polarizGrid.cc
/// @author Sam DeLuca

#include <protocols/qsar/scoring_grid/polarizGrid.hh>
#include <protocols/qsar/scoring_grid/polarizGridCreator.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <core/conformation/Residue.hh>
#include <basic/database/open.hh>
#include <basic/Tracer.hh>

#include <core/pose/Pose.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace qsar {
namespace scoring_grid {

static basic::Tracer polarizGridTracer("protocols.ligand_docking.scoring_grid.polarizGrid");

std::string polarizGridCreator::keyname() const
{
	return polarizGridCreator::grid_name();
}

GridBaseOP polarizGridCreator::create_grid(utility::tag::TagPtr const tag) const
{
	GridBaseOP polariz_grid= new polarizGrid();

	polariz_grid->parse_my_tag(tag);

	return polariz_grid;
}

std::string polarizGridCreator::grid_name()
{
	return "polarizGrid";
}

polarizGrid::polarizGrid() : GridBase("polarizGrid",1.0), radius_(3.0)
{
	//TODO: extract this to a method
	utility::io::izstream infile;

	infile.open(basic::database::full_name("scoring/qsar/TPSA.txt",true),std::ifstream::in);
	std::string line;
	do
	{
		getline(infile,line);
		if(line.size() > 0 && line[0] != '#')
		{
			std::vector<std::string> split_line(utility::string_split(line, ' '));

			if(split_line.size() != 2)
			{
				utility_exit_with_message("TPSA data table must have exactly 2 fields per line, check your input");
			}
			core::Real value(utility::from_string(split_line[1],core::Real(0.0)));
			std::pair<std::string, core::Real> polarization_entry(split_line[0],value);
			TPSA_map_.insert(polarization_entry);
		}

	}while(!infile.eof());
	infile.close();
}

polarizGrid::polarizGrid(core::Real weight) : GridBase("polarizGrid",weight), radius_(3.0)
{
	utility::io::izstream infile;

	infile.open(basic::database::full_name("scoring/qsar/TPSA.txt",true),std::ifstream::in);
	std::string line;
	do
	{
		getline(infile,line);
		if(line.size() > 0 && line[0] != '#')
		{
			std::vector<std::string> split_line(utility::string_split(line, ' '));

			if(split_line.size() != 2)
			{
				utility_exit_with_message("TPSA data table must have exactly 2 fields per line, check your input");
			}
			core::Real value(utility::from_string(split_line[1],core::Real(0.0)));
			std::pair<std::string, core::Real> polarization_entry(split_line[0],value);
			TPSA_map_.insert(polarization_entry);
		}

	}while(!infile.eof());
	infile.close();
}

void
polarizGrid::parse_my_tag(utility::tag::TagPtr const tag){
	if (!tag->hasOption("weight")){
		utility_exit_with_message("Could not make PolarizGrid: you must specify a weight when making a new grid");
	}
	set_weight( tag->getOption<core::Real>("weight") );
}

void polarizGrid::refresh(core::pose::Pose const & pose, core::Vector const & )
{
	for(core::Size residue_index=1; residue_index <= pose.total_residue(); ++residue_index)
	{
		core::conformation::Residue const & residue(pose.residue(residue_index));
		core::Vector neighbor_atom_coords(residue.nbr_atom_xyz());

		core::Real polarizability(get_polarizability(residue));
		if(polarizability == -1)
		{
			polarizGridTracer << "WARNING: residue " << residue.name3() << " is not in polarization table."<< std::endl;;

			//utility_exit_with_message("unknown residue");
		}
		else
		{

			this->set_sphere(neighbor_atom_coords,radius_,polarizability);
		}
	}
}


void polarizGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, core::Size const & )
{
	refresh(pose,center);
}

void polarizGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, utility::vector1<core::Size> )
{
	refresh(pose,center);
}

core::Real polarizGrid::get_polarizability(core::conformation::Residue residue)
{
	std::string name(residue.name3());
	utility::trim(name," ");
	std::map<std::string, core::Real>::iterator table_iterator(TPSA_map_.find(name));
	if(table_iterator == TPSA_map_.end())
	{
		return -1;
	}
	else
	{
		return table_iterator->second;
	}
}

}
}
}
