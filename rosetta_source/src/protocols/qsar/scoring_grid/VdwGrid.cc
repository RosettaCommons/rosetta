// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   /git/src/protocols/qsar/scoring_grid/VdwGrid.cc
/// @author Sam DeLuca

#include <protocols/qsar/scoring_grid/VdwGrid.hh>
#include <protocols/qsar/scoring_grid/VdwGridCreator.hh>
#include <protocols/qsar/qsarMap.hh>
#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/AtomGraph.hh>
#include <core/conformation/find_neighbors.hh>
#include <basic/database/open.hh>
#include <numeric/interpolation/util.hh>
#include <utility/tag/Tag.hh>

#include <core/pose/Pose.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace qsar {
namespace scoring_grid {

std::string VdwGridCreator::keyname() const
{
	return VdwGridCreator::grid_name();
}

GridBaseOP VdwGridCreator::create_grid(utility::tag::TagPtr const tag) const
{
	GridBaseOP vdw_grid= new VdwGrid();

	vdw_grid->parse_my_tag(tag);

	return vdw_grid;
}

std::string VdwGridCreator::grid_name()
{
	return "VdwGrid";
}

VdwGrid::VdwGrid() : GridBase("VdwGrid",1.0), cutoff_(10.0)
{
	std::string lj_file(basic::database::full_name("qsar/lj_table.txt"));
	lj_spline_ = numeric::interpolation::spline_from_file(lj_file,0.01);
}

VdwGrid::VdwGrid(core::Real weight) : GridBase ("VdwGrid",weight), cutoff_(10.0)
{
	std::string lj_file(basic::database::full_name("qsar/lj_table.txt"));
	lj_spline_ = numeric::interpolation::spline_from_file(lj_file,0.01);
}

void
VdwGrid::parse_my_tag(utility::tag::TagPtr const tag){
	if (!tag->hasOption("weight")){
		utility_exit_with_message("Could not make VdwGrid: you must specify a weight when making a new grid");
	}
	set_weight( tag->getOption<core::Real>("weight") );
}


void VdwGrid::refresh(core::pose::Pose const & pose, core::Vector const &  )
{
	// loop through all the atoms in the pose
	// get the VDW radius of the atom
	// for each square within cutoff of the atom, update the score
	// continue

	this->fill_with_value(cutoff_);


	for(core::Size residue_index = 1; residue_index <= pose.n_residue(); ++residue_index)
	{
		core::conformation::Residue residue = pose.residue(residue_index);
		for(core::Size atom_index = 1; atom_index <= residue.natoms();++atom_index)
		{
			core::id::AtomID atom_id(atom_index,residue_index);
			core::Vector xyz(pose.xyz(atom_id));
			this->set_distance_sphere_for_atom(xyz,cutoff_);
		}
	}
}

void VdwGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, core::Size const & )
{
	refresh(pose,center);
}

void VdwGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, utility::vector1<core::Size> )
{
	refresh(pose,center);
}


core::Real VdwGrid::score(core::conformation::Residue const & residue, core::Real const max_score, qsarMapOP qsar_map)
{
	numeric::interpolation::spline::InterpolatorOP interpolator(lj_spline_.get_interpolator());
	core::Real score = 0.0;

	for(core::Size atom_index = 1; atom_index <= residue.natoms() && score < max_score; ++atom_index )
	{
		//qsarPointOP qsar_info(qsar_map->get_point(atom_index,"VdwGrid"));
		//if(qsar_info != 0)
		//{
		core::Vector const & atom_coord(residue.xyz(atom_index));
		core::Real const & radius(residue.atom_type(atom_index).lj_radius());
		if(this->get_grid().is_in_grid(atom_coord.x(),atom_coord.y(),atom_coord.z()))
		{
			//core::grid::CartGrid<core::Real>::GridPt grid_point(atom_coord.x(),atom_coord.y(),atom_coord.z());
			core::Real max_radius = this->get_point(atom_coord.x(),atom_coord.y(),atom_coord.z());
			core::Real spline_score = 0.0;
			core::Real spline_score_deriv = 0.0;

			if(max_radius <= lj_spline_.get_lbx())
			{
				spline_score = lj_spline_.get_lby();
			}else if(max_radius >= lj_spline_.get_ubx())
			{
				spline_score = 0.0;
			}else
			{
				interpolator->interpolate(max_radius-radius,spline_score,spline_score_deriv);
			}

			score += spline_score;
		}
		//}
	}
	return score;
}



}
}
}
