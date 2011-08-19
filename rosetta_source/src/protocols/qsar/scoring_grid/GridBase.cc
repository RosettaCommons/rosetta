// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/scoring_grid/GridBase.cc
/// @author Sam DeLuca

#include <protocols/qsar/scoring_grid/GridBase.hh>
#include <core/conformation/Residue.hh>
#include <protocols/qsar/qsarMap.hh>
#include <core/grid/CartGrid.hh>
#include <protocols/ligand_docking/grid_functions.hh>
#include <basic/Tracer.hh>
//#include <protocols/qsar/qsarTypeManager.hh>

namespace protocols {
namespace qsar {
namespace scoring_grid {

static basic::Tracer GridBaseTracer("protocols.qsar.scoring_grid.GridBase");


GridBase::GridBase(std::string type, core::Real weight) : type_(type), weight_(weight)
{

}

void GridBase::initialize(core::Vector const & center, core::Real width,core::Real resolution )
{
	center_=center;

	core::Size num_pts = static_cast<core::Size>(width/resolution);
	//core::Size num_pts = 160;
	//core::Real const resolution = 0.1;
	core::Real const grid_halfwidth = width / 2.0;
	grid_.setBase(
			center.x() - grid_halfwidth,
			center.y() - grid_halfwidth,
			center.z() - grid_halfwidth
	);
	grid_.setDimensions(num_pts,num_pts,num_pts,resolution,resolution,resolution);
	grid_.setupZones();
	grid_.zero();
}

core::grid::CartGrid<core::Real> const & GridBase::get_grid()
{
	return grid_;
}

void GridBase::set_type(std::string type)
{
	type_ = type;
}

std::string GridBase::get_type()
{
	return type_;
}

void GridBase::set_weight(core::Real weight)
{
	weight_ = weight;
}

core::Real GridBase::get_weight()
{
	return weight_;
}

void GridBase::set_center(core::Vector center)
{
	center_ = center;
}

core::Vector GridBase::get_center()
{
	return center_;
}

core::Real GridBase::get_point(core::Real x, core::Real y, core::Real z)
{
	 return grid_.getValue(x,y,z);
}

numeric::xyzVector<core::Size> GridBase::get_dimensions()
{
	int x_size(0);
	int y_size(0);
	int z_size(0);

	grid_.getNumberOfPoints(x_size,y_size,z_size);
	return numeric::xyzVector<core::Size>(x_size,y_size,z_size);
}

core::Real GridBase::score(core::conformation::Residue const & residue, core::Real const max_score,qsarMapOP qsar_map)
{
	core::Real score = 0.0;
	//GridBaseTracer << "map size is: " << qsar_map->size() <<std::endl;
	for(core::Size atom_index = 1; atom_index <= residue.nheavyatoms() && score < max_score;++atom_index)
	{
		qsarPointOP qsar_info(qsar_map->get_point(atom_index,type_));

		if(qsar_info != 0)
		{
			core::Vector const & atom = residue.xyz(atom_index);
			if(grid_.is_in_grid(atom.x(),atom.y(), atom.z()))
			{

				core::Real grid_score = grid_.getValue(atom.x(),atom.y(),atom.z());

				score = score+ grid_score*qsar_info->get_value();
			}
		}

	}
	return score*weight_;
}

void GridBase::grid_to_kin(utility::io::ozstream & out, core::Real min_val, core::Real max_val, core::Size stride)
{
	int x_size(0);
	int y_size(0);
	int z_size(0);
	grid_.getNumberOfPoints(x_size, y_size, z_size);
	for(int x_index =0; x_index < x_size;x_index +=stride)
	{
		for(int y_index = 0; y_index < y_size; y_index += stride)
		{
			for(int z_index = 0; z_index < z_size; z_index += stride)
			{
				core::grid::CartGrid<core::Real>::GridPt grid_point(x_index,y_index,z_index);
				core::Vector box_counter = grid_.coords(grid_point);
				core::Real value = grid_.getValue(grid_point);
				if(min_val <= value && value <= max_val)
				{
					out << '{' << x_index << ' ' << y_index << ' ' << z_index << "}U "
						<< box_counter.x() << ' ' << box_counter.y() << ' ' << box_counter.z() << '\n';
				}
			}
		}
	}
}

/*
void GridBase::grid_rotamer_trials(core::pose::Pose &  pose, core::Size residue_id, int const min_score)
{
	utility::vector1<core::conformation::ResidueOP> conformers;
	ligand_docking::rotamers_for_trials(pose,residue_id,conformers);
	if(conformers.empty())
		return;
	int best_score = 9999;
	int best_conformer = 0;
	for(core::Size current_conformer = 1; current_conformer <= conformers.size(); ++current_conformer )
	{
		core::conformation::ResidueOP residue = conformers[current_conformer];
		int const new_score = this->score(*residue,best_score,qsa);
		if(new_score < best_score && new_score >= min_score)
		{
			best_score=new_score;
			best_conformer = current_conformer;
			if(best_score == min_score)
				break;
		}
	}

	if(best_conformer > 0)
	{
		TR << "best fit is conformer " << best_conformer << "with score" << best_score << std::endl;
		pose.replace_residue(residue_id,*conformers[best_conformer],false);
	}

}
*/

void GridBase::set_sphere(core::Vector const & coords, core::Real radius, core::Real value)
{
	//TR <<"making sphere of radius " << radius << "and value " << value <<std::endl;
	core::Real radius2 = radius*radius;
	int x_count(0);
	int y_count(0);
	int z_count(0);
	grid_.getNumberOfPoints(x_count,y_count,z_count);
	core::Vector vector_radius (radius);
	core::grid::CartGrid<core::Real>::GridPt grid_min = grid_.gridpt(coords - radius);
	core::grid::CartGrid<core::Real>::GridPt grid_max = grid_.gridpt(coords + radius);
	//TR <<"min: " <<grid_min.x() << " "<<grid_min.y() << " " <<grid_min.z() <<std::endl;
	//TR <<"max: "<<grid_max.x() << " "<<grid_max.y() << " " <<grid_max.z() <<std::endl;
	//TR <<"counts: " << x_count <<" "<< y_count << " " << z_count <<std::endl;
	for(int x_index = std::max(0,grid_min.x()); x_index <= std::min(x_count-1,grid_max.x());++x_index)
	{
		for(int y_index = std::max(0,grid_min.y());y_index <= std::min(y_count-1,grid_max.y());++y_index)
		{
			for(int z_index = std::max(0,grid_min.z()); z_index <= std::min(z_count-1,grid_max.z());++z_index)
			{
				core::grid::CartGrid<core::Real>::GridPt point(x_index, y_index, z_index);
				core::Vector box_center = grid_.coords(point);
				//TR <<x_index << " "<<y_index << " " <<z_index <<std::endl;
				if(box_center.distance_squared(coords) <= radius2)
				{
					grid_.setValue(point, value);
				}
			}
		}
	}
	//TR << "done making sphere "<<std::endl;
}



void GridBase::diffuse_ring(core::Vector const & coords, core::Real radius, core::Real width, core::Real magnitude)
{
	core::Real radius_squared = radius*radius;
	core::Real half_width = width/2;
	core::Real half_width_squared = half_width*half_width;
	int x_count(0);
	int y_count(0);
	int z_count(0);
	grid_.getNumberOfPoints(x_count,y_count,z_count);
	core::Vector vector_radius (radius);
	core::grid::CartGrid<core::Real>::GridPt grid_min = grid_.gridpt(coords - radius);
	core::grid::CartGrid<core::Real>::GridPt grid_max = grid_.gridpt(coords + radius);
	for(int x_index = std::max(0,grid_min.x()); x_index <= std::min(x_count-1,grid_max.x());++x_index)
	{
		for(int y_index = std::max(0,grid_min.y());y_index <= std::min(y_count-1,grid_max.y());++y_index)
		{
			for(int z_index = std::max(0,grid_min.z()); z_index <= std::min(z_count-1,grid_max.z());++z_index)
			{
				core::grid::CartGrid<core::Real>::GridPt point(x_index,y_index,z_index);
				core::Vector origin = grid_.coords(point);
				core::Real distance_squared = origin.distance_squared(coords);
				//core::Real rad_distance_squared = abs(radius2-distance_squared);

				//if((distance_squared < radius2 + half_width) || (distance_squared > radius2-half_width))
				//if(distance_squared <radius2)
				if((distance_squared >radius_squared-half_width_squared && distance_squared < radius_squared)|| (distance_squared <radius_squared+half_width_squared && distance_squared > radius_squared))
				{
					//TR << "inserting value " <<std::endl;
					grid_.setValue(point,magnitude);
				}else
				{

					//TR<<  "d^2 is " << distance_squared <<std::endl;
					core::Real distance_from_ring = 0.0;
					if(distance_squared < radius_squared)
					{
						distance_from_ring = (radius_squared-half_width_squared)-distance_squared;
					}else
					{
						distance_from_ring = distance_squared-(radius_squared+half_width_squared);
					}

					if(distance_from_ring <= magnitude)
					{
						grid_.setValue(point,magnitude-distance_from_ring);
					}
					else
					{
					}

				}

			}
		}
	}
}

void GridBase::set_point(core::Vector const & coords, core::Real value)
{
	grid_.setValue(coords,value);
}

void GridBase::dump_BRIX(std::string const & prefix)
{
	//std::string grid_type_string(qsar::qsarTypeManager::name_from_qsar_type(type_));
	grid_.write_to_BRIX(prefix + type_ + ".omap");
}

}
}
}
