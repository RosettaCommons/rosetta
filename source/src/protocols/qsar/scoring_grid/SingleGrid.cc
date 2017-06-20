// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/scoring_grid/GridBase.cc
/// @author Sam DeLuca

#include <protocols/qsar/scoring_grid/SingleGrid.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/UltraLightResidue.hh>
#include <core/grid/CartGrid.hh>
#include <core/pose/util.hh>
#include <core/conformation/Conformation.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/tools/make_vector.hh>
#include <utility/json_spirit/json_spirit_value.h>

#include <numeric/xyz.json.hh>
#include <numeric/xyzVector.io.hh>

#include <algorithm>

namespace protocols {
namespace qsar {
namespace scoring_grid {

static THREAD_LOCAL basic::Tracer TR( "protocols.qsar.scoring_grid.SingleGrid" );


SingleGrid::SingleGrid(std::string type) : type_(type),chain_('A')
{

}

SingleGrid::~SingleGrid()
{

}

void SingleGrid::initialize(core::Vector const & center, core::Real width,core::Real resolution )
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

utility::json_spirit::Value SingleGrid::serialize() const
{
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;

	Pair type_record("type",Value(type_));
	Pair center_record("center",numeric::serialize(center_));
	Pair chain_record("chain",Value(chain_));
	Pair grid_record("grid_data",grid_.serialize());

	return Value(utility::tools::make_vector(type_record,center_record,chain_record,grid_record));

}

void SingleGrid::deserialize(utility::json_spirit::mObject data)
{


	type_ = data["type"].get_str();
	center_ = numeric::deserialize<core::Real>(data["center"].get_array());

	grid_.deserialize(data["grid_data"].get_obj());

}

void SingleGrid::set_chain(char chain)
{
	chain_ = chain;
}

char SingleGrid::get_chain()
{
	return chain_;
}

core::grid::CartGrid<core::Real> const & SingleGrid::get_grid() const
{
	return grid_;
}

void SingleGrid::set_type(std::string type)
{
	type_ = type;
}

std::string SingleGrid::get_type() const
{
	return type_;
}

void SingleGrid::set_center(core::Vector center)
{
	center_ = center;
}

core::Vector SingleGrid::get_center() const
{
	return center_;
}

core::Real SingleGrid::get_max_value() const
{
	return grid_.getMaxValue();
}

core::Real SingleGrid::get_min_value() const
{
	return grid_.getMinValue();
}

core::Real SingleGrid::get_point(core::Real x, core::Real y, core::Real z) const
{
	if ( grid_.is_in_grid(x,y,z) ) {
		return grid_.getValue(x,y,z);
	} else {
		return 0.0;
	}
}

core::Real SingleGrid::get_point(core::Vector coords) const
{
	return grid_.getValue(coords);
}


numeric::xyzVector<core::Size> SingleGrid::get_dimensions()
{
	core::Size x_size(0);
	core::Size y_size(0);
	core::Size z_size(0);

	grid_.getNumberOfPoints(x_size,y_size,z_size);
	return numeric::xyzVector<core::Size>(x_size,y_size,z_size);
}

core::Vector SingleGrid::get_pdb_coords(int x, int y, int z)
{
	core::grid::CartGrid<core::Real>::GridPt gridpt(x,y,z);
	return get_pdb_coords(gridpt);
}

core::Vector SingleGrid::get_pdb_coords(core::grid::CartGrid<core::Real>::GridPt gridpt)
{
	return grid_.coords(gridpt);
}

core::Real SingleGrid::score(core::conformation::UltraLightResidue const & residue, core::Real const max_score, qsarMapCOP) const
{
	core::Real score = 0.0;
	//TR << "map size is: " << qsar_map->size() <<std::endl;
	for ( core::Size atom_index = 1; atom_index <= residue.natoms() && score < max_score; ++atom_index ) {
		//TODO qsar map is broken, comment it out until it works right
		//qsarPointOP qsar_info(qsar_map->get_point(atom_index,type_));

		//if(qsar_info != 0)
		//{
		core::Vector const & atom = residue[atom_index];
		if ( grid_.is_in_grid(atom.x(),atom.y(), atom.z()) ) {

			core::Real grid_score = grid_.getValue(atom.x(),atom.y(),atom.z());

			score = score+ grid_score; //*qsar_info->get_value();
		}
		//}

	}
	return score;
}

core::Real SingleGrid::atom_score(core::conformation::UltraLightResidue const & residue, core::Size atomno, qsarMapCOP) const
{
	core::Vector const & atom = residue[atomno];
	if ( grid_.is_in_grid(atom.x(),atom.y(), atom.z()) ) {

		core::Real grid_score = grid_.getValue(atom.x(),atom.y(),atom.z());
		return grid_score;
	} else {
		return 0;
	}
}

core::Real SingleGrid::score(core::conformation::Residue const & residue, core::Real const max_score,qsarMapCOP /*qsar_map*/) const
{
	core::Real score = 0.0;
	//TR << "map size is: " << qsar_map->size() <<std::endl;
	for ( core::Size atom_index = 1; atom_index <= residue.natoms() && score < max_score; ++atom_index ) {
		//TODO qsar map is broken, comment it out until it works right
		//qsarPointOP qsar_info(qsar_map->get_point(atom_index,type_));

		//if(qsar_info != 0)
		//{
		core::Vector const & atom = residue.xyz(atom_index);
		if ( grid_.is_in_grid(atom.x(),atom.y(), atom.z()) ) {

			core::Real grid_score = grid_.getValue(atom.x(),atom.y(),atom.z());

			score = score+ grid_score; //*qsar_info->get_value();
		}
		//}

	}
	return score;
}

core::Real SingleGrid::atom_score(core::conformation::Residue const & residue, core::Size atomno, qsarMapCOP /*qsar_map*/) const
{
	core::Vector const & atom = residue.xyz(atomno);
	if ( grid_.is_in_grid(atom.x(),atom.y(), atom.z()) ) {

		core::Real grid_score = grid_.getValue(atom.x(),atom.y(),atom.z());
		return grid_score;
	} else {
		return 0;
	}

}


std::list<std::pair<core::Vector, core::Real> >
SingleGrid::get_point_value_list_within_range(
	core::Real lower_bound,
	core::Real upper_bound,
	core::Size stride) const
{
	std::list<std::pair<core::Vector, core::Real> > point_list;
	core::Size x_size(0);
	core::Size y_size(0);
	core::Size z_size(0);
	grid_.getNumberOfPoints(x_size, y_size, z_size);
	for ( core::Size x_index =0; x_index < x_size; x_index+= stride ) {
		for ( core::Size y_index = 0; y_index < y_size; y_index += stride ) {
			for ( core::Size z_index = 0; z_index < z_size; z_index+= stride ) {
				core::grid::CartGrid<core::Real>::GridPt grid_point(x_index,y_index,z_index);
				core::Vector pdb_coords(grid_.coords(grid_point));
				core::Real value = grid_.getValue(grid_point);
				//std::cout << value <<std::endl;
				if ( value >= lower_bound && value <= upper_bound ) {
					std::pair<core::Vector,core::Real> data(pdb_coords,value);
					point_list.push_back(data);
				}
			}
		}
	}
	return point_list;
}

void SingleGrid::grid_to_kin(utility::io::ozstream & out, core::Real min_val, core::Real max_val, core::Size stride)
{
	core::Size x_size(0);
	core::Size y_size(0);
	core::Size z_size(0);
	grid_.getNumberOfPoints(x_size, y_size, z_size);
	for ( core::Size x_index =0; x_index < x_size; x_index +=stride ) {
		for ( core::Size y_index = 0; y_index < y_size; y_index += stride ) {
			for ( core::Size z_index = 0; z_index < z_size; z_index += stride ) {
				core::grid::CartGrid<core::Real>::GridPt grid_point(x_index,y_index,z_index);
				core::Vector box_counter = grid_.coords(grid_point);
				core::Real value = grid_.getValue(grid_point);
				if ( min_val <= value && value <= max_val ) {
					out << '{' << x_index << ' ' << y_index << ' ' << z_index << "}U "
						<< box_counter.x() << ' ' << box_counter.y() << ' ' << box_counter.z() << '\n';
				}
			}
		}
	}
}


bool SingleGrid::is_in_grid(core::conformation::UltraLightResidue const & residue) const
{
	for ( core::Size atom_index = 1; atom_index <= residue.natoms(); ++atom_index ) {
		core::Vector atom_coords = residue[atom_index];
		if ( !grid_.is_in_grid(atom_coords.x(), atom_coords.y(), atom_coords.z()) ) {
			return false;
		}
	}
	return true;
}

bool SingleGrid::is_in_grid(core::conformation::Residue const & residue) const
{
	for ( core::Size atom_index = 1; atom_index <= residue.natoms(); ++atom_index ) {
		core::Vector atom_coords = residue.xyz(atom_index);
		if ( !grid_.is_in_grid(atom_coords.x(), atom_coords.y(), atom_coords.z()) ) {
			return false;
		}
	}
	return true;
}


void SingleGrid::fill_with_value(core::Real value)
{
	grid_.setFullOccupied(value);
}

void SingleGrid::set_sphere(
	core::Vector const & coords,
	core::Real radius,
	core::Real value
){
	set_ring(coords, 0, radius, value);
}

void SingleGrid::set_ring(
	core::Vector const & coords,
	core::Real inner_radius,
	core::Real outer_radius,
	core::Real value
){
	//TR << "Making ring with radii " << inner_radius << "/" << outer_radius << " and value " << value << " at " << coords << std::endl;
	core::Real inner_radius2 = inner_radius*inner_radius;
	core::Real outer_radius2 = outer_radius*outer_radius;
	core::Size x_count(0);
	core::Size y_count(0);
	core::Size z_count(0);
	grid_.getNumberOfPoints(x_count,y_count,z_count);
	debug_assert( x_count >0 && y_count > 0 && z_count > 0 );

	core::Vector vector_radius (outer_radius);
	core::grid::CartGrid<core::Real>::GridPt grid_min = grid_.gridpt(coords - outer_radius);
	core::grid::CartGrid<core::Real>::GridPt grid_max = grid_.gridpt(coords + outer_radius);
	//TR <<"min: " <<grid_min.x() << " "<<grid_min.y() << " " <<grid_min.z() <<std::endl;
	//TR <<"max: "<<grid_max.x() << " "<<grid_max.y() << " " <<grid_max.z() <<std::endl;
	//TR <<"counts: " << x_count <<" "<< y_count << " " << z_count <<std::endl;
	for ( int x_index = std::max(0,grid_min.x()); x_index <= std::min(int(x_count)-1,grid_max.x()); ++x_index ) {
		for ( int y_index = std::max(0,grid_min.y()); y_index <= std::min(int(y_count)-1,grid_max.y()); ++y_index ) {
			for ( int z_index = std::max(0,grid_min.z()); z_index <= std::min(int(z_count)-1,grid_max.z()); ++z_index ) {
				core::grid::CartGrid<core::Real>::GridPt point(x_index, y_index, z_index);
				core::Vector box_center = grid_.coords(point);
				//TR <<x_index << " "<<y_index << " " <<z_index <<std::endl;
				if (
						box_center.distance_squared(coords) >= inner_radius2
						&& box_center.distance_squared(coords) <= outer_radius2
						) {
					grid_.setValue(point, value);
				}
			}
		}
	}
	//TR << "done making sphere "<<std::endl;
}

void SingleGrid::set_distance_sphere_for_atom(core::Real const & atom_shell,core::Vector const & coords,core::Real cutoff)
{
	core::Real cutoff2 = cutoff*cutoff;
	core::Size x_count(0);
	core::Size y_count(0);
	core::Size z_count(0);
	grid_.getNumberOfPoints(x_count,y_count,z_count);
	debug_assert( x_count >0 && y_count > 0 && z_count > 0 );

	core::Vector vector_radius (cutoff);
	core::grid::CartGrid<core::Real>::GridPt grid_min = grid_.gridpt(coords - cutoff);
	core::grid::CartGrid<core::Real>::GridPt grid_max = grid_.gridpt(coords + cutoff);
	for ( int x_index = std::max(0,grid_min.x()); x_index <= std::min(int(x_count)-1,grid_max.x()); ++x_index ) {
		for ( int y_index = std::max(0,grid_min.y()); y_index <= std::min(int(y_count)-1,grid_max.y()); ++y_index ) {
			for ( int z_index = std::max(0,grid_min.z()); z_index <= std::min(int(z_count)-1,grid_max.z()); ++z_index ) {
				core::grid::CartGrid<core::Real>::GridPt point(x_index, y_index, z_index);
				core::Vector box_center(grid_.coords(point));
				core::Real distance2 = box_center.distance_squared(coords);
				if ( distance2 <= cutoff2 ) {
					core::Real distance = sqrt(distance2);
					core::Real current_value = grid_.getValue(point);
					if ( distance - atom_shell <= current_value ) {
						grid_.setValue(point,distance - atom_shell);
					}
				}
			}
		}
	}
}

void SingleGrid::set_score_sphere_for_atom(numeric::interpolation::spline::InterpolatorCOP lj_spline,core::Vector const & coords, core::Real cutoff)
{
	core::Real cutoff2 = cutoff*cutoff;
	core::Size x_count(0);
	core::Size y_count(0);
	core::Size z_count(0);
	grid_.getNumberOfPoints(x_count,y_count,z_count);
	debug_assert( x_count >0 && y_count > 0 && z_count > 0 );

	core::Vector vector_radius (cutoff);
	core::grid::CartGrid<core::Real>::GridPt grid_min = grid_.gridpt(coords - cutoff);
	core::grid::CartGrid<core::Real>::GridPt grid_max = grid_.gridpt(coords + cutoff);
	for ( int x_index = std::max(0,grid_min.x()); x_index <= std::min(int(x_count)-1,grid_max.x()); ++x_index ) {
		for ( int y_index = std::max(0,grid_min.y()); y_index <= std::min(int(y_count)-1,grid_max.y()); ++y_index ) {
			for ( int z_index = std::max(0,grid_min.z()); z_index <= std::min(int(z_count)-1,grid_max.z()); ++z_index ) {
				core::grid::CartGrid<core::Real>::GridPt point(x_index, y_index, z_index);
				core::Vector box_center(grid_.coords(point));
				core::Real distance2 = box_center.distance_squared(coords);
				if ( distance2 <= cutoff2 ) {
					core::Real distance = sqrt(distance2);
					core::Real current_value = grid_.getValue(point);
					core::Real spline_score = 0.0;
					core::Real spline_score_deriv = 0.0;

					lj_spline->interpolate(distance,spline_score,spline_score_deriv);
					if ( spline_score <= current_value ) {
						grid_.setValue(point,spline_score);
					}

				}
			}
		}
	}

}

void SingleGrid::diffuse_ring(core::Vector const & coords, core::Real radius, core::Real width, core::Real magnitude)
{
	core::Real radius_squared = radius*radius;
	core::Real half_width = width/2;
	core::Real half_width_squared = half_width*half_width;
	core::Size x_count(0);
	core::Size y_count(0);
	core::Size z_count(0);
	grid_.getNumberOfPoints(x_count,y_count,z_count);
	debug_assert( x_count >0 && y_count > 0 && z_count > 0 );

	core::Vector vector_radius (radius);
	core::grid::CartGrid<core::Real>::GridPt grid_min = grid_.gridpt(coords - radius);
	core::grid::CartGrid<core::Real>::GridPt grid_max = grid_.gridpt(coords + radius);
	for ( int x_index = std::max(0,grid_min.x()); x_index <= std::min(int(x_count)-1,grid_max.x()); ++x_index ) {
		for ( int y_index = std::max(0,grid_min.y()); y_index <= std::min(int(y_count)-1,grid_max.y()); ++y_index ) {
			for ( int z_index = std::max(0,grid_min.z()); z_index <= std::min(int(z_count)-1,grid_max.z()); ++z_index ) {
				core::grid::CartGrid<core::Real>::GridPt point(x_index,y_index,z_index);
				core::Vector origin = grid_.coords(point);
				core::Real distance_squared = origin.distance_squared(coords);
				//core::Real rad_distance_squared = std::abs(radius2-distance_squared);

				//if((distance_squared < radius2 + half_width) || (distance_squared > radius2-half_width))
				//if(distance_squared <radius2)
				if ( (distance_squared >radius_squared-half_width_squared && distance_squared < radius_squared)|| (distance_squared <radius_squared+half_width_squared && distance_squared > radius_squared) ) {
					//TR << "inserting value " <<std::endl;
					grid_.setValue(point,magnitude);
				} else {

					//TR<<  "d^2 is " << distance_squared <<std::endl;
					core::Real distance_from_ring = 0.0;
					if ( distance_squared < radius_squared ) {
						distance_from_ring = (radius_squared-half_width_squared)-distance_squared;
					} else {
						distance_from_ring = distance_squared-(radius_squared+half_width_squared);
					}

					if ( distance_from_ring <= std::abs(magnitude) ) {
						grid_.setValue(point,magnitude+distance_from_ring);
					} else { }

				}

			}
		}
	}
}

void SingleGrid::set_point(core::Vector const & coords, core::Real value)
{
	grid_.setValue(coords,value);
}

void SingleGrid::dump_BRIX(std::string const & prefix) const
{
	//std::string grid_type_string(qsar::qsarTypeManager::name_from_qsar_type(type_));
	grid_.write_to_BRIX(prefix + type_ + ".omap");
}

/// @brief Print a brief summary about this grid to the provided output stream
void SingleGrid::show( std::ostream & out ) const {
	out << "SingleGrid of type " << type_ << " for chain " << chain_ << " with center " << center_ << " between " << grid_.getBase() << " and " << grid_.getTop() << std::endl;
}

}
}
}
