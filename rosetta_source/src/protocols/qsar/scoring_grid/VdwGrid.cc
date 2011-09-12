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
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/AtomGraph.hh>
#include <core/conformation/find_neighbors.hh>
#include <basic/database/open.hh>
#include <numeric/interpolation/util.hh>
#include <utility/tag/Tag.hh>

namespace protocols {
namespace qsar {
namespace scoring_grid {

std::string VdwGridCreator::keyname() const
{
	return VdwGridCreator::grid_name();
}

GridBaseOP VdwGridCreator::create_grid(utility::tag::TagPtr const tag) const
{
	if (tag->hasOption("weight")){
		return new VdwGrid( tag->getOption<core::Real>("weight") );
	}else{
		return new VdwGrid();
	}
	// This is impossible
	return NULL;
}

std::string VdwGridCreator::grid_name()
{
	return "VdwGrid";
}

VdwGrid::VdwGrid() : GridBase("VdwGrid",0.0)
{
	std::string lj_file(basic::database::full_name("qsar/lj_table.txt"));
	lj_spline_ = numeric::interpolation::spline_from_file(lj_file,0.01);
}

VdwGrid::VdwGrid(core::Real weight) : GridBase ("VdwGrid",weight)
{
	std::string lj_file(basic::database::full_name("qsar/lj_table.txt"));
	lj_spline_ = numeric::interpolation::spline_from_file(lj_file,0.01);
}

void VdwGrid::refresh(core::pose::Pose const & pose, core::Vector const &  )
{
	core::conformation::AtomGraphOP atom_graph(new core::conformation::AtomGraph);
	core::PointPosition query_point(0.0,0.0,0.0);
	core::Size index_vertex(core::conformation::annotated_atom_graph_from_conformation(pose.conformation(),atom_graph,query_point));
	numeric::xyzVector<core::Size> grid_size(this->get_dimensions());
	for(core::Size x_index = 0; x_index < grid_size.x(); ++x_index)
	{
		for(core::Size y_index = 0; y_index < grid_size.y(); ++y_index)
		{
			for(core::Size z_index = 0; z_index < grid_size.z(); ++z_index)
			{
				core::Vector current_point(x_index,y_index,z_index);
				atom_graph->drop_all_edges();
				atom_graph->get_vertex(index_vertex).data().xyz() = current_point;

				core::Size neighbor_id = core::conformation::get_nearest_neighbor<core::conformation::AtomGraphVertexData,core::conformation::AtomGraphEdgeData>(atom_graph,index_vertex,20.0);
				if(neighbor_id == index_vertex)
				{
					this->set_point(current_point,20.0);
					continue;
				}
				core::conformation::AtomGraphVertexData neighbor_node_data  = atom_graph->get_vertex(neighbor_id).data();
				core::conformation::AtomGraphEdgeData neighbor_edge_data = atom_graph->get_edge(index_vertex,neighbor_id)->data();

				core::Real dsq_to_center = neighbor_edge_data.dsq();
				core::Real radius_squared = neighbor_node_data.atom_radius_squared();
				core::Real real_dsq = dsq_to_center - radius_squared;
				if(real_dsq < 0.0)
				{
					this->set_point(current_point,0.0);
					continue;
				}
				else
				{
					this->set_point(current_point,std::sqrt(real_dsq));
				}
			}
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
		qsarPointOP qsar_info(qsar_map->get_point(atom_index,"VdwGrid"));
		if(qsar_info != 0)
		{
			core::Vector const & atom_coord(residue.xyz(atom_index));
			core::Real const & radius(residue.atom_type(atom_index).lj_radius());
			if(this->get_grid().is_in_grid(atom_coord.x(),atom_coord.y(),atom_coord.z()))
			{
				//core::grid::CartGrid<core::Real>::GridPt grid_point(atom_coord.x(),atom_coord.y(),atom_coord.z());
				core::Real max_radius = this->get_point(atom_coord.x(),atom_coord.y(),atom_coord.z());
				core::Real spline_score = 0.0;
				core::Real spline_score_deriv = 0.0;
				interpolator->interpolate(max_radius-radius,spline_score,spline_score_deriv);

				score += spline_score;
			}
		}
	}
	return score;
}

}
}
}
