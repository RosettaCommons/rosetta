// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/qsar/scoring_grid/SingleGrid.hh
/// @author Sam DeLuca


#ifndef INCLUDED_protocols_qsar_scoring_grid_SingleGrid_hh
#define INCLUDED_protocols_qsar_scoring_grid_SingleGrid_hh

#include <protocols/qsar/scoring_grid/SingleGrid.fwd.hh>
#include <protocols/qsar/scoring_grid/GridBase.hh>
#include <protocols/qsar/qsarMap.fwd.hh>

#include <core/grid/CartGrid.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.fwd.hh>

#include <numeric/interpolation/spline/SplineGenerator.hh>
#include <numeric/xyzVector.hh>

#include <list>

#ifdef WIN32
#include <protocols/qsar/qsarMap.hh>
#endif


namespace protocols {
namespace qsar {
namespace scoring_grid {

class SingleGrid : public GridBase
{
public:

	SingleGrid(std::string type);
	virtual ~SingleGrid();
	/// @brief initialize a grid of zeros with a given centerpoint, width and resolution (in angstroms).
	virtual void initialize(core::Vector const & center, core::Real width, core::Real resolution);
	/// @brief set the chain around which to calculate the grid
	virtual void set_chain(char chain);
	/// @brief get the chain around which the grid is calculated
	char get_chain();
	/// @brief populate the grid with values based on a passed pose
	virtual void refresh(core::pose::Pose const & pose, core::Vector const & center, core::Size const & ligand_chain_id_to_exclude)=0;
	/// @brief populate the grid with values based on a passed pose
	virtual void refresh(core::pose::Pose const & pose, core::Vector const & center,utility::vector1<core::Size> ligand_chain_ids_to_exclude)=0;
	/// @brief populate the grid with values based on a passed pose
	virtual void refresh(core::pose::Pose const & pose, core::Vector const & center)=0;
	// virtual void reset(){};
	/// @brief setup a grid based on RosettaScripts input
	virtual void parse_my_tag(utility::tag::TagCOP tag)=0;
	/// @brief serialize the SingleGrid to a json_spirit object
	virtual utility::json_spirit::Value serialize();
	/// @brief deserialize a json_spirit object to a SingleGrid
	virtual void deserialize(utility::json_spirit::mObject data);
	/// @brief return a constant reference to the grid
	core::grid::CartGrid<core::Real> const &  get_grid();
	/// @brief set the grid type
	void set_type(std::string type);
	/// @brief return the grids type
	virtual std::string get_type();
	/// @brief set the center of the grid
	void set_center(core::Vector center);
	/// @brief get the center of the grid
	core::Vector get_center();
	/// @brief get the max score value in the grid
	core::Real get_min_value() const;
	/// @brief get the minimum score value in the grid
	core::Real get_max_value() const;
	/// @brief get the value of a single point in the grid based on pdb coordinates
	core::Real get_point(core::Real x, core::Real y, core::Real z);
	/// @brief get the value of a single point in the grid based on pdb coordinates
	core::Real get_point(core::Vector coords);
	/// @brief get dimensions of the grid
	numeric::xyzVector<core::Size> get_dimensions();
	/// @brief get the pdb coordinates based on grid point coordinates
	core::Vector get_pdb_coords(int x, int y, int z);
	/// @brief get the pdb coordinates based on grid point coordinates
	core::Vector get_pdb_coords(core::grid::CartGrid<core::Real>::GridPt gridpt);
	/// @brief return the current score of an UltraLightResidue using the current grid
	virtual core::Real score(core::conformation::UltraLightResidue const & residue, core::Real const max_score, qsarMapOP qsar_map);
	/// @brief return the current score of an atom using the current grid
	virtual core::Real atom_score(core::conformation::UltraLightResidue const & residue, core::Size atomno, qsarMapOP qsar_map);
	/// @brief return the current score of a residue using the current grid
	virtual core::Real score(core::conformation::Residue const & residue, core::Real const max_score, qsarMapOP qsar_map);
	/// @brief return the current score of an atom using the current grid
	virtual core::Real atom_score(core::conformation::Residue const & residue, core::Size atomno, qsarMapOP qsar_map);
	void grid_to_kin(utility::io::ozstream & out, core::Real min_val, core::Real max_val, core::Size stride);
	//void grid_rotamer_trials(core::pose::Pose &  pose, core::Size residue_id, int const min_score);
	/// @brief check to see if residue is in grid
	virtual bool is_in_grid(core::conformation::UltraLightResidue const & residue);
	/// @brief check to see if residue is in grid
	virtual bool is_in_grid(core::conformation::Residue const & residue);

	std::list<std::pair<core::Vector, core::Real> > get_point_value_list_within_range(core::Real lower_bound, core::Real upper_bound,core::Size stride);

	virtual void dump_BRIX(std::string const & prefix);

	//Various mathematical functions for assigning values to the grid go here
	void set_sphere(core::Vector const & coords, core::Real radius, core::Real value);
	void set_ring(
		core::Vector const & coords,
		core::Real inner_radius,
		core::Real outer_radius,
		core::Real value
	);
	void diffuse_ring(core::Vector const & coords, core::Real radius, core::Real width, core::Real magnitude);
	//void set_distance_sphere( core::Vector const & coords,core::Real cutoff);
	void set_point(core::Vector const & coords, core::Real value);
	void set_distance_sphere_for_atom(core::Real const & atom_shell, core::Vector const & coords,core::Real cutoff);
	void set_score_sphere_for_atom(numeric::interpolation::spline::InterpolatorOP lj_spline,core::Vector const & coords, core::Real cutoff);
	/// @fill the entire grid with some value
	void fill_with_value(core::Real);
private:
	core::grid::CartGrid<core::Real> grid_;
	std::string type_;
	core::Vector center_;
	char chain_;

};


}
}
}

#endif /* GRIDBASE_HH_ */
