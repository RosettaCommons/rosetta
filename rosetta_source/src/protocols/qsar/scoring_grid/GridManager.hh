// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   /git/src/protocols/qsar/scoring_grid/GridManager.hh
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_qsar_scoring_grid_GridManager_hh
#define INCLUDED_protocols_qsar_scoring_grid_GridManager_hh

#include <map>
#include <utility/pointer/ReferenceCount.hh>
#include <protocols/qsar/qsarMap.hh>
#include <protocols/qsar/scoring_grid/GridBase.hh>
#include <protocols/qsar/scoring_grid/GridManager.fwd.hh>
#include <protocols/jd2/Job.fwd.hh>

namespace protocols {
namespace qsar {
namespace scoring_grid {

class GridManager: public utility::pointer::ReferenceCount
{
public:
	static GridManager* get_instance();
	///@brief set width and resolution (must be done before initialization
	void set_dimensions(core::Real width, core::Real resolution);
	///@brief make a new grid given the name of a grid type, and insert it into the grid manager
	void make_new_grid(std::string grid_name);
	///@brief insert a grid pointer into the grid manager
	void insert_grid(GridBaseOP grid);
	///@brief set the qsar_map
	void set_qsar_map(qsarMapOP qsar_map);
	///@brief is a qsar map attached to the grid manager?
	bool is_qsar_map_attached();
	///@brief given a grid type, return a pointer to the grid
	GridBaseOP get_grid(std::string const & grid_name);
	///@brief get a list of grid names
	utility::vector1<std::string> get_grid_names();
	///@brief return the total score of a residue on the grid
	core::Real total_score(core::conformation::Residue const & residue);
	///@brief return the total score of a chain on the grid
	core::Real total_score(core::pose::Pose const & pose, core::Size const chain_id);
	///@brief recalculate all grids for a pose.  This must be called if the backbone or sidechain conformations change!
	void update_grids(core::pose::Pose const & pose, core::Vector const & center, core::Size const & ligand_chain_id_to_exclude);
	///@brief recalculate all grids for a pose.  This must be called if the backbone or sidechain conformations change!
	void update_grids(core::pose::Pose const & pose, core::Vector const & center,utility::vector1<core::Size> ligand_chain_ids_to_exclude);
	///@brief recalculate all grids for a pose.  This must be called if the backbone or sidechain conformations change!
	void update_grids(core::pose::Pose const & pose, core::Vector const & center);
	///@brief recalculate all grids for a pose, but only if the tag has changed
	void update_grids(core::pose::Pose const & pose, core::Vector const & center,std::string const & tag);
	///@brief initialize all grids and fill with 0s given a center point
	void initialize_all_grids(core::Vector const & center);
	///@brief return the number of grids in the manager
	core::Size size();
	///@brief get a map of cached scores
	std::map<std::string, core::Real> get_cached_scores();
	///@brief append all cached scores to a current job
	void append_cached_scores(jd2::JobOP job);
	///@brief write all grids out using the BRIX format
	void write_grids(std::string prefix);

private:
	GridManager();
	GridManager(GridManager const &);
	GridManager const & operator = (GridManager const & );
	//GridManager(core::Real width,core::Real resolution);
	static GridManager * instance_;

	std::map<std::string,GridBaseOP> grid_map_;
	std::map<std::string,core::Real> score_map_;
	std::map<std::string,core::Real> weight_map_;
	std::string last_tag_;
	core::Real width_;
	core::Real resolution_;
	qsar::qsarMapOP qsar_map_;
	bool initialized_;
};

}
}
}

#endif /* GRIDMANAGER_HH_ */
