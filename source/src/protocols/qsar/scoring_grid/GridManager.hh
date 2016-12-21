// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/qsar/scoring_grid/GridManager.hh
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_qsar_scoring_grid_GridManager_hh
#define INCLUDED_protocols_qsar_scoring_grid_GridManager_hh


#include <protocols/qsar/scoring_grid/GridManager.fwd.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <protocols/qsar/qsarMap.fwd.hh>
#include <protocols/qsar/scoring_grid/GridBase.fwd.hh>
#include <protocols/qsar/scoring_grid/ScoreNormalization.fwd.hh>

#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/UltraLightResidue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/json_spirit/json_spirit_writer.h>
#include <utility/json_spirit/json_spirit_reader.h>

// C++ Headers
#include <string>
#include <map>

namespace protocols {
namespace qsar {
namespace scoring_grid {

typedef std::map<std::string,GridBaseOP> GridMap;
typedef std::map<std::string,core::Real> ScoreMap;

class GridManager : public utility::SingletonBase< GridManager >
{
public:
	friend class utility::SingletonBase< GridManager >;

	/// @brief reset GridManager to the default settings
	void reset();
	/// @brief set width (must be done before initialization)
	void set_width(core::Real width);
	/// @brief set resolution (must be done before initialization)
	void set_resolution(core::Real resolution);
	/// @brief set normalization function
	void set_normalization_function(std::string norm_function_name);
	/// @brief set ligand chain (must be done before initialization)
	void set_chain(char chain);
	/// @brief make a new grid from grid tag, and insert it into the grid manager
	void make_new_grid(utility::tag::TagCOP tag);
	/// @brief insert a grid pointer into the grid manager
	void insert_grid(std::string const & name,GridBaseOP const grid);
	/// @brief set the qsar_map
	void set_qsar_map(qsarMapOP qsar_map);
	/// @brief is a qsar map attached to the grid manager?
	bool is_qsar_map_attached();
	/// @brief is normalization applied during scoring?
	bool is_normalization_enabled();
	/// @brief given a grid type, return a pointer to the grid
	GridBaseOP get_grid(std::string const & grid_name);
	/// @brief get a list of grid names
	utility::vector1<std::string> get_grid_names();

	///@brief return the ideal average score of multiple residues on the grid
	core::Real ideal_score(utility::vector1<core::conformation::UltraLightResidue> & residues);
	///@brief return the ideal score of a residue on the grid
	core::Real ideal_score(core::conformation::UltraLightResidue const & residue);
	///@brief return the average score of multiple residues on the grid 
	core::Real total_score(utility::vector1<core::conformation::UltraLightResidue> & residues);
	///@brief return the total score of a residue on the grid
	core::Real total_score(core::conformation::UltraLightResidue const & residue);
	/// @brief return the total score of a residue on the grid
	core::Real total_score(core::conformation::Residue const & residue);
	/// @brief return the total score of a chain on the grid
	core::Real total_score(core::pose::Pose const & pose, core::Size const chain_id);
	/// @brief get a map of scoring terms and scores for each term given a residue and atom number
	std::map<std::string,core::Real> atom_score(core::pose::Pose const & pose, core::conformation::Residue const & residue, core::Size atomindex );
	/// @brief recalculate all grids for a pose.  This must be called if the backbone or sidechain conformations change!
	void update_grids(core::pose::Pose const & pose, core::Vector const & center, core::Size const & ligand_chain_id_to_exclude);
	/// @brief recalculate all grids for a pose.  This must be called if the backbone or sidechain conformations change!
	void update_grids(core::pose::Pose const & pose, core::Vector const & center,utility::vector1<core::Size> ligand_chain_ids_to_exclude);
	///@brief recalculate all grids for a pose.  This must be called if the backbone or sidechain conformations change!
	void update_grids(core::pose::Pose const & pose, core::Vector const & center, bool multi = false);
//	void reset_grids();
	///@brief initialize all grids and fill with 0s given a center point
	void initialize_all_grids(core::Vector const & center);
	/// @brief return the number of grids in the manager
	core::Size size();
	/// @brief get a map of cached scores
	ScoreMap get_cached_scores();
	/// @brief append all cached scores to a current job
	void append_cached_scores(jd2::JobOP job);
	/// @brief write all grids out using the BRIX format
	void write_grids(std::string prefix);
	///@brief check if all atoms in all the ligands are in grid 
	bool is_in_grid(utility::vector1<core::conformation::UltraLightResidue> const & residues);
    ///@brief check to see if all atoms in the ligand are in the grid
    bool is_in_grid(core::conformation::UltraLightResidue const & residue);
    ///@brief check to see if all atoms in the ligand are in the grid
    bool is_in_grid(core::conformation::Residue const & residue);

private:

	GridManager();
	GridManager(GridManager const &) = delete;
	GridManager const & operator = (GridManager const & ) = delete;

	/// @brief serialize the current map to a JSON object.  There is no public interface for this because the grid manager takes care of it on its own
	utility::json_spirit::Value serialize();
	/// @brief deserialize the JSON object to a map.  There is no public interface for this because the grid manager takes care of it on its own
	void deserialize(utility::json_spirit::mArray data);

private:

	std::map<std::string,GridMap> grid_map_cache_;

	std::map<std::string,core::Real> grid_weights_;

	GridMap grid_map_;
	ScoreMap score_map_;
	std::string last_tag_;
	core::Real width_;
	core::Real resolution_;
	qsar::qsarMapOP qsar_map_;
	bool initialized_;
	char chain_;
	ScoreNormalizationOP norm_function_;

};

}
}
}

#endif /* GRIDMANAGER_HH_ */
