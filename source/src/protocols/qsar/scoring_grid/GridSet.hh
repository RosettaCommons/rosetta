// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/qsar/scoring_grid/GridSet.hh
/// @brief A set of related grids
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_protocols_qsar_scoring_grid_GridSet_hh
#define INCLUDED_protocols_qsar_scoring_grid_GridSet_hh

#include <protocols/qsar/scoring_grid/GridSet.fwd.hh>

#include <protocols/qsar/scoring_grid/GridBase.fwd.hh>
#include <protocols/qsar/qsarMap.fwd.hh>
#include <protocols/qsar/scoring_grid/ScoreNormalization.fwd.hh>

#include <core/conformation/UltraLightResidue.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/types.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/json_spirit/json_spirit_writer.h>
#include <utility/json_spirit/json_spirit_reader.h>

#include <map>

namespace protocols {
namespace qsar {
namespace scoring_grid {

///@brief A set of related grids
class GridSet : public utility::pointer::ReferenceCount {

public:

	typedef std::map< std::string, core::Real > ScoreMap;

	GridSet();
	GridSet(GridSet const & src);

	virtual ~GridSet();

	GridSetOP
	clone() const;

	/// @brief set width (must be done before initialization)
	void width(core::Real width) { width_ = width; }
	/// @brief set resolution (must be done before initialization)
	void resolution(core::Real resolution) { resolution_ = resolution; }
	/// @brief set ligand chain (must be done before initialization)
	void chain(char chain) { chain_ = chain; }

	/// @brief get grid width
	core::Real width() const { return width_; }
	/// @brief get grid resoultion
	core::Real resolution() const { return resolution_; }
	/// @brief get ligand chain
	char chain() const { return chain_; }

	/// @brief set normalization function
	void set_normalization_function(std::string norm_function_name);
	/// @brief is normalization applied during scoring?
	bool is_normalization_enabled() const;
	/// @brief set the qsar_map
	void set_qsar_map(qsarMapCOP qsar_map);
	/// @brief is a qsar map attached to the grid manager?
	bool is_qsar_map_attached() const;

	/// @brief Reinitialize the included grids with the given information
	void
	reinitialize( core::pose::Pose const & pose, core::Vector const & center);

	GridBaseCOP get_grid( std::string const & name ) const;

	/// @brief Add a grid to the GridSet under the given name.
	void add_grid( std::string const & name, GridBaseOP grid, core::Real weight=1.0 );

	/// @brief make a new grid from grid tag, and insert it into the grid manager
	void make_new_grid(utility::tag::TagCOP tag);

	/// @brief Returns true if we already have a grid with the given name.
	bool has_grid( std::string const & name ) const;

	/// @brief get a list of grid names
	utility::vector1<std::string> get_grid_names() const;

	/// @brief return the number of grids
	core::Size size() const;

	///@brief return the average score of multiple residues on the grid
	core::Real average_score(utility::vector1<core::conformation::UltraLightResidue> & residues) const;

	///@brief return the total score of a residue on the grid
	core::Real total_score(core::conformation::UltraLightResidue const & residue) const;
	/// @brief return the total score of a residue on the grid
	core::Real total_score(core::conformation::Residue const & residue) const;
	/// @brief return the total score of a chain on the grid
	core::Real total_score(core::pose::Pose const & pose, core::Size const chain_id) const;
	/// @brief return the total score for a set of residues on the grid
	core::Real total_score(core::pose::Pose const & pose, utility::vector1< core::Size > const & residues) const;

	// @brief Get a map of (unweighted, unnormalized) scores for each grid type
	ScoreMap grid_scores( core::conformation::Residue const & residue ) const;
	// @brief Get a map of (unweighted, unnormalized) scores for each grid type
	ScoreMap grid_scores( core::pose::Pose const & pose, utility::vector1< core::Size > const & residues ) const;

	/// @brief get a map of (unweighted, unnormalized) scoring terms and scores for each term given a residue and atom number
	ScoreMap atom_score(core::pose::Pose const & pose, core::conformation::Residue const & residue, core::Size atomindex ) const;

	///@brief check if all atoms in all the ligands are in grid
	bool is_in_grid(utility::vector1<core::conformation::UltraLightResidue> const & residues) const;
	///@brief check to see if all atoms in the ligand are in the grid
	bool is_in_grid(core::conformation::UltraLightResidue const & residue) const;
	///@brief check to see if all atoms in the ligand are in the grid
	bool is_in_grid(core::conformation::Residue const & residue) const;

	/// @brief Return a string representing the settings of this GridSet
	/// @details This should encapsulate everything that doesn't change with a call to reinitialize()
	/// The string is not meant to be interpretable.
	std::string hash_fingerprint() const;

	/// @brief write all grids out using the BRIX format
	void write_grids(std::string prefix) const;

	/// @brief serialize the current map to a JSON object.
	utility::json_spirit::Value serialize() const;
	/// @brief deserialize the JSON object to a map.
	void deserialize(utility::json_spirit::mArray data);

private:

	/// @brief Return a single weighted, normalized score for all items
	/// Templated to allow for common implementation between Residue and UltraLightResidue
	/// @details Yes, raw pointer. This is for speed reasons (don't have to copy objects)
	/// and because this is a class-internal implementation detail. Do not publically expose this interface.
	template< class GridScorable >
	core::Real
	total_score( utility::vector1< GridScorable const * > const & items ) const;

	/// @brief Normalize a score based on the passed items;
	/// @details Yes, raw pointer. This is to interface with the other functions. Do not publically expose this interface.
	core::Real
	normalize( core::Real total_score, utility::vector1< core::conformation::UltraLightResidue const* > const & items ) const;

	/// @brief Normalize a score based on the passed items;
	/// @details Yes, raw pointer. This is to interface with the other functions. Do not publically expose this interface.
	core::Real
	normalize( core::Real total_score, utility::vector1< core::conformation::Residue const* > const & items ) const;

	/// @brief Return a map of summed unweighted, unormalized scores for each grid.
	/// Templated to allow for common implementation between Residue and UltraLightResidue
	/// @details Yes, raw pointer. This is for speed reasons (don't have to copy objects)
	/// and because this is a class-internal implementation detail. Do not publically expose this interface.
	template< class GridScorable >
	ScoreMap
	grid_scores( utility::vector1< GridScorable const* > const & items ) const;

private:

	typedef std::map<std::string,GridBaseOP> MappingType;

	MappingType grids_;
	std::map< std::string, core::Real > grid_weights_;

	core::Real width_ = 40;
	core::Real resolution_ = 0.25;
	char chain_ = 'X';

	qsar::qsarMapCOP qsar_map_;
	ScoreNormalizationCOP norm_function_;
};


} //protocols
} //qsar
} //scoring_grid



#endif //INCLUDED_protocols_qsar_scoring_grid_GridSet_hh





