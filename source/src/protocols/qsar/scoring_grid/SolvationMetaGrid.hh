// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/ligand_docking/scoring_grid/SolvationMetaGrid.hh
/// @author Sam DeLuca


#ifndef INCLUDED_protocols_qsar_scoring_grid_SolvationMetaGrid_hh
#define INCLUDED_protocols_qsar_scoring_grid_SolvationMetaGrid_hh

#include <protocols/qsar/scoring_grid/SolvationMetaGrid.fwd.hh>

#include <protocols/qsar/scoring_grid/SolvationGrid.hh>
#include <protocols/qsar/scoring_grid/GridBase.hh>

#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#include <map>

namespace protocols {
namespace qsar {
namespace scoring_grid {

class SolvationMetaGrid : public GridBase
{
public:
	SolvationMetaGrid();
	SolvationMetaGrid( SolvationMetaGrid const & other);
	~SolvationMetaGrid() override;
	/// @brief Make a copy of the grid, respecting the subclassing.
	GridBaseOP clone() const override;
	/// @brief initialize a grid of zeros with a given centerpoint, width and resolution (in angstroms).
	void initialize(core::Vector const & center, core::Real width, core::Real resolution) override;
	/// @brief populate the grid with values based on a passed pose
	void refresh(core::pose::Pose const & pose, core::Vector const & center, core::Size const & ligand_chain_id_to_exclude) override;
	/// @brief populate the grid with values based on a passed pose
	void refresh(core::pose::Pose const & pose, core::Vector const & center,utility::vector1<core::Size> ligand_chain_ids_to_exclude) override;
	/// @brief populate the grid with values based on a passed pose
	void refresh(core::pose::Pose const & pose, core::Vector const & center) override;
	/// @setup a grid based on RosettaScripts input
	void parse_my_tag(utility::tag::TagCOP tag) override;
	/// @brief return the current score of an UltraLightResidue using the current grid
	core::Real score(core::conformation::UltraLightResidue const & residue, core::Real const max_score, qsarMapOP qsar_map) const override;
	/// @brief return the current score of an atom using the current grid
	core::Real atom_score(core::conformation::UltraLightResidue const & residue, core::Size atomno, qsarMapOP qsar_map) const override;
	/// @brief return the current score of a residue using the current grid
	core::Real score(core::conformation::Residue const & residue, core::Real const max_score, qsarMapOP qsar_map) const override;
	/// @brief return the current score of an atom using the current grid
	core::Real atom_score(core::conformation::Residue const & residue, core::Size atomno, qsarMapOP qsar_map) const override;
	/// @brief get the type of the grid
	std::string get_type() const override;
	/// @brief set the chain the grid applies to
	void set_chain(char chain) override;
	/// @brief output a BRIX formatted grid.  This really does not work well but is being left for legacy purposes
	void dump_BRIX(std::string const & prefix) const override;
	/// @brief Serialize the GridBase object into a json_spirit Value
	utility::json_spirit::Value serialize() const override;
	/// @brief deserialize a json spirit Value into a GridBase object
	void deserialize(utility::json_spirit::mObject data) override;
	/// @brief determine if all residue atoms are in a grid
	bool is_in_grid(core::conformation::UltraLightResidue const & residue) const override;
	/// @brief determine if all residue atoms are in a grid
	bool is_in_grid(core::conformation::Residue const & residue) const override;

	static std::string grid_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief Print a brief summary about this grid to the provided output stream
	void show( std::ostream & out ) const override;

private:
	std::string type_;
	std::map<core::ShortSize,SingleGridOP> grid_map_;
};

}
}
}


#endif /* SOLVATIONMETAGRID_HH_ */
