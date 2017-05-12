// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/ligand_docking/scoring_grid/SolvationGrid.hh
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_qsar_scoring_grid_SolvationGrid_hh
#define INCLUDED_protocols_qsar_scoring_grid_SolvationGrid_hh

#include <protocols/qsar/scoring_grid/SolvationGrid.fwd.hh>
#include <protocols/qsar/scoring_grid/SingleGrid.hh>

#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#include <list>

namespace protocols {
namespace qsar {
namespace scoring_grid {

class SolvationGrid : public SingleGrid
{
public:

	SolvationGrid();
	~SolvationGrid() override;
	/// @brief Make a copy of the grid, respecting the subclassing.
	GridBaseOP clone() const override;
	void refresh(core::pose::Pose const & pose, core::Vector const & center, core::Size const & ) override;
	void refresh(core::pose::Pose const & pose, core::Vector const & center) override;
	void refresh(core::pose::Pose const & pose, core::Vector const & center, utility::vector1<core::Size> ) override;

	void parse_my_tag(utility::tag::TagCOP tag) override;

	/// @brief serialize the Interpolator to a json_spirit object
	utility::json_spirit::Value serialize() const override;
	/// @brief deserialize a json_spirit object to a Interpolator
	void deserialize(utility::json_spirit::mObject data) override;

	void set_probe_atom_type(core::ShortSize const & atom_type);

	static std::string grid_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	core::ShortSize probe_atom_type_;



};

}
}
}

#endif
