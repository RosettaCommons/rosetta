// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/qsar/scoring_grid/ElectrostaticMetaGrid.hh
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_qsar_scoring_grid_ElectrostaticMetaGrid_HH
#define INCLUDED_protocols_qsar_scoring_grid_ElectrostaticMetaGrid_HH

#include <protocols/qsar/scoring_grid/ElectroStaticMetaGrid.fwd.hh>
#include <protocols/qsar/scoring_grid/GridBase.hh>
#include <protocols/qsar/scoring_grid/ChargeGrid.hh>

#include <map>

namespace protocols {
namespace qsar {
namespace scoring_grid {

class ElectroStaticMetaGrid : public GridBase {

public:

	ElectroStaticMetaGrid();
	virtual ~ElectroStaticMetaGrid() {}

	virtual void initialize(core::Vector const & center, core::Real width, core::Real resolution);
	/// @brief populate the grid with values based on a passed pose
	virtual void refresh(core::pose::Pose const & pose, core::Vector const & center, core::Size const & ligand_chain_id_to_exclude);
	/// @brief populate the grid with values based on a passed pose
	virtual void refresh(core::pose::Pose const & pose, core::Vector const & center,utility::vector1<core::Size> ligand_chain_ids_to_exclude);
	/// @brief populate the grid with values based on a passed pose
	virtual void refresh(core::pose::Pose const & pose, core::Vector const & center);
	/// @setup a grid based on RosettaScripts input
	virtual void parse_my_tag(utility::tag::TagPtr const tag);
	/// @brief return the current score of a residue using the current grid
	virtual core::Real score(core::conformation::Residue const & residue, core::Real const max_score, qsarMapOP qsar_map);
	/// @brief get the type of the grid
	virtual std::string get_type();
	/// @brief set the chain the grid applies to
	virtual void set_chain(char chain);

	virtual void dump_BRIX(std::string const & prefix);
private:
	void add_new_charge_grid(core::Real const & charge);

private:
	utility::vector1<core::Real> charges_;
	std::map<core::Real, ChargeGrid> charge_grid_map_;
	std::string type_;
	char chain_;
	core::Real weight_;
};

}
}
}

#endif /* ELECTROSTATICMETAGRID_HH_ */
