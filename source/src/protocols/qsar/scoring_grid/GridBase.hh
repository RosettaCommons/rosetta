// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/qsar/scoring_grid/GridBase.hh
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_qsar_scoring_grid_GridBase_HH
#define INCLUDED_protocols_qsar_scoring_grid_GridBase_HH

#include <protocols/qsar/scoring_grid/GridBase.fwd.hh>
#include <protocols/qsar/qsarMap.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/UltraLightResidue.fwd.hh>
#include <utility/json_spirit/json_spirit_value.h>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

namespace protocols {
namespace qsar {
namespace scoring_grid {

class GridBase:  public utility::pointer::ReferenceCount
{
public:
	GridBase() {}
	virtual ~GridBase() {}

	/// @brief Make a copy of the grid, respecting the subclassing.
	/// @details Note that due to the heavy use of modification to reset positions
	/// you need to be sure to do deep copying on anything that can be changed
	/// with refresh/initialize.
	virtual GridBaseOP clone() const =0;

	/// @brief initialize a grid of zeros with a given centerpoint, width and resolution (in angstroms).
	virtual void initialize(core::Vector const & center, core::Real width, core::Real resolution)=0;
	/// @brief populate the grid with values based on a passed pose
	virtual void refresh(core::pose::Pose const & pose, core::Vector const & center, core::Size const & ligand_chain_id_to_exclude)=0;
	/// @brief populate the grid with values based on a passed pose
	virtual void refresh(core::pose::Pose const & pose, core::Vector const & center,utility::vector1<core::Size> ligand_chain_ids_to_exclude)=0;
	/// @brief populate the grid with values based on a passed pose
	virtual void refresh(core::pose::Pose const & pose, core::Vector const & center)=0;
	// virtual void reset(){};
	/// @setup a grid based on RosettaScripts input
	virtual void parse_my_tag(utility::tag::TagCOP tag)=0;

	/// @brief return the current scoer of an UltraLightResidue using the current grid
	virtual core::Real score(core::conformation::UltraLightResidue const & residue, core::Real const max_score, qsarMapOP qsar_map) const = 0;
	/// @brief return the current score of an atom using the current grid
	virtual core::Real atom_score(core::conformation::UltraLightResidue const & residue, core::Size atomno, qsarMapOP qsar_map) const = 0;
	/// @brief return the current score of a residue using the current grid
	virtual core::Real score(core::conformation::Residue const & residue, core::Real const max_score, qsarMapOP qsar_map) const = 0;
	/// @brief return the current score of an atom using the current grid
	virtual core::Real atom_score(core::conformation::Residue const & residue, core::Size atomno, qsarMapOP qsar_map) const = 0;
	/// @brief get the type of the grid
	virtual std::string get_type() const = 0;
	/// @brief set the chain the grid applies to
	virtual void set_chain(char chain) = 0;
	/// @brief output a BRIX formatted grid.  This really does not work well but is being left for legacy purposes
	virtual void dump_BRIX(std::string const & prefix) const = 0;
	/// @brief Serialize the GridBase object into a json_spirit Value
	virtual utility::json_spirit::Value serialize() const = 0;
	/// @brief deserialize a json spirit Value into a GridBase object
	virtual void deserialize(utility::json_spirit::mObject data) = 0;
	/// @brief determine if all residue atoms are in a grid
	virtual bool is_in_grid(core::conformation::UltraLightResidue const & residue) const =0;
	/// @brief determine if all residue atoms are in a grid
	virtual bool is_in_grid(core::conformation::Residue const & residue) const =0;

	/// @brief Print a brief summary about this grid to the provided output stream
	virtual void show( std::ostream & out ) const = 0;
};

}
}
}

#endif /* GRIDBASE_HH_ */
