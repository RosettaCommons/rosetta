// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/ligand_docking/scoring_grid/AtrGrid.hh
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_qsar_scoring_grid_AtrGrid_hh
#define INCLUDED_protocols_qsar_scoring_grid_AtrGrid_hh

#include <protocols/qsar/scoring_grid/AtrGrid.fwd.hh>
#include <protocols/qsar/scoring_grid/SingleGrid.hh>

#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace qsar {
namespace scoring_grid {

class AtrGrid : public SingleGrid
{
public:

	AtrGrid();
	virtual void refresh(core::pose::Pose const & pose, core::Vector const & center, core::Size const & ligand_chain_id_to_exclude);
	virtual void refresh(core::pose::Pose const & pose, core::Vector const & center);
	virtual void refresh(core::pose::Pose const & pose, core::Vector const & center, utility::vector1<core::Size> ligand_chain_ids_to_exclude);
	/// @brief serialize the SingleGrid to a json_spirit object
	virtual utility::json_spirit::Value serialize();
	/// @brief deserialize a json_spirit object to a SingleGrid
	virtual void deserialize(utility::json_spirit::mObject data);
	void parse_my_tag(utility::tag::TagCOP tag);

	static std::string grid_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	void set_protein_rings( core::conformation::Residue const & rsd);
	void set_ligand_rings(
		core::conformation::Residue const & rsd,
		utility::vector1<core::Size> ligand_chain_ids_to_exclude
	);

	core::Real inner_radius_;
	core::Real outer_radius_;
	core::Real bb_; // score for clashes with a backbone atom
	core::Real sc_; // score for clashes with a side-chain atom
	core::Real ligand_; // score for clashes with a non-protein atom
};

}
}
}

#endif /* ATRGRID_CC_ */
