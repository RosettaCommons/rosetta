// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/ResfileReader.hh
/// @brief  header of classes for resfile options
/// @author Gordon Lemmon

#ifndef INCLUDED_protocols_ligand_docking_FinalMinimizer_hh
#define INCLUDED_protocols_ligand_docking_FinalMinimizer_hh

// Unit Headers
#include <protocols/ligand_docking/FinalMinimizer.fwd.hh>
#include <protocols/ligand_docking/MoveMapBuilder.fwd.hh>

// Package Headers
//// Project Headers
#include <protocols/moves/Mover.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/minimization_packing/MinMover.fwd.hh>

//// Scripter Headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

// Utility Headers


///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace ligand_docking {

/// @brief
class FinalMinimizer : public protocols::moves::Mover
{
public:
	FinalMinimizer();
	FinalMinimizer(core::scoring::ScoreFunctionOP score_fxn, MoveMapBuilderOP movemap_builder, bool remove_constraints=false);
	~FinalMinimizer() override;
	FinalMinimizer(FinalMinimizer const & that);

	void apply( core::pose::Pose & pose ) override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &
	) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

public: //Citation manager functions

	/// @brief Provide the citation.
	void provide_citation_info(basic::citation_manager::CitationCollectionList & ) const override;

private:
	core::scoring::ScoreFunctionOP score_fxn_;
	MoveMapBuilderOP movemap_builder_;
	/// @brief If true, remove any backbone constraints added during minimization
	bool remove_bb_constraints_ = false;
	/// @brief Use Cartesian minimization
	bool cartesian_ = false;

	/// @brief Get the minimization submover
	/// If backbone is true, make sure that things are set up for backbone minimization
	protocols::minimization_packing::MinMoverOP const
	get_final_min_mover(core::pose::Pose const & pose, bool backbone = true) const;
};

} //namespace ligand_docking
} //namespace protocols

#endif
