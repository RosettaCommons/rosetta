// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/InterfaceScoreCalculator.hh
/// @brief
/// @author Gordon Lemmon

#ifndef INCLUDED_protocols_ligand_docking_InterfaceScoreCalculator_hh
#define INCLUDED_protocols_ligand_docking_InterfaceScoreCalculator_hh

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/ligand_docking/InterfaceScoreCalculator.fwd.hh>

// Scripter Headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

//// Project Headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/qsar/scoring_grid/ScoreNormalization.hh>
#include <protocols/qsar/scoring_grid/GridSet.fwd.hh>


#include <utility/vector1.hh>


///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace ligand_docking {

class InterfaceScoreCalculator : public protocols::moves::Mover
{
public:

	typedef std::map< std::string, core::Real > StringRealMap;

	InterfaceScoreCalculator();
	~InterfaceScoreCalculator() override;
	InterfaceScoreCalculator(InterfaceScoreCalculator const & that);

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;
	// XRW TEMP  std::string get_name() const override;

	void chains(std::vector<std::string> const & chains);

	void score_fxn(core::scoring::ScoreFunctionOP const & score_fxn);

	void grid_set_prototype(protocols::qsar::scoring_grid::GridSetCOP grid_prototype);

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	void apply(core::pose::Pose & pose) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:

	StringRealMap
	get_scores(
		core::pose::Pose & pose
	) const;

	StringRealMap
	get_ligand_docking_scores(
		core::pose::Pose const & after
	) const;

	StringRealMap
	get_ligand_docking_scores(
		char chain,
		core::pose::Pose const & after
	) const;

private:
	utility::vector1<std::string> chains_;
	core::pose::PoseOP native_;
	core::scoring::ScoreFunctionOP score_fxn_;
	protocols::qsar::scoring_grid::ScoreNormalizationOP normalization_function_;
	bool compute_grid_scores_;
	protocols::qsar::scoring_grid::GridSetCOP grid_set_prototype_;
	std::string prefix_;

};

} //namespace ligand_docking
} //namespace protocols

#endif
