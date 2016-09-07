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
#include <protocols/jd2/Job.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/qsar/scoring_grid/ScoreNormalization.hh>


#include <utility/vector1.hh>


///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace ligand_docking {

class InterfaceScoreCalculator : public protocols::moves::Mover
{
public:
	InterfaceScoreCalculator();
	~InterfaceScoreCalculator() override;
	InterfaceScoreCalculator(InterfaceScoreCalculator const & that);

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;
	std::string get_name() const override;

	void chains(std::vector<std::string> const & chains);

	void score_fxn(core::scoring::ScoreFunctionOP const & score_fxn);

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	void apply(core::pose::Pose & pose) override;

private:
	utility::vector1<std::string> chains_;
	core::pose::PoseOP native_;
	core::scoring::ScoreFunctionOP score_fxn_;
	protocols::qsar::scoring_grid::ScoreNormalizationOP normalization_function_;
	bool compute_grid_scores_;
	std::string prefix_;

	void
	add_scores_to_job(
		core::pose::Pose & pose,
		protocols::jd2::JobOP job
	) const;

	void
	append_ligand_docking_scores(
		core::pose::Pose const & after,
		protocols::jd2::JobOP job
	) const;

	void
	append_ligand_docking_scores(
		core::Size jump_id,
		core::pose::Pose const & after,
		protocols::jd2::JobOP job
	) const;

};

} //namespace ligand_docking
} //namespace protocols

#endif
