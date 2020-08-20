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

	void chains(std::vector<std::string> const & chains);

	void set_native(core::pose::PoseOP setting ) { native_ = setting; }

	void native_ensemble_best( bool setting ) { native_ensemble_best_ = setting; }

	void score_fxn(core::scoring::ScoreFunctionOP const & score_fxn);

	void grid_set_prototype(protocols::qsar::scoring_grid::GridSetCOP grid_prototype);

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &
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

public: //Citation info:

	///@brief Does this mover provide information about how to cite it?
	///@details Returns true.
	bool mover_provides_citation_info() const override;

	/// @brief Provide the citation.
	/// @returns A vector of citation collections.  This allows the mover to provide citations for itself and for any modules that it invokes.
	/// @details Also provides citations for movers called by the InterfaceScoreCalculator.
	utility::vector1< basic::citation_manager::CitationCollectionCOP > provide_citation_info() const override;

	/// @brief Provide a list of authors and their e-mail addresses, as strings.
	/// @returns A list of pairs of (author, e-mail address).  This mover IS published, so it returns nothing for itself,
	/// but can return  information for preselection filters and movers.
	utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP > provide_authorship_info_for_unpublished() const override;

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
	bool native_ensemble_best_;
	core::scoring::ScoreFunctionOP score_fxn_;
	protocols::qsar::scoring_grid::ScoreNormalizationOP normalization_function_;
	bool compute_grid_scores_;
	protocols::qsar::scoring_grid::GridSetCOP grid_set_prototype_;
	std::string prefix_;
	//Compute mem_interface_delta_X by translating the ligand in the membrane.
	bool score_in_mem_;

};

} //namespace ligand_docking
} //namespace protocols

#endif
