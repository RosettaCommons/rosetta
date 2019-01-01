// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/nmr/ParaNMRScoreMover.hh
/// @brief   Evaluate pose with paramagnetic NMR data.
///          Print table of experimental and predicted values to file.
///          Write scores and q-factors to scorefile.
/// @details last Modified: 11/15/18
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_protocols_nmr_ParaNMRScoreMover_HH
#define INCLUDED_protocols_nmr_ParaNMRScoreMover_HH

// Unit headers
#include <protocols/nmr/ParaNMRScoreMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Package headers
#include <core/scoring/nmr/pcs/PCSData.fwd.hh>
#include <core/scoring/nmr/rdc/RDCData.fwd.hh>
#include <core/scoring/nmr/pre/PREData.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

// Utility headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <string>

namespace protocols {
namespace nmr {

class ParaNMRScoreMover : public protocols::moves::Mover {

public: // Methods

	/// @brief Default constructor
	ParaNMRScoreMover();

	/// @brief Copy constructor
	ParaNMRScoreMover( ParaNMRScoreMover const & other );

	/// @brief Assignment operator
	ParaNMRScoreMover & operator=( ParaNMRScoreMover const & rhs );

	/// @brief Destructor
	~ParaNMRScoreMover() override;

	/// @brief Get the name of this mover
	std::string get_name() const override;

	static std::string mover_name();

	/// @brief Score pose with para NMR data and print scores and Q-factors to scorefile.
	void apply( core::pose::Pose & pose ) override;

	/// @brief Return a clone of the Mover object.
	protocols::moves::MoverOP clone() const override;

	/// @brief Generates a new Mover object freshly created with the default ctor.
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief Parse tags of XML script
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & //pose//
	) override;

	/// @brief Create XML schema definition for ParaNMRScoreMover
	static void provide_xml_schema(utility::tag::XMLSchemaDefinition & xsd);

	core::scoring::ScoreFunctionOP get_scorefunction() const { return sfxn_; }
	void set_scorefunction(core::scoring::ScoreFunctionOP sfxn) { sfxn_ = sfxn; }

	void write_calc_values(bool setting) { write_calc_values_ = setting; };
	void write_tensor_info(bool setting) { write_tensor_info__ = setting; };
	void set_verbose_scorefile(bool setting) { verbose_ = setting; };

private: // Methods

	/// @brief Write PCS scores and Q-factors to scorefile
	void add_pcs_scores_to_scorefile(
		core::pose::Pose & pose,
		core::scoring::nmr::pcs::PCSData & pcs_data,
		core::Real const pcs_wght
	) const;

	/// @brief Write RDC scores and Q-factors to scorefile
	void add_rdc_scores_to_scorefile(
		core::pose::Pose & pose,
		core::scoring::nmr::rdc::RDCData & rdc_data,
		core::Real const rdc_wght )
	const;

	/// @brief Write PRE scores and Q-factors to scorefile
	void add_pre_scores_to_scorefile(
		core::pose::Pose & pose,
		core::scoring::nmr::pre::PREData & pre_data,
		core::Real const pre_wght
	) const;

private: // Data

	/// @brief Scorefunction to be used by ParaNMRScoreMover.
	core::scoring::ScoreFunctionOP sfxn_;

	/// @brief Print separate scores and Q-factors for each tagging site, lanthanide and alignment medium.
	bool verbose_;

	/// @brief Save a table of experimental vs. predicted values to file.
	bool write_calc_values_;

	/// @brief Save tensor and/or spinlabel info to file.
	bool write_tensor_info__;

};

} // namespace nmr
} // namespace protocols

#endif // INCLUDED_protocols_nmr_ParaNMRScoreMover_HH
