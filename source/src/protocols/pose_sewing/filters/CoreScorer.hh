// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/filters/CoreScorer.hh
/// @brief a filter that evaluates pairwise MotifScores
/// @author frankdt (frankdt@email.unc.edu)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_pose_sewing_filters_CoreScorer_hh
#define INCLUDED_protocols_pose_sewing_filters_CoreScorer_hh

// Unit headers
#include <protocols/pose_sewing/filters/CoreScorer.fwd.hh>
#include <protocols/filters/Filter.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <map>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
//#include <utility/tag/Tag.fwd.hh> //transcluded from Filter.hh
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Filter.hh

namespace protocols {
namespace pose_sewing {
namespace filters {

struct CSSettings {

	core::Real score_cutoff_;
	core::Size link_cutoff_;
	core::Size window_width_;

};

///@brief a filter that evaluates pairwise MotifScores
class CoreScorer : public protocols::filters::Filter {

public:
	CoreScorer();

	// destructor (important for properly forward-declaring smart-pointer members)
	~CoreScorer() override;

	/// @brief returns true if the structure passes the filter, false otherwise
	bool
	apply( core::pose::Pose const & pose ) const override;

	/// @brief required for reporting score values
	core::Real
	report_sm( core::pose::Pose const & pose ) const override;

	/// @brief allows printing data to a stream
	void
	report( std::ostream & os, core::pose::Pose const & pose ) const override;

	//void
	//populate_secstruct_shift_array(core::pose::Pose const & pose);

	//core::Size
	//get_secstruct_shift(core::Size N_resnum) const;

public:
	std::string
	name() const override;

	static
	std::string
	class_name();

	void
	set_selector(core::select::residue_selector::ResidueSelectorCOP selector);

	void
	set_score_cutoff(core::Real score_cutoff);

	void
	set_distance_cutoff(core::Real distance_cutoff);

	void
	set_window_width(core::Size window_width);

	void
	set_link_cutoff(core::Size link_cutoff);

	void
	set_sum(bool sum);

	void
	set_distance_mode(bool distance_mode);

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief parse XML tag (to use this Filter in Rosetta Scripts)
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::filters::FilterOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::filters::FilterOP
	clone() const override;

private:
	core::Real score_cutoff_;
	core::Real distance_cutoff_;
	core::Size link_cutoff_ = 3;
	core::Size window_width_;
	core::select::residue_selector::ResidueSelectorCOP selector_;
	std::map < char, CSSettings > ss_settings_;
	bool sum_ = false;
	bool distance_mode_ = false;
};

} //protocols
} //pose_sewing
} //filters

#endif //INCLUDED_protocols_pose_sewing_filters_CoreScorer_hh
