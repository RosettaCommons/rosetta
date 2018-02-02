// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/fold_from_loops/ScorePoseSegmentFromResidueSelectorFilter.hh
/// @brief  Applies the scoring function only on the selected segment of the pose.
/// @author Jaume Bonet (jaume.bonet@gmail.com)

#ifndef INCLUDED_protocols_fold_from_loops_filters_ScorePoseSegmentFromResidueSelectorFilter_hh
#define INCLUDED_protocols_fold_from_loops_filters_ScorePoseSegmentFromResidueSelectorFilter_hh

#include <protocols/fold_from_loops/filters/ScorePoseSegmentFromResidueSelectorFilter.fwd.hh>

#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <list>

#include <utility/vector1.hh>


namespace protocols {
namespace fold_from_loops {
namespace filters {

class ScorePoseSegmentFromResidueSelectorFilter : public protocols::filters::Filter
{
public:
	ScorePoseSegmentFromResidueSelectorFilter();
	~ScorePoseSegmentFromResidueSelectorFilter() override;

	inline protocols::filters::FilterOP clone() const override {
		return protocols::filters::FilterOP( new ScorePoseSegmentFromResidueSelectorFilter( *this ) );
	};
	inline protocols::filters::FilterOP fresh_instance() const override {
		return protocols::filters::FilterOP( new ScorePoseSegmentFromResidueSelectorFilter() );
	};

	void residue_selector( core::select::residue_selector::ResidueSelectorCOP const &  select ) { residue_select_ = select; }
	void residue_selector( core::select::residue_selector::ResidueSelector const &  select ) { residue_select_ = select.clone(); }
	core::scoring::ScoreFunctionOP scorefxn() const { return scorefxn_; };
	void scorefxn( core::scoring::ScoreFunctionOP const & scorefxn ) { scorefxn_ = scorefxn; };
	std::string score_id_name() const { return score_id_name_; };
	void score_id_name( std::string name ) { score_id_name_ = name; };
	bool in_context() const { return in_context_; };
	void in_context( bool pick ) { in_context_ = pick; };

	bool apply( core::pose::Pose const & pose ) const override;
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	core::Real compute( core::pose::Pose const & pose ) const;

	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data_map, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;

	std::string
	name() const override {
		return class_name();
	};

	static
	std::string
	class_name() {
		return "ScorePoseSegmentFromResidueSelectorFilter";
	};

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	inline static core::select::residue_selector::ResidueSelectorCOP default_selector() {
		return core::select::residue_selector::ResidueSelectorCOP( new core::select::residue_selector::TrueResidueSelector );
	}
	inline static core::scoring::ScoreFunctionOP default_scorefxn() { return core::scoring::get_score_function(); }

private:
	std::string score_id_name_;
	core::select::residue_selector::ResidueSelectorCOP residue_select_;
	core::scoring::ScoreFunctionOP scorefxn_;
	bool in_context_;

};

}
} // fold_from_loops
} // protocols

#endif //INCLUDED_protocols_protein_interface_design_filters_ScorePoseSegmentFilter_HH_
