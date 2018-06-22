// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/InterfaceHydrophobicResidueContactsFilter.hh
/// @brief Counts the number of hydrophobic residues on the target that have at least a certain hydrophobic score
/// @author Brian Coventry (bcov@uw.edu)

#ifndef INCLUDED_protocols_simple_filters_InterfaceHydrophobicResidueContactsFilter_hh
#define INCLUDED_protocols_simple_filters_InterfaceHydrophobicResidueContactsFilter_hh

#include <protocols/simple_filters/InterfaceHydrophobicResidueContactsFilter.fwd.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/scoring/ScoreType.hh>


#include <set>

namespace protocols {
namespace simple_filters {

class InterfaceHydrophobicResidueContactsFilter : public filters::Filter
{
public:
	// Finally the defaults are defined in exactly once
	static const core::Size DEFAULT_THRESHOLD;
	static const std::string DEFAULT_APOLAR_RES;
	static const core::Real DEFAULT_SCORE_CUT;
	static const std::set<core::scoring::ScoreType> SCORE_TERMS;

	// InterfaceHydrophobicResidueContactsFilter() : Filter( "InterfaceHydrophobicResidueContacts" ) {}
	InterfaceHydrophobicResidueContactsFilter(
		core::Size const hydrophobic_residue_contacts_threshold = DEFAULT_THRESHOLD,
		core::select::residue_selector::ResidueSelectorCOP target_selector = nullptr,
		core::select::residue_selector::ResidueSelectorCOP binder_selector = nullptr,
		core::scoring::ScoreFunctionCOP scorefxn = nullptr,
		core::Real score_cut = DEFAULT_SCORE_CUT,
		std::string const & apolar_res = DEFAULT_APOLAR_RES

	) : Filter( "InterfaceHydrophobicResidueContacts" ),
		hydrophobic_residue_contacts_threshold_( hydrophobic_residue_contacts_threshold ),
		target_selector_( target_selector ),
		binder_selector_( binder_selector ),
		score_cut_( score_cut )
	{
		set_scorefxn( scorefxn );
		set_apolar_res( apolar_res );
	}



	bool apply( core::pose::Pose const & pose ) const override;
	filters::FilterOP clone() const override {
		return filters::FilterOP( new InterfaceHydrophobicResidueContactsFilter( *this ) );
	}
	filters::FilterOP fresh_instance() const override{
		return filters::FilterOP( new InterfaceHydrophobicResidueContactsFilter() );
	}
	InterfaceHydrophobicResidueContactsFilter &
	operator=( InterfaceHydrophobicResidueContactsFilter const & ot ) {
		if ( this != & ot ) {
			hydrophobic_residue_contacts_threshold_ = ot.hydrophobic_residue_contacts_threshold_;
			target_selector_ = ot.target_selector_;
			binder_selector_ = ot.binder_selector_;
			scorefxn_ = ot.scorefxn_;
			score_cut_ = ot.score_cut_;
			apolar_res_set_ = ot.apolar_res_set_;
		}
		return *this;
	}
	InterfaceHydrophobicResidueContactsFilter( InterfaceHydrophobicResidueContactsFilter const & init ) :
		Filter( init ) {
		*this = init; // the default class operator = should copy everything sufficiently
	}
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	core::Size compute( core::pose::Pose const & pose ) const;
	~InterfaceHydrophobicResidueContactsFilter() override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	void set_scorefxn( core::scoring::ScoreFunctionCOP scorefxn );
	void set_apolar_res( std::string const & apolar_res );
	void sanity_check() const;


private:
	core::Size hydrophobic_residue_contacts_threshold_;
	core::select::residue_selector::ResidueSelectorCOP target_selector_;
	core::select::residue_selector::ResidueSelectorCOP binder_selector_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::Real score_cut_;
	std::set<std::string> apolar_res_set_;
};

}
}

#endif
