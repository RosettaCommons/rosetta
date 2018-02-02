// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/fold_from_loops/RmsdFromResidueSelectorFilter.hh
/// @brief  Evaluate RMSD between two poses allowing to select the regions to compare in each pose through ResidueSelector
/// @author Jaume Bonet (jaume.bonet@gmail.com)

#ifndef INCLUDED_protocols_fold_from_loops_filters_RmsdFromResidueSelectorFilter_hh
#define INCLUDED_protocols_fold_from_loops_filters_RmsdFromResidueSelectorFilter_hh

#include <protocols/fold_from_loops/filters/RmsdFromResidueSelectorFilter.fwd.hh>

#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <list>

#include <utility/vector1.hh>


namespace protocols {
namespace fold_from_loops {
namespace filters {

class RmsdFromResidueSelectorFilter : public protocols::filters::Filter
{
public:
	RmsdFromResidueSelectorFilter();
	~RmsdFromResidueSelectorFilter() override;

	inline protocols::filters::FilterOP clone() const override {
		return protocols::filters::FilterOP( new RmsdFromResidueSelectorFilter( *this ) );
	};
	inline protocols::filters::FilterOP fresh_instance() const override {
		return protocols::filters::FilterOP( new RmsdFromResidueSelectorFilter() );
	};
	void reference_pose( core::pose::PoseOP const & ref ) {
		reference_pose_->detached_copy( *ref );
		reference_pose_->transfer_constraint_set( *ref );
	}
	void reference_pose( core::pose::Pose const & ref )  {
		reference_pose_->detached_copy( ref );
		reference_pose_->transfer_constraint_set( ref );
	}
	void threshold( core::Real threshold ) { threshold_ = threshold; }
	void CA_only( bool pick ) { CA_only_ = pick; }
	void GDT( bool pick ) { gdt_ = pick; }
	void superimpose( bool pick ) { superimpose_ = pick; }
	void count_residues( bool pick ) { count_residues_ = pick; }
	void reference_selector( core::select::residue_selector::ResidueSelectorCOP const &  select ) { reference_select_ = select; }
	void reference_selector( core::select::residue_selector::ResidueSelector const &  select ) { reference_select_ = select.clone(); }
	void query_selector( core::select::residue_selector::ResidueSelectorCOP const &  select ) { query_select_ = select; }
	void query_selector( core::select::residue_selector::ResidueSelector const &  select ) { query_select_ = select.clone(); }
	void selectors(
		core::select::residue_selector::ResidueSelectorCOP const &  reference,
		core::select::residue_selector::ResidueSelectorCOP const &  query) {
		reference_selector(reference);
		query_selector(query);
	};
	void selectors(
		core::select::residue_selector::ResidueSelector const &  reference,
		core::select::residue_selector::ResidueSelector const &  query) {
		reference_selector(reference);
		query_selector(query);
	};
	bool apply( core::pose::Pose const & pose ) const override;
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	core::Real compute( core::pose::Pose const & pose ) const;

	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data_map, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & reference_pose ) override;

	std::string
	name() const override {
		return class_name();
	};

	static
	std::string
	class_name() {
		return "RmsdFromResidueSelectorFilter";
	};

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	inline static core::Real default_rmsd_threshold() { return 5.0; }
	inline static bool default_ca_selection() { return true; }
	inline static bool default_gdt_selection() { return false; }
	inline static bool default_superimpose_selection() { return true; }
	inline static bool default_count_residues() { return true; }

private:
	core::select::residue_selector::ResidueSelectorCOP reference_select_;
	core::select::residue_selector::ResidueSelectorCOP query_select_;
	core::Real threshold_;
	core::pose::PoseOP reference_pose_;
	core::Real sensitivity_;
	bool CA_only_;
	bool gdt_;
	bool superimpose_;
	bool count_residues_;

};

}
} // fold_from_loops
} // protocols

#endif //INCLUDED_protocols_protein_interface_design_filters_RmsdFilter_HH_
