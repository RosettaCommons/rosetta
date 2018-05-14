// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/BuriedUnsatHbondFilter.hh
/// @brief filter to report buried polar atoms that are not satisfied by hydrogen bonds
/// @details significantly updated in 2017 to have more generous definition of h-bonds
///  (previously, legit h-bonds were excluded because of sfxn exceptions); users can now choose
///  between legacy SASA and VSASA for burial; poses with more than 3 chains now supported; the
///  way unsats are counted and reported is now different (before, all unsats were counted as equal,
///  which is misleading); users can choose different reporting behaviours; the Filter is now Symmetry
///  compatible, and users can pass sym_dof_names as in Jacob Bale's earlier modified filter; users can
///  still use ddG-style behavior if desired; legacy options can be restored by setting legacy=true,
///  but this is only recommended for benchmarking purposes
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)
/// @author Scott Boyken (sboyken@gmail.com), major updates and refactoring

#ifndef INCLUDED_protocols_simple_filters_BuriedUnsatHbondFilter_hh
#define INCLUDED_protocols_simple_filters_BuriedUnsatHbondFilter_hh

#include <protocols/simple_filters/BuriedUnsatHbondFilter.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.hh>

namespace protocols {
namespace simple_filters {

/// @brief filters based on an upper bound # of buried unsatisfied polar residues
class BuriedUnsatHbondFilter : public filters::Filter
{
public:
	// @brief default constructor
	BuriedUnsatHbondFilter();
	// @brief constructor with options
	BuriedUnsatHbondFilter(core::Size const upper_threshold,
		core::select::residue_selector::ResidueSelectorOP residue_selector=NULL,
		bool const use_legacy_options=false,
		bool const generous_hbonds=true,
		bool const use_vasa=true,
		bool const ignore_surface_res=false,
		bool const ignore_bb_heavy_unsats=false,
		bool const use_sc_neighbors=false,
		bool const only_interface=false,
		bool const use_ddG_style=false,
		bool const use_hbnet_behavior=false
	);
	// @brief copy constructor
	BuriedUnsatHbondFilter( BuriedUnsatHbondFilter const & rval );
	// @brief make clone
	filters::FilterOP clone() const override { return filters::FilterOP( new BuriedUnsatHbondFilter( *this ) ); }
	// @brief make fresh instance
	filters::FilterOP fresh_instance() const override { return filters::FilterOP( new BuriedUnsatHbondFilter() ); }
	bool apply( core::pose::Pose const & pose ) const override;
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	core::Real compute( core::pose::Pose const & pose ) const;

	~BuriedUnsatHbondFilter() override;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;
	void task_factory( core::pack::task::TaskFactoryOP tf );
	core::pack::task::TaskFactoryOP task_factory() const;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	// getters and setters
	void set_generous_hbonds_( bool generous_hbonds ){ generous_hbonds_ = generous_hbonds; };
	void set_use_vsasa( bool use_vsasa ){ use_vsasa_ = use_vsasa; };
	void set_ignore_surface_res( bool ignore_surface_res ){ ignore_surface_res_ = ignore_surface_res; };
	void set_ignore_bb_heavy_unsats( bool ignore_bb_heavy_unsats ){ ignore_bb_heavy_unsats_ = ignore_bb_heavy_unsats; };
	void set_use_ddG_style( bool use_ddG_style ){ use_ddG_style_ = use_ddG_style; };
	void set_ddG_style_dont_recalc_surface( bool ddG_style_dont_recalc_surface ) { ddG_style_dont_recalc_surface_ = ddG_style_dont_recalc_surface; }
	void set_print_out_info_to_pdb( bool print_out_info_to_pdb ){ print_out_info_to_pdb_ = print_out_info_to_pdb; };
	void set_use_reporter_behavior( bool use_reporter_behavior ){ use_reporter_behavior_ = use_reporter_behavior; };
	void set_use_hbnet_behavior( bool use_hbnet_behavior ){ use_hbnet_behavior_ = use_hbnet_behavior; };
	void set_report_all_unsats( bool report_all_unsats ){ report_all_unsats_ = report_all_unsats; };
	void set_report_all_heavy_atom_unsats( bool report_all_heavy_atom_unsats ){ report_all_heavy_atom_unsats_ = report_all_heavy_atom_unsats; };
	void set_report_sc_heavy_atom_unsats( bool report_sc_heavy_atom_unsats ){ report_sc_heavy_atom_unsats_ = report_sc_heavy_atom_unsats; };
	void set_report_bb_heavy_atom_unsats( bool report_bb_heavy_atom_unsats ){ report_bb_heavy_atom_unsats_ = report_bb_heavy_atom_unsats; };
	void set_report_nonheavy_unsats( bool report_nonheavy_unsats ){ report_nonheavy_unsats_ = report_nonheavy_unsats; };

	void set_probe_radius( core::Real probe_radius ){ probe_radius_ = probe_radius; };
	void set_burial_cutoff( core::Real burial_cutoff ){ burial_cutoff_ = burial_cutoff; };
	void set_residue_surface_cutoff( core::Real residue_surface_cutoff ){ residue_surface_cutoff_ = residue_surface_cutoff; };
	void set_upper_threshold( core::Size upper_threshold ){ upper_threshold_ = upper_threshold; };
	void set_jump_num( core::Size jump_num ){ jump_num_ = jump_num; };

	void sym_dof_names( std::string const & sym_dofs ) { sym_dof_names_ = sym_dofs; };
	std::string sym_dof_names() const { return sym_dof_names_; };

	void set_residue_selector( core::select::residue_selector::ResidueSelectorOP residue_selector )
	{
		residue_selector_ = residue_selector;
	};

	void set_use_sc_neighbors( bool use_sc_neighbors )
	{
		use_sc_neighbors_ = use_sc_neighbors;
		update_cutoffs();
	};

private:

	void set_legacy_options() {
		generous_hbonds_ = false;
		legacy_counting_ = true;
		use_vsasa_ = false;
		ignore_surface_res_ = false;
		ignore_bb_heavy_unsats_ = false;
		use_sc_neighbors_ = false;
		only_interface_ = false;
		use_ddG_style_ = true;
		use_reporter_behavior_ = false;
		use_hbnet_behavior_ = false;
		//jump_num_ = 1; // defaults to 1 anyway, and legacy behavior is that users can define jump number
		// get set to legacy defaults by Calculator based on the above options
		//probe_radius_ =
		//burial_cutoff_;
		//residue_surface_cutoff_;
	}

	void update_cutoffs()
	{
		//need check here in cases where called from code (non-XML)
		if ( use_sc_neighbors_ ) {
			// SASA uses atom cutoff; us_sc_neighbor uses residue cutoff
			if ( burial_cutoff_ < 1.0 ) burial_cutoff_ = 4.4; // need better solution but this catches cases where users forgot to update cutoff for sc_neighbor (rather than sasa) values
			if ( residue_surface_cutoff_ > 7.0 ) burial_cutoff_ = 2.0; // need better solution but this catches cases where users forgot to update cutoff for sc_neighbor (rather than sasa) values
		}
	};

private:
	std::string name_of_sasa_calc_;
	std::string sym_dof_names_;
	bool legacy_options_;
	bool generous_hbonds_;
	bool legacy_counting_;
	bool use_vsasa_;
	bool ignore_surface_res_;
	bool ignore_bb_heavy_unsats_;
	bool use_sc_neighbors_;
	bool only_interface_;
	bool use_ddG_style_;
	bool ddG_style_dont_recalc_surface_;
	bool print_out_info_to_pdb_;
	bool use_reporter_behavior_;
	bool use_hbnet_behavior_;
	bool report_all_unsats_;
	bool report_all_heavy_atom_unsats_;
	bool report_sc_heavy_atom_unsats_;
	bool report_bb_heavy_atom_unsats_;
	bool report_nonheavy_unsats_;
	core::Real probe_radius_;
	core::Real burial_cutoff_;
	core::Real residue_surface_cutoff_;
	core::Size upper_threshold_;
	core::Size jump_num_;
	core::select::residue_selector::ResidueSelectorOP residue_selector_;
	core::scoring::ScoreFunctionCOP sfxn_;
	core::pack::task::TaskFactoryOP task_factory_; // dflt NULL; only residues defined as packable by the taskoperations will be tested for burial
};

}
}
#endif
