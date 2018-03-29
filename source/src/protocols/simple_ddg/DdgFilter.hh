// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_ddg/DdgFilter.hh
/// @brief definition of filter class DdgFilter.
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)
/// Augmented for PB support by Sachko Honda (honda@apl.washington.edu) December 2012.

#ifndef INCLUDED_protocols_simple_filters_DdgFilter_hh
#define INCLUDED_protocols_simple_filters_DdgFilter_hh

#include <protocols/simple_ddg/DdgFilter.fwd.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/TaskFactory.hh>

namespace protocols {
namespace simple_ddg {

class DdgFilter : public filters::Filter
{
public:
	DdgFilter();
	DdgFilter( core::Real const ddg_threshold,
		core::scoring::ScoreFunctionCOP scorefxn,
		core::Size const rb_jump=1,
		core::Size const repeats=1);
	bool apply( core::pose::Pose const & pose ) const override;
	filters::FilterOP clone() const override;
	filters::FilterOP fresh_instance() const override;

	core::Size repeats() const;
	void repeats( core::Size const repeats );
	void repack( bool const repack );
	bool repack() const;
	void repack_bound( bool rpb ) { repack_bound_ = rpb; }
	bool repack_bound() const { return repack_bound_; }
	void relax_bound( bool rlb ) { relax_bound_ = rlb; }
	bool relax_bound() const { return relax_bound_; }
	void translate_by( core::Real const translate_by );
	core::Real translate_by() const;
	void task_factory( core::pack::task::TaskFactoryOP task_factory ) { task_factory_ = task_factory; }
	core::pack::task::TaskFactoryOP task_factory() const { return task_factory_; }
	void use_custom_task( bool uct ) { use_custom_task_ = uct; }
	bool use_custom_task() const { return use_custom_task_; }
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	core::Real compute( core::pose::Pose const & pose ) const;
	~DdgFilter() override;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;
	void relax_mover( protocols::moves::MoverOP m );
	protocols::moves::MoverOP relax_mover() const;
	void filter( protocols::filters::FilterOP m );
	protocols::filters::FilterOP filter() const;

	void extreme_value_removal( bool const b ) { extreme_value_removal_ = b; }
	bool extreme_value_removal() const{ return extreme_value_removal_; }

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:

	// initialize PB related features
	void initPB();

	core::Real ddg_threshold_; //dflt -15
	core::scoring::ScoreFunctionOP scorefxn_; //dflt NULL
	core::Size rb_jump_; // dflt 1
	core::pack::task::TaskFactoryOP task_factory_;
	bool use_custom_task_;
	bool repack_bound_; //dflt true; Do you want to repack in the bound state (ddG). Avoid redundant packing if already packed beforing calling the filter.
	bool relax_bound_; //dflt false; Do you want to relax in the bound state (ddG). Avoid redundant relax if already relaxed before calling the filter.
	utility::vector1<core::Size> chain_ids_;
	core::Size repeats_;//average of how many repeats? defaults to 1
	bool repack_; //dflt true; Do you want to repack in the bound and unbound states (ddG) or merely compute the dG
	protocols::moves::MoverOP relax_mover_; //dflt NULL; in the unbound state, prior to taking the energy, should we do any relaxation
	protocols::filters::FilterOP filter_; //dflt NULL; use a filter instead of the scorefunction

	/// is PB enabled?
	bool pb_enabled_;

	/// translation distance in A
	core::Real translate_by_; //dflt set to 1000. Default resets to 100 for RosettaScripts with a scorefxn containing PB_elec.
	bool extreme_value_removal_; //dflt false; if true, and repeats >= 3 removes the lowest and highest value computed during repeats. This should lower noise levels
};


}
}
#endif
