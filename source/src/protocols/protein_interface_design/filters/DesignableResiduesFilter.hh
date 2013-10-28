// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/DesignableResiduesFilter.hh
/// @brief Reports to Tracer which residues are designable in a taskfactory
/// @author Sarel Fleishman (sarelf@uw.edu)

#ifndef INCLUDED_protocols_protein_interface_design_filters_DesignableResiduesFilter_hh
#define INCLUDED_protocols_protein_interface_design_filters_DesignableResiduesFilter_hh


// Project Headers
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/protein_interface_design/filters/DesignableResiduesFilter.fwd.hh>

#include <core/pack/task/TaskFactory.fwd.hh>
#include <utility/vector1.hh>

// Unit headers

namespace protocols {
namespace protein_interface_design{
namespace filters {

class DesignableResiduesFilter : public protocols::filters::Filter

{
private:
	typedef protocols::filters::Filter parent;
public:
	/// @brief default ctor
	DesignableResiduesFilter();
	///@brief Constructor with a single target residue
	virtual bool apply( core::pose::Pose const & pose ) const;
	virtual void report( std::ostream & out, core::pose::Pose const & pose ) const;
	virtual core::Real report_sm( core::pose::Pose const & pose ) const;
	virtual protocols::filters::FilterOP clone() const;
	virtual protocols::filters::FilterOP fresh_instance() const;
	core::Size compute( core::pose::Pose const & pose ) const;
	virtual ~DesignableResiduesFilter();
	core::pack::task::TaskFactoryOP task_factory() const;
	void task_factory( core::pack::task::TaskFactoryOP task_factory );
	void parse_my_tag( utility::tag::TagCOP const tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );
	void parse_def( utility::lua::LuaObject const & def,
					utility::lua::LuaObject const & score_fxns,
					utility::lua::LuaObject const & tasks );
	core::Size lower_threshold() const;
	core::Size upper_threshold() const;
	bool packable() const;
	bool designable() const;
	void lower_threshold( core::Size const l );
	void upper_threshold( core::Size const u );
	void packable( bool const p );
	void designable( bool const d );
private:
	core::pack::task::TaskFactoryOP task_factory_;
	core::Size lower_threshold_, upper_threshold_; // how many design positions should be allowed for a passing design  
	bool packable_, designable_; // which sort of residues to report (packable or designable or both)
};

} // filters
} //protein_interface_design
} // protocols

#endif //INCLUDED_protocols_Filters_DesignableResiduesFilter_HH_

