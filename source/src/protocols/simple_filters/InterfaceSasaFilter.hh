// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/InterfaceSasaFilter.hh
/// @brief definition of filter class InterfaceSasaFilter.
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_simple_filters_InterfaceSasaFilter_hh
#define INCLUDED_protocols_simple_filters_InterfaceSasaFilter_hh

#include <utility/vector1.hh>
#include <protocols/simple_filters/InterfaceSasaFilter.fwd.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/scoring/ScoreType.hh>

namespace protocols {
namespace simple_filters {

class InterfaceSasaFilter : public filters::Filter
{
public:
	InterfaceSasaFilter();
	InterfaceSasaFilter( core::Real const lower_threshold, bool const hydrophobic=false, bool const polar=false, core::Real const upper_threshold=100000000.0, std::string const sym_dof_names="" );

	bool apply( core::pose::Pose const & pose ) const;
	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const & pose ) const; // which residue numbers are neighbors
	filters::FilterOP clone() const;
	filters::FilterOP fresh_instance() const;

	virtual ~InterfaceSasaFilter();
	void jump( core::Size const jump );
	void add_jump( core::Size const jump );
	void jumps( utility::vector1<core::Size> const jumps );

	void sym_dof_names( std::string const sym_dof_names );
	void sym_dof_names( utility::vector1<std::string> const sym_dof_names );
	void add_sym_dof_name( std::string const sym_dof_name );

	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
	void parse_def( utility::lua::LuaObject const & def,
		utility::lua::LuaObject const & score_fxns,
		utility::lua::LuaObject const & tasks );
private:
	core::Real lower_threshold_;
	core::Real upper_threshold_;
	bool hydrophobic_, polar_; /// count only hydrophobics? polars?
	utility::vector1<core::Size> jumps_; // dflt 1; across which jumps to compute sasa
	utility::vector1<std::string> sym_dof_names_; // dflt 1; sym_dof_names for jumps across which to compute sasa
};

}
}

#endif
