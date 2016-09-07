// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/protein_interface_design/DockDesignFilterFactory.hh
/// @brief
/// @author ashworth

#ifndef INCLUDED_protocols_protein_interface_design_DockDesignFilterFactory_hh
#define INCLUDED_protocols_protein_interface_design_DockDesignFilterFactory_hh

// Unit Headers
#include <protocols/protein_interface_design/DockDesignFilterFactory.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

// c++ headers
#include <map>

#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {

class DockDesignFilterFactory : public utility::pointer::ReferenceCount
{
public:
	typedef utility::tag::Tag Tag;
	typedef utility::tag::TagCOP TagCOP;
	typedef core::pose::Pose Pose;

public:
	DockDesignFilterFactory();
	~DockDesignFilterFactory() override;

	/// @brief add a Filter prototype, using its default type name as the map key
	void add_type( protocols::filters::FilterOP );
	/// @brief add a Filter prototype, using an arbitrary type name as the map key
	void add_type( std::string const &, protocols::filters::FilterOP );
	/// @brief return new Filter by key lookup in dock_design_filter_map_
	protocols::filters::FilterOP newFilter( std::string const & );
	/// @brief return new Filter by Tag parsing
	protocols::filters::FilterOP newFilter(
		TagCOP,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		moves::Movers_map const &,
		Pose const & );

private:
	protocols::filters::Filters_map dock_design_filter_map_;

};

} //namespace protein_interface_design
} //namespace protocols

#endif
