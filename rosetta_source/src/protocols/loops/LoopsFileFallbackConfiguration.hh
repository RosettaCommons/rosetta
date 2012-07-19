// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/looops/LoopsFileFallbackConfiguration.hh
/// @author Brian D. Weitzner brian.weitzner@gmail.com

#ifndef INCLUDED_protocols_loops_loops_file_fallback_configuration_HH
#define INCLUDED_protocols_loops_loops_file_fallback_configuration_HH

// Unit Header
#include <protocols/loops/LoopsFileFallbackConfiguration.fwd.hh>
#include <basic/resource_manager/FallbackConfiguration.hh>
#include <basic/resource_manager/types.hh>

namespace protocols {
namespace loops {

class LoopsFileFallbackConfiguration : public basic::resource_manager::FallbackConfiguration {

public:
	LoopsFileFallbackConfiguration();
	
	virtual
	basic::resource_manager::ResourceTag const & get_resource_tag_from_description(
		basic::resource_manager::ResourceDescription const & desc
	) const;
	
	virtual
	basic::resource_manager::LocatorTag const & get_locator_tag_from_description(
		basic::resource_manager::ResourceDescription const & desc
	) const;
	
	virtual
	basic::resource_manager::LocatorID const & get_locator_id_from_description(
		basic::resource_manager::ResourceDescription const & desc
	) const;
	
	virtual
	basic::resource_manager::LoaderType const & get_loader_type_from_description(
		basic::resource_manager::ResourceDescription const & desc
	) const;
	
	virtual
	basic::resource_manager::ResourceOptionsTag const & get_resource_options_tag_from_description(
		basic::resource_manager::ResourceDescription const & desc
	) const;

private:
	static basic::resource_manager::LoaderType loader_type_;
	static basic::resource_manager::LocatorID locator_id_;
	static basic::resource_manager::LocatorTag locator_tag_;
	static basic::resource_manager::ResourceTag resource_tag_;
	static basic::resource_manager::ResourceOptionsTag resource_options_tag_;
	
private:
	basic::resource_manager::LocatorID get_loops_filename_from_options() const;
	
};

} // namespace loops
} // namespace protocols

#endif // INCLUDED_protocols_loops_loops_file_fallback_configuration_HH
