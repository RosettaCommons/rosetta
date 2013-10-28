// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/operation/ResFilterFactory.hh
/// @brief
/// @author ashworth

#ifndef INCLUDED_core_pack_task_operation_ResFilterFactory_hh
#define INCLUDED_core_pack_task_operation_ResFilterFactory_hh

// Unit Headers
#include <core/pack/task/operation/ResFilterFactory.fwd.hh>

// Package Headers
#include <core/pack/task/operation/ResFilter.fwd.hh>
#include <core/pack/task/operation/ResFilterCreator.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.hh>

// c++ headers
#include <string>
#include <map>

#include <utility/vector0.hh>


namespace core {
namespace pack {
namespace task {
namespace operation {

class ResFilterFactory // singleton; no need to derive from RefCount : public utility::pointer::ReferenceCount
{
public:
	//typedef utility::pointer::ReferenceCount parent;
	typedef std::map< std::string, ResFilterCreatorOP > ResFilterCreatorMap;
	typedef utility::tag::Tag Tag;
	typedef utility::tag::TagCOP TagCOP;

public:
	static ResFilterFactory * get_instance();
	void factory_register( ResFilterCreatorOP );

	///@brief add a prototype, using its default type name as the map key
	void add_creator( ResFilterCreatorOP );
	bool has_type( std::string const & ) const;

///@brief return new ResFilter by key lookup in filter_map_ (new ResFilter parses Tag if provided)
	ResFilterOP newResFilter( std::string const &, TagCOP = new Tag ) const;

private:
	ResFilterFactory();
	virtual ~ResFilterFactory();

	static ResFilterFactory * instance_;
	ResFilterCreatorMap filter_creator_map_;

};

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core

#endif
