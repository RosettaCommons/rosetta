// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/operation/ResLvlTaskOperationFactory.hh
/// @brief
/// @author ashworth

#ifndef INCLUDED_core_pack_task_operation_ResLvlTaskOperationFactory_hh
#define INCLUDED_core_pack_task_operation_ResLvlTaskOperationFactory_hh

// Unit Headers
#include <core/pack/task/operation/ResLvlTaskOperationFactory.fwd.hh>

// Package Headers
#include <core/pack/task/operation/ResLvlTaskOperation.fwd.hh>
#include <core/pack/task/operation/ResLvlTaskOperationCreator.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

// c++ headers
// AUTO-REMOVED #include <string>
#include <map>

#ifdef PYROSETTA
	#include <utility/tag/Tag.hh>
#endif


namespace core {
namespace pack {
namespace task {
namespace operation {

// singleton class
class ResLvlTaskOperationFactory
{
public:
	typedef std::map< std::string, ResLvlTaskOperationCreatorOP > RLTOC_Map;
	typedef utility::tag::Tag Tag;
	typedef utility::tag::TagCOP TagCOP;

public:
	static ResLvlTaskOperationFactory * get_instance();
	void factory_register( ResLvlTaskOperationCreatorOP );

	///@brief add a prototype, using its default type name as the map key
	void add_creator( ResLvlTaskOperationCreatorOP );
	bool has_type( std::string const & ) const;

///@brief return new ResLvlTaskOperation by key lookup in rlto_map_ (new ResLvlTaskOperation parses Tag if provided)
	ResLvlTaskOperationOP newRLTO( std::string const & ) const;

private:
	ResLvlTaskOperationFactory();
	virtual ~ResLvlTaskOperationFactory();

	static ResLvlTaskOperationFactory * instance_;
	RLTOC_Map rltoc_map_;

};

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core

#endif
