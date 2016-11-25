// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/operation/TaskOperationFactory.hh
/// @brief
/// @author ashworth

#ifndef INCLUDED_core_pack_task_operation_TaskOperationFactory_hh
#define INCLUDED_core_pack_task_operation_TaskOperationFactory_hh

// Unit Headers
#include <core/pack/task/operation/TaskOperationFactory.fwd.hh>

// Package Headers
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <core/pack/task/operation/TaskOperationCreator.fwd.hh>
#include <core/pack/task/operation/ResLvlTaskOperationCreator.fwd.hh>
#include <core/pack/task/operation/ResFilterCreator.fwd.hh>

// Basic headers
#include <basic/datacache/DataMap.fwd.hh>

// Utility Headers
#include <utility/SingletonBase.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.fwd.hh>

// c++ headers
#include <string>
#include <map>

#include <utility/vector0.hh>

namespace core {
namespace pack {
namespace task {
namespace operation {

// singleton class
class TaskOperationFactory: public utility::SingletonBase< TaskOperationFactory >
{
public:
	friend class utility::SingletonBase< TaskOperationFactory >;

	typedef utility::vector1< TaskOperationOP > TaskOperationOPs;
	typedef std::map< std::string, TaskOperationCreatorOP > TaskOperationCreatorMap;
	typedef utility::tag::Tag Tag;
	typedef utility::tag::TagOP TagOP;
	typedef utility::tag::TagCOP TagCOP;

public:
	void factory_register( TaskOperationCreatorOP );

	/// @brief add a prototype, using its default type name as the map key
	void add_creator( TaskOperationCreatorOP );
	bool has_type( std::string const & ) const;
	/// @brief pass through to child factories
	void add_creator( ResLvlTaskOperationCreatorOP );
	void add_creator( ResFilterCreatorOP );
	/// @brief return new TaskOperation by key lookup in task_operation_creator_map_ (new TaskOperation parses Tag if provided)
	TaskOperationOP newTaskOperation(
		std::string const &,
		basic::datacache::DataMap & datamap,
		TagCOP = utility::tag::TagCOP( utility::tag::TagOP( new utility::tag::Tag() ) )
	) const;
	/// @brief fills vector with new TaskOperations from nested "TASKOPERATIONS" TagCOP
	void newTaskOperations( TaskOperationOPs &, basic::datacache::DataMap & datamap, TagCOP ) const;
	/// @brief fills vector with new TaskOperations from xml-like tag file
	void newTaskOperations( TaskOperationOPs &, basic::datacache::DataMap & datamap, std::string const & ) const;


	/// @brief The %TaskOperationFactory is the point of entry for the definition of the XML Schemas
	/// for every TaskOperation that may be instantiated from a file. It is  responsible for defining
	/// an xs:group named "task_operation" listing each of the task-operation-complex types that
	/// may be initialized using the %TaskOperationFactory and to iterate across each of the
	/// TaskOperationCreator s it contains asking them for the XML schema of the TaskOperation they
	/// are responsible for creating.
	void define_task_op_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;

	/// @brief Read access to the map of creator names to creators -- for unit testing purposes only
	TaskOperationCreatorMap const & creator_map() const;

	static std::string task_operation_xml_schema_group_name();

private:
	// private constructor/destructor
	TaskOperationFactory();
	virtual ~TaskOperationFactory();

	TaskOperationFactory( TaskOperationFactory const &) = delete;
	TaskOperationFactory & operator=( TaskOperationFactory const &) = delete;

private:

	TaskOperationCreatorMap task_operation_creator_map_;

};

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core

#endif
