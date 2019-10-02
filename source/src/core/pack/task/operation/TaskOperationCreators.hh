// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/operation/TaskOperationCreators.hh
/// @brief  Declaration for the class that connects TaskOperations with the TaskOperationFactory
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author ashworth

#ifndef INCLUDED_core_pack_task_operation_TaskOperationCreators_hh
#define INCLUDED_core_pack_task_operation_TaskOperationCreators_hh

#include <core/pack/task/operation/TaskOperationCreator.hh>

#include <core/pack/task/operation/TaskOperation.fwd.hh>

#include <string>


namespace core {
namespace pack {
namespace task {
namespace operation {

class RestrictYSDesignCreator : public TaskOperationCreator {
public:
	TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class PreventRepackingCreator : public TaskOperationCreator {
public:
	TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class PreserveCBetaCreator : public TaskOperationCreator {
public:
	TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class AppendRotamerSetCreator : public TaskOperationCreator {
public:
	TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class AppendResidueRotamerSetCreator : public TaskOperationCreator {
public:
	TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class AppendRotamerCreator : public TaskOperationCreator {
public:
	TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class SetRotamerCouplingsCreator : public TaskOperationCreator {
public:
	TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class SetRotamerLinksCreator : public TaskOperationCreator {
public:
	TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};


class ReadResfileCreator : public TaskOperationCreator {
public:
	TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class ReadResfileAndObeyLengthEventsCreator : public TaskOperationCreator {
public:
	TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};


class IncludeCurrentCreator : public TaskOperationCreator {
public:
	TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};


class InitializeExtraRotsFromCommandlineCreator : public TaskOperationCreator {
public:
	TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class InitializeFromCommandlineCreator : public TaskOperationCreator {
public:
	TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class UseMultiCoolAnnealerCreator : public TaskOperationCreator {
public:
	TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class InitializeFromOptionCollectionCreator : public TaskOperationCreator {
public:
	TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class ExtraRotamersGenericCreator : public TaskOperationCreator {
public:
	TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class RotamerExplosionCreator : public TaskOperationCreator {
public:
	TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class RestrictAbsentCanonicalAASCreator : public TaskOperationCreator {
public:
	TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};


class RestrictToSpecifiedBaseResidueTypesCreator : public TaskOperationCreator {
	TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class ProhibitSpecifiedBaseResidueTypesCreator : public TaskOperationCreator {
	TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class RestrictToResiduePropertiesCreator : public TaskOperationCreator {
	TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class ProhibitResiduePropertiesCreator : public TaskOperationCreator {
	TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class RestrictResidueToRepackingCreator : public TaskOperationCreator {
public:
	TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class RestrictToRepackingCreator : public TaskOperationCreator {
public:
	TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class DisallowIfNonnativeCreator : public TaskOperationCreator {
public:
	TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class ExtraRotamersCreator : public TaskOperationCreator {
public:
	TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class ExtraChiCutoffCreator : public TaskOperationCreator {
public:
	TaskOperationOP create_task_operation() const override;
	std::string keyname() const override { return "ExtraChiCutoff"; }
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core

#endif
