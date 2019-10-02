// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/operation/ResFilterCreators.hh
/// @brief  Declaration for the class that connects ResFilters with the ResFilterFactory
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author ashworth

#ifndef INCLUDED_core_pack_task_operation_ResFilterCreators_hh
#define INCLUDED_core_pack_task_operation_ResFilterCreators_hh

#include <core/pack/task/operation/ResFilterCreator.hh>

#include <core/pack/task/operation/ResFilter.fwd.hh>

#include <string>


namespace core {
namespace pack {
namespace task {
namespace operation {

class AnyResFilterCreator : public ResFilterCreator {
public:
	ResFilterOP create_res_filter() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class AllResFilterCreator : public ResFilterCreator {
public:
	ResFilterOP create_res_filter() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class NoResFilterCreator : public ResFilterCreator {
public:
	ResFilterOP create_res_filter() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class ResidueTypeFilterCreator : public ResFilterCreator {
public:
	ResFilterOP create_res_filter() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class ResidueHasPropertyCreator : public ResFilterCreator {
public:
	ResFilterOP create_res_filter() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class ResidueLacksPropertyCreator : public ResFilterCreator {
public:
	ResFilterOP create_res_filter() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class ResiduePDBInfoHasLabelCreator : public ResFilterCreator {
public:
	ResFilterOP create_res_filter() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class ResiduePDBInfoLacksLabelCreator : public ResFilterCreator {
public:
	ResFilterOP create_res_filter() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class ResidueName3IsCreator : public ResFilterCreator {
public:
	ResFilterOP create_res_filter() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class ResidueName3IsntCreator : public ResFilterCreator {
public:
	ResFilterOP create_res_filter() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class ResidueIndexIsCreator : public ResFilterCreator {
public:
	ResFilterOP create_res_filter() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class ResidueIndexIsntCreator : public ResFilterCreator {
public:
	ResFilterOP create_res_filter() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class ResiduePDBIndexIsCreator : public ResFilterCreator {
public:
	ResFilterOP create_res_filter() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class ResiduePDBIndexIsntCreator : public ResFilterCreator {
public:
	ResFilterOP create_res_filter() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class ChainIsCreator : public ResFilterCreator {
public:
	ResFilterOP create_res_filter() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class ChainIsntCreator : public ResFilterCreator {
public:
	ResFilterOP create_res_filter() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};


} //namespace operation
} //namespace task
} //namespace pack
} //namespace core

#endif
