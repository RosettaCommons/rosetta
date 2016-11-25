// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/rdf/RDFFunctionFactory.hh
/// @brief  headers for RDFFunctionFactory
/// @author Sam DeLuca


#ifndef INCLUDED_protocols_ligand_docking_rdf_RDFFunctionFactory_hh
#define INCLUDED_protocols_ligand_docking_rdf_RDFFunctionFactory_hh

// Unit Headers
#include <protocols/ligand_docking/rdf/RDFFunctionFactory.fwd.hh>

// Project Headers
#include <protocols/ligand_docking/rdf/RDFBase.fwd.hh>
#include <protocols/ligand_docking/rdf/RDFFunctionCreator.fwd.hh>

// Platform Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/factory/WidgetRegistrator.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/SingletonBase.hh>

// C++ Headers
#include <map>

#include <utility/vector1.hh>


namespace protocols {
namespace ligand_docking {
namespace rdf {

/// Create RDFFunctions
class RDFFunctionFactory : public utility::SingletonBase< RDFFunctionFactory > {
public:
	friend class utility::SingletonBase< RDFFunctionFactory >;

	// Private constructor to make it singleton managed
	RDFFunctionFactory();
	RDFFunctionFactory(const RDFFunctionFactory & src) = delete;

	RDFFunctionFactory &
	operator=( RDFFunctionFactory const & ) = delete;

public:

	// Warning this is not called because of the singleton pattern
	virtual ~RDFFunctionFactory();

	void factory_register( RDFFunctionCreatorCOP creator );
	RDFBaseOP get_rdf_function( std::string const & type_name );

	/// @brief convienence header for use with RosettaScripts parse_my_tag interface
	RDFBaseOP
	get_rdf_function(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data);

	utility::vector1<std::string> get_all_function_names();

	void define_rdf_function_group( utility::tag::XMLSchemaDefinition & xsd );
	static std::string rdf_function_group_name();
	static std::string rdf_function_ct_namer( std::string );
	static void xsd_type_definition_w_attributes( utility::tag::XMLSchemaDefinition & xsd, std::string name, utility::tag::AttributeList & attlist, std::string description );

private:

private:

	typedef std::map< std::string, RDFFunctionCreatorCOP > RDFFunctionCreatorMap;
	RDFFunctionCreatorMap types_;
};


/// @brief This templated class will register an instance of an
/// RDFFunctionCreator (class T) with the
/// RDFFunctionFactory.  It will ensure that no
/// FeaturesReporterCreator is registered twice, and, centralizes this
/// registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class RDFFunctionRegistrator : public utility::factory::WidgetRegistrator< RDFFunctionFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< RDFFunctionFactory, T > parent;
public:
	RDFFunctionRegistrator() : parent() {}
};


}
} // namespace
} // namespace

#endif
