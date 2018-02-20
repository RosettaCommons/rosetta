// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/parser/DataLoader.cc
/// @brief  Implementation of the XML parser's DataLoader base class (ctor & dstor)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <protocols/parser/ResLvlTaskOperationLoader.hh>
#include <protocols/parser/StandardLoaderCreators.hh>

// Project Headers
#include <core/pack/task/operation/ResLvlTaskOperation.hh>
#include <core/pack/task/operation/ResLvlTaskOperationFactory.hh>
#include <core/pack/task/operation/parsing_utilities.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

#include <basic/datacache/DataMap.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace parser {

static basic::Tracer TR( "protocols.jd2.parser.ResLvlTaskOperationLoader" );

ResLvlTaskOperationLoader::ResLvlTaskOperationLoader() = default;
ResLvlTaskOperationLoader::~ResLvlTaskOperationLoader() = default;

void ResLvlTaskOperationLoader::load_data(
	core::pose::Pose const &,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) const
{
	using namespace core::pack::task::operation;

	for ( utility::tag::TagCOP subtag : tag->getTags() ) {
		std::string const & type( subtag->getName() );
		if ( ! subtag->hasOption("name") ) {
			utility_exit_with_message( "Can't create unnamed ResLvlTaskOperation (type: " + type + ")" );
		}
		std::string const & name( subtag->getOption<std::string>("name") );
		if ( data.has( rlto_datamap_category(), name ) ) {
			TR.Error << "ResLvlTaskOperation of name \"" << name
				<< "\" (with type " << type << ") already exists. \n" << subtag << std::endl;
			utility_exit_with_message("Duplicate definition of ResLvlTaskOperation with name " + name);
		}
		ResLvlTaskOperationOP new_rlto( ResLvlTaskOperationFactory::get_instance()->newRLTO( type, subtag ) );
		runtime_assert( new_rlto != nullptr );
		data.add( core::pack::task::operation::rlto_datamap_category(), name, new_rlto );
		TR << "Defined ResLvlTaskOperation named \"" << name << "\" of type " << type << std::endl;
	}
	TR.flush();

}

std::string ResLvlTaskOperationLoader::loader_name() { return "RESIDUE_LEVEL_TASK_OPERATIONS"; }

std::string ResLvlTaskOperationLoader::res_lvl_task_op_loader_ct_namer( std::string const & element_name )
{
	return "res_lvl_task_op_loader_" + element_name + "_type";
}

void ResLvlTaskOperationLoader::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	using namespace core::pack::task::operation;

	ResLvlTaskOperationFactory::get_instance()->define_res_lvl_task_op_xml_schema( xsd );

	XMLSchemaSimpleSubelementList res_lvl_task_op_subelements;
	res_lvl_task_op_subelements.add_group_subelement( & ResLvlTaskOperationFactory::res_lvl_task_op_xml_schema_group_name );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.element_name( loader_name() )
		.description( "Load a set of operations that apply to the packing data for single residues -- the ResidueLevelTask objects held as an array in the PackerTask. These ResLvlTaskOperations must be combined with ResidueSelectors in order to be applied to particular residues." )
		.complex_type_naming_func( res_lvl_task_op_loader_ct_namer )
		.set_subelements_repeatable( res_lvl_task_op_subelements )
		.write_complex_type_to_schema( xsd );
}


DataLoaderOP
ResLvlTaskOperationLoaderCreator::create_loader() const { return DataLoaderOP( new ResLvlTaskOperationLoader ); }

std::string
ResLvlTaskOperationLoaderCreator::keyname() const { return ResLvlTaskOperationLoader::loader_name(); }

ResLvlTaskOperationLoaderCreator::DerivedNameFunction
ResLvlTaskOperationLoaderCreator::schema_ct_naming_function() const
{
	return & ResLvlTaskOperationLoader::res_lvl_task_op_loader_ct_namer;
}

void
ResLvlTaskOperationLoaderCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ResLvlTaskOperationLoader::provide_xml_schema( xsd );
}


} //namespace parser
} //namespace protocols
