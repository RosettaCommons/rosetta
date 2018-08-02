// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/operation/OperateOnResidueSubset.cc
/// @brief  Class, much like the OperateOnResidueSubset task operation, to apply a particular
///         residue-level task operation on the residues identified by a ResidueSelector.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com) - Generalized for any vector1< bool >, flip_subset.

// Unit Headers
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/pack/task/operation/OperateOnResidueSubsetCreator.hh>


#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/ResLvlTaskOperation.hh>
#include <core/pack/task/operation/ResLvlTaskOperationFactory.hh>
#include <core/pack/task/operation/parsing_utilities.hh>
#include <core/pack/task/operation/task_op_schemas.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh>

#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

// Utility Headers
#include <utility>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace task {
namespace operation {

using basic::t_info;
using basic::t_debug;
static basic::Tracer TR( "core.pack.task.operation.OperateOnResidueSubset", t_info );

OperateOnResidueSubset::OperateOnResidueSubset() = default;

OperateOnResidueSubset::OperateOnResidueSubset(
	ResLvlTaskOperationCOP rlto,
	core::select::residue_selector::ResidueSelectorCOP selector,
	bool flip_subset
)
: parent(),
	residue_selector_(std::move( selector )),
	flip_subset_(flip_subset)
{
	ops_.push_back( std::move(rlto) );
}

OperateOnResidueSubset::OperateOnResidueSubset(
	ResLvlTaskOperationCOP rlto,
	utility::vector1< bool > const & subset
)
: parent(),
	user_provided_subset_( subset )
{
	ops_.push_back( std::move(rlto) );
}

OperateOnResidueSubset::OperateOnResidueSubset( OperateOnResidueSubset const & src )
: TaskOperation( src )
{
	*this = src;
}

OperateOnResidueSubset &
OperateOnResidueSubset::operator = ( OperateOnResidueSubset const & ) = default;

OperateOnResidueSubset::~OperateOnResidueSubset() = default;

TaskOperationOP OperateOnResidueSubset::clone() const
{
	return TaskOperationOP( new OperateOnResidueSubset( *this ) );
}

void
OperateOnResidueSubset::apply( Pose const & pose, PackerTask & ptask ) const
{
	Size const nres( pose.size() );
	runtime_assert( nres == ptask.total_residue() );

	core::select::residue_selector::ResidueSubset subset;

	if ( user_provided_subset_.size() > 0 ) {
		subset = user_provided_subset_;
	} else {
		subset = residue_selector_->apply( pose );
	}


	// Take the opposite of what the selection has chosen.
	if ( flip_subset_ ) {
		subset.flip();
	}

	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		if ( subset[ ii ] ) {
			for ( auto const & op : ops_ ) {
				op->apply( ptask.nonconst_residue_task( ii ) );
			}
		}
	}
}

void OperateOnResidueSubset::op( ResLvlTaskOperationCOP op_in )
{
	ops_.clear();
	ops_.push_back( op_in );
}

void OperateOnResidueSubset::append_op( ResLvlTaskOperationCOP op_in )
{
	ops_.push_back( op_in );

}

void OperateOnResidueSubset::flip_subset( bool flip_subset ){
	flip_subset_ = flip_subset;
}

void OperateOnResidueSubset::selector( core::select::residue_selector::ResidueSelectorCOP rs_in )
{
	residue_selector_ = rs_in;
}

void OperateOnResidueSubset::subset(utility::vector1<bool> const & subset_residues){
	user_provided_subset_ = subset_residues;
}


/// @brief tag parsing for factory construction of this class and its children
/*!
Example Tag syntax for parser as of October 2013

<OperateOnResidueSubset name=disable1to10>
<PreventRepackingRLT/>
<Index resnums=1-10/>
</OperateOnResidueSubset>

*/
void OperateOnResidueSubset::parse_tag( TagCOP tag , DataMap & datamap )
{
	using namespace core::select::residue_selector;

	std::string selector_name;
	ResidueSelectorCOP selector;
	if ( tag->hasOption( "selector_logic" ) ) {
		try {
			selector = parse_residue_selector_logic_string( datamap, tag );
		} catch ( utility::excn::Exception & e ) {
			std::ostringstream oss;
			oss << "Error in parsing selector_logic attribute of OperateOnResidueSubset tag";
			if ( tag->hasOption( "name" ) ) {
				oss << " with name \"" << tag->getOption< std::string >( "name" );
			}
			oss << "\nGenerated the following error message:\n" << e.msg();
			throw CREATE_EXCEPTION( utility::excn::Exception, oss.str() );
		}
	} else if ( tag->hasOption( "selector" ) ) {
		selector_name = tag->getOption< std::string >( "selector" );
		try {
			selector = datamap.get_ptr< ResidueSelector const >( "ResidueSelector", selector_name );
		} catch ( utility::excn::Exception & e ) {
			std::string error_message = "Failed to find ResidueSelector named '" + selector_name + "' from the Datamap from OperateOnResidueSubset::parse_tag\n" + e.msg();
			throw CREATE_EXCEPTION(utility::excn::Exception,  error_message );
		}
	}

	if ( tag->hasOption( "residue_level_operations" ) ) {
		parse_residue_level_task_operations( tag, datamap, ops_ );
	}

	utility::vector0< TagCOP > const & subtags( tag->getTags() );
	for ( auto const & subtag : subtags ) {
		std::string const & type( subtag->getName() );
		ResLvlTaskOperationFactory * rltof = ResLvlTaskOperationFactory::get_instance();
		if ( rltof && rltof->has_type( type ) ) {
			ResLvlTaskOperationOP rltop = rltof->newRLTO( type );
			rltop->parse_tag( subtag );
			ops_.push_back( rltop );
			TR(t_debug) << "Adding ResLvlTaskOperation of type " << type << std::endl;
			continue;
		}
		ResidueSelectorFactory * res_selector_factory = ResidueSelectorFactory::get_instance();
		if ( res_selector_factory && res_selector_factory->has_type( type ) ) {
			if ( selector ) {
				std::string error_message = "Error from OperateOnResidueSubset::parse_tag: multiple residue selectors given; only one allowed\n";
				throw CREATE_EXCEPTION(utility::excn::Exception,  error_message );
			}
			selector = res_selector_factory->new_residue_selector( subtag->getName(), subtag, datamap );
			TR(t_debug) << "using ResidueSelector of type " << type << std::endl;
			continue;
		}
		TR.Warning << type << " is not known to either the ResLvelTaskOperationFactory or the ResidueSelectorFactory. Ignorning it." << std::endl;
	}

	if ( ops_.empty() ) {
		std::string error_message = "Error from OperateOnResidueSubset::parse_tag.  No residue-level task operation provided\n";
		error_message += "The residue level task operation(s) must be provided either as a subtag of the OperateOnResidueSubset tag or with the residue_level_operations attribute\n";
		throw CREATE_EXCEPTION(utility::excn::Exception,  error_message );
	}
	if ( ! selector ) {
		std::string error_message = "Error from OperateOnResidueSubset::parse_tag.  No residue selector provided\n";
		error_message += "The residue selector may be provided either as a subtag of the OperateOnResidueSubset tag\n";
		error_message += "or through the 'selector' option.\n";
		throw CREATE_EXCEPTION(utility::excn::Exception,  error_message );
	}

	residue_selector_ = selector;
}


std::string OperateOnResidueSubset::keyname() { return "OperateOnResidueSubset"; }

/// @details The XSD says that the residue selector must appear before the RLTO.
void OperateOnResidueSubset::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	ResidueSelectorFactory::get_instance()->define_residue_selector_xml_schema( xsd );
	ResLvlTaskOperationFactory::get_instance()->define_res_lvl_task_op_xml_schema( xsd );

	AttributeList attlist;
	attlist + optional_name_attribute();
	attlist + XMLSchemaAttribute(
		"selector", xs_string,
		"Residue selector that indicates to which residues the operation will be applied.");
	select::residue_selector::attributes_for_parse_residue_selector_logic_string( xsd, attlist );
	attributes_for_parse_residue_level_operations( attlist );

	XMLSchemaComplexTypeGenerator ct_gen;
	XMLSchemaSimpleSubelementList rs_subelement, rlto_subelement;
	// the ResidueSelector subelement is not required; it can be provided by name through the datamap; thus
	// the min_occurs for this subelement is set to 0
	rs_subelement.add_group_subelement( & ResidueSelectorFactory::residue_selector_xml_schema_group_name );
	rlto_subelement.add_group_subelement( & ResLvlTaskOperationFactory::res_lvl_task_op_xml_schema_group_name );
	ct_gen.element_name( keyname() )
		.description(
		"OperateOnResidueSubset is a TaskOperation that applies one or more ResLevelTaskOperations "
		"to the residues indicated by a ResidueSelector. The ResidueSelector can either be previously "
		"defined, or constructed as a subtag beneath the main tag, or can be composed using boolean "
		"AND, OR, and ! (not) operations using the selector_logic attribute." )
		.complex_type_naming_func( & complex_type_name_for_task_op )
		.add_ordered_subelement_set_as_optional( rs_subelement )
		.add_ordered_subelement_set_as_repeatable( rlto_subelement )
		.add_attributes( attlist )
		.write_complex_type_to_schema( xsd );
}

TaskOperationOP OperateOnResidueSubsetCreator::create_task_operation() const
{
	return TaskOperationOP( new OperateOnResidueSubset );
}

std::string OperateOnResidueSubsetCreator::keyname() const { return OperateOnResidueSubset::keyname(); }

void OperateOnResidueSubsetCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	OperateOnResidueSubset::provide_xml_schema( xsd );
}

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core
