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
#include <core/pack/task/operation/task_op_schemas.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>

#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

// Utility Headers
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
static THREAD_LOCAL basic::Tracer TR( "core.pack.task.operation.OperateOnResidueSubset", t_info );

OperateOnResidueSubset::OperateOnResidueSubset()
: parent(),
	op_(/* 0 */),
	residue_selector_(/* 0 */),
	flip_subset_(false)
{}

OperateOnResidueSubset::OperateOnResidueSubset(
	ResLvlTaskOperationCOP rlto,
	core::select::residue_selector::ResidueSelectorCOP selector,
	bool flip_subset
)
: parent(),
	op_( rlto ),
	residue_selector_( selector ),
	flip_subset_(flip_subset)
{}

OperateOnResidueSubset::OperateOnResidueSubset(
	ResLvlTaskOperationCOP rlto,
	utility::vector1< bool > const & subset
)
: parent(),
	op_( rlto ),
	user_provided_subset_( subset )
{}

OperateOnResidueSubset::OperateOnResidueSubset( OperateOnResidueSubset const & src )
: TaskOperation( src )
{
	*this = src;
}

OperateOnResidueSubset &
OperateOnResidueSubset::operator = ( OperateOnResidueSubset const & src )
{
	if ( this != & src ) {
		op_ = src.op_;
		residue_selector_ = src.residue_selector_;
		user_provided_subset_ = src.user_provided_subset_;
		flip_subset_ = src.flip_subset_;
	}
	return *this;
}

OperateOnResidueSubset::~OperateOnResidueSubset() {}

TaskOperationOP OperateOnResidueSubset::clone() const
{
	return TaskOperationOP( new OperateOnResidueSubset( *this ) );
}

void
OperateOnResidueSubset::apply( Pose const & pose, PackerTask & ptask ) const
{
	Size const nres( pose.total_residue() );
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

	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		if ( subset[ ii ] ) op_->apply( ptask.nonconst_residue_task( ii ) );
	}
}

void OperateOnResidueSubset::op( ResLvlTaskOperationCOP op_in )
{
	op_ = op_in;
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

	ResLvlTaskOperationCOP rltop;

	std::string selector_name;
	ResidueSelectorCOP selector;
	if ( tag->hasOption( "selector" ) ) {
		selector_name = tag->getOption< std::string >( "selector" );
		try {
			selector = datamap.get_ptr< ResidueSelector const >( "ResidueSelector", selector_name );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::string error_message = "Failed to find ResidueSelector named '" + selector_name + "' from the Datamap from OperateOnResidueSubset::parse_tag\n" + e.msg();
			throw utility::excn::EXCN_Msg_Exception( error_message );
		}
	}

	utility::vector0< TagCOP > const & subtags( tag->getTags() );
	for ( utility::vector0< TagCOP >::const_iterator
			subtag( subtags.begin() ), end( subtags.end() );
			subtag != end; ++subtag ) {
		std::string const type( (*subtag)->getName() );
		ResLvlTaskOperationFactory * rltof = ResLvlTaskOperationFactory::get_instance();
		if ( rltof && rltof->has_type( type ) ) {
			if ( rltop ) {
				std::string error_message = "Error from OperateOnResidueSubset::parse_tag: multiple residue-level task operations given; only one allowed\n";
				throw utility::excn::EXCN_Msg_Exception( error_message );
			}
			rltop = rltof->newRLTO( type );
			rltop->parse_tag( *subtag );
			TR(t_debug) << "using ResLvlTaskOperation of type " << type << std::endl;
			continue;
		}
		ResidueSelectorFactory * res_selector_factory = ResidueSelectorFactory::get_instance();
		if ( res_selector_factory && res_selector_factory->has_type( type ) ) {
			if ( selector ) {
				std::string error_message = "Error from OperateOnResidueSubset::parse_tag: multiple residue selectors given; only one allowed\n";
				throw utility::excn::EXCN_Msg_Exception( error_message );
			}
			selector = res_selector_factory->new_residue_selector( (*subtag)->getName(), *subtag, datamap );
			TR(t_debug) << "using ResidueSelector of type " << type << std::endl;
			continue;
		}
		TR.Warning << "WARNING: " << type << " is not known to either the ResLvelTaskOperationFactory or the ResidueSelectorFactory. Ignorning it." << std::endl;
	}

	if ( ! rltop ) {
		std::string error_message = "Error from OperateOnResidueSubset::parse_tag.  No residue-level task operation provided\n";
		error_message += "The residue level task operation must be provided as a subtag of the OperateOnResidueSubset tag\n";
		throw utility::excn::EXCN_Msg_Exception( error_message );
	}
	if ( ! selector ) {
		std::string error_message = "Error from OperateOnResidueSubset::parse_tag.  No residue selector provided\n";
		error_message += "The residue selector may be provided either as a subtag of the OperateOnResidueSubset tag\n";
		error_message += "or through the 'selector' option.\n";
		throw utility::excn::EXCN_Msg_Exception( error_message );
	}

	op_ = rltop;
	residue_selector_ = selector;
}


std::string OperateOnResidueSubset::keyname() { return "OperateOnResidueSubset"; }

/// @details The XSD says that the residue selector must appear before the RLTO.
void OperateOnResidueSubset::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	ResidueSelectorFactory::get_instance()->define_residue_selector_xml_schema( xsd );
	ResLvlTaskOperationFactory::get_instance()->define_res_lvl_task_op_xml_schema( xsd );

	using namespace utility::tag;
	XMLSchemaComplexTypeGenerator ct_gen;
	XMLSchemaSimpleSubelementList subelements;
	// the ResidueSelector subelement is not required; it can be provided by name through the datamap; thus
	// the min_occurs for this subelement is set to 0
	subelements.add_group_subelement( & ResidueSelectorFactory::residue_selector_xml_schema_group_name, 0 /*min_occurs*/ );
	subelements.add_group_subelement( & ResLvlTaskOperationFactory::res_lvl_task_op_xml_schema_group_name );
	ct_gen.element_name( keyname() )
		.complex_type_naming_func( & complex_type_name_for_task_op )
		.set_subelements_single_appearance_required_and_ordered( subelements )
		.add_attribute( optional_name_attribute() )
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
