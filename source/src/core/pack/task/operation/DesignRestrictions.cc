// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/task/operation/DesignRestrictions.cc
/// @brief Combine a set of ResidueSelectors and ResidueLevelTaskOperations concisely, using boolean logic to join selectors concisely
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#include <core/pack/task/operation/DesignRestrictions.hh>
#include <core/pack/task/operation/DesignRestrictionsCreator.hh>

#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/ResLvlTaskOperation.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/task_op_schemas.hh>
#include <core/pack/task/operation/parsing_utilities.hh>
#include <core/pack/task/PackerTask.hh>

// Project headers
#include <core/chemical/AA.hh>
#include <core/pose/Pose.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

static basic::Tracer TR( "core.pack.task.operation.DesignRestrictions" );

namespace core {
namespace pack {
namespace task {
namespace operation {

using namespace core::pack::task::operation;

DesignRestrictions::DesignRestrictions():
	TaskOperation()
{}

DesignRestrictions::~DesignRestrictions() {}

TaskOperationOP
DesignRestrictions::clone() const {
	return TaskOperationOP( new DesignRestrictions( *this ) );
}

DesignRestrictions::DesignRestrictions( DesignRestrictions const & src ):
	TaskOperation(src),
	actions_( src.actions_ )
{}

void
DesignRestrictions::add_selector_rlto_pair(
	select::residue_selector::ResidueSelectorCOP selector,
	ResLvlTaskOperationCOP rlto
)
{
	actions_.push_back( { selector, { rlto } } );
}


void
DesignRestrictions::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & dm )
{
	using namespace select::residue_selector;
	for ( auto const & subtag : tag->getTags() ) {
		ResidueSelectorCOP selector;
		if ( subtag->hasOption( "selector_logic" ) ) {
			try {
				selector = parse_residue_selector_logic_string( dm, subtag );
			} catch ( utility::excn::Exception & e ) {
				std::ostringstream oss;
				oss << "Error in parsing selector_logic attribute of Action subtag of DesignRestriction tag";
				if ( tag->hasOption( "name" ) ) {
					oss << " with name \"" << tag->getOption< std::string >( "name" );
				}
				oss << "\nGenerated the following error message:\n" << e.msg();
				throw CREATE_EXCEPTION( utility::excn::Exception, oss.str() );
			}
		} else if ( subtag->hasOption( "residue_selector" ) ) {
			selector = parse_residue_selector( subtag, dm );
		} else {
			std::ostringstream oss;
			oss << "Error in parsing DesignRestriction task operation Action subtag: one of selector_logic or residue_selector attributes must be provided\n";
			throw CREATE_EXCEPTION( utility::excn::Exception, oss.str() );
		}

		std::list< ResLvlTaskOperationCOP > rltos;
		if ( subtag->hasOption( "aas" ) ) {
			auto restrict_aas( utility::pointer::make_shared< RestrictAbsentCanonicalAASRLT >() );
			std::string aas = subtag->getOption< std::string >( "aas" );
			for ( char aa : aas ) {
				char upper_aa = toupper( aa );
				if ( ! chemical::oneletter_code_specifies_aa( upper_aa ) ) {
					std::ostringstream oss;
					oss << "Neither '" << aa << "' nor its upper case equivalent ('" << upper_aa << "') are";
					oss << " valid identifiers for the twenty canonical amino acids.";
					throw CREATE_EXCEPTION( utility::excn::Exception, oss.str() );
				}
				if ( ! chemical::is_canonical_L_aa_or_gly( chemical::aa_from_oneletter_code( upper_aa ) ) ) {
					std::ostringstream oss;
					oss << "The one letter code '" << aa << "' (upper case converted to '" << upper_aa << "')";
					oss << " does not specify one of the twenty canonical AAs and cannot be used in the";
					oss << " DesignRestrictions \"aas\" attribute";
					throw CREATE_EXCEPTION( utility::excn::Exception, oss.str() );
				}
			}
			restrict_aas->aas_to_keep( aas );
			rltos.push_back( restrict_aas );
		}
		if ( subtag->hasOption( "residue_level_operations" ) ) {
			pack::task::operation::parse_residue_level_task_operations( subtag, dm, rltos );
		}

		if ( rltos.empty() ) {
			std::ostringstream oss;
			oss << "One or both of the \"aas\" or \"residue_level_operation\" attributes must be provided to the Action subtag of the DesignRestrictions task operation";
			throw CREATE_EXCEPTION( utility::excn::Exception, oss.str() );
		}
		actions_.push_back( { selector, rltos } );
	}
}



void
DesignRestrictions::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	for ( auto const & selector_and_rlto : actions_ ) {
		utility::vector1< bool > selection = selector_and_rlto.first->apply( pose );
		for ( Size ii = 1; ii <= selection.size(); ++ii ) {
			if ( selection[ ii ] ) {
				for ( auto const & rlto : selector_and_rlto.second ) {
					rlto->apply( task.nonconst_residue_task( ii ) );
				}
			}
		}
	}
}

std::string
DesignRestrictions::keyname() {
	return "DesignRestrictions";
}

void
DesignRestrictions::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	using namespace select::residue_selector;
	AttributeList action_attlist;
	attributes_for_parse_residue_selector_default_option_name( action_attlist );
	attributes_for_parse_residue_selector_logic_string( xsd, action_attlist );
	action_attlist
		+ XMLSchemaAttribute( "aas", xs_string, "A list of the canonical L amino acids that are to be"
		" allowed for the selected residues. Either upper or lower case letters are accepted. Note"
		" that if an amino acid is allowed for a residue that is selected by one action, but disallowed"
		" for that residue selected by a second action, the amino acid will not be allowed." );
	attributes_for_parse_residue_level_operations( action_attlist );

	XMLSchemaSimpleSubelementList action_subelement;
	action_subelement.add_simple_subelement( "Action", action_attlist, "Select a set of residues and then perform a set of actions on those residues" );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & complex_type_name_for_task_op )
		.element_name( keyname() )
		.description( "Concise way to combine residue selectors and residue-level task operations. The Action subtags"
		" allow you to specify several modifications to a packer task simultaneously. Using the"
		" 'selector_logic' attribute, you can combine as many previously-defined ResidueSelectors as you"
		" would like using AND, OR, ! (not), and parentheses or you can choose a particular, previously defined ResidueSelector"
		" using the 'residue_selector' attribute. One of the two must be provided. Then you may specify the set of"
		" canonical amino acids you'd like to restrict the indicated set of residues using the 'aas' attribute and/or"
		" by listing a set of previously-declared ResidueLevelTaskOperations." )
		.set_subelements_repeatable( action_subelement )
		.add_attribute( optional_name_attribute() )
		.write_complex_type_to_schema( xsd );
}

TaskOperationOP
DesignRestrictionsCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new DesignRestrictions );
}

std::string
DesignRestrictionsCreator::keyname() const
{
	return DesignRestrictions::keyname();
}

void
DesignRestrictionsCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DesignRestrictions::provide_xml_schema( xsd );
}


} //core
} //pack
} //task
} //operation
