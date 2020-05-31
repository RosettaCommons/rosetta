// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file --path--/--class--.cc
/// @brief --brief--
/// @author --name-- (--email--)

// Unit headers
#include <--path--/--class--.hh>
#include <--path--/--class--Creator.hh>

// Core headers
#include <core/pose/Pose.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// Citation Manager
#include <utility/vector1.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

static basic::Tracer TR( "--namespace_dot--.--class--" );

--namespace--

	/////////////////////
	/// Constructors  ///
	/////////////////////

/// @brief Default constructor
--class--::--class--():
	protocols::moves::Mover( --class--::mover_name() )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
--class--::~--class--(){}

////////////////////////////////////////////////////////////////////////////////
	/// Mover Methods ///
	/////////////////////

/// @brief Apply the mover
void
--class--::apply( core::pose::Pose& ){

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
--class--::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
--class--::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap&
) {

}
void --class--::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "DOCUMENTATION STRING", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
--class--::fresh_instance() const
{
	return utility::pointer::make_shared< --class-- >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
--class--::clone() const
{
	return utility::pointer::make_shared< --class-- >( *this );
}

std::string --class--::get_name() const {
	return mover_name();
}

std::string --class--::mover_name() {
	return "--class--";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
--class--Creator::create_mover() const
{
	return utility::pointer::make_shared< --class-- >();
}

std::string
--class--Creator::keyname() const
{
	return --class--::mover_name();
}

void --class--Creator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	--class--::provide_xml_schema( xsd );
}

/// @brief Indicate that this mover is unpublished.
bool
--class--::mover_is_unpublished() const {
	return true;
}

/// @brief Provide authorship information for an unpublished Rosetta module.
utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP >
--class--::provide_authorship_info_for_unpublished() const {
	using namespace basic::citation_manager;
	UnpublishedModuleInfoOP moduleinfo =
		utility::pointer::make_shared< UnpublishedModuleInfo >( "--class--", CitedModuleType::Mover );
	moduleinfo->add_author( "--name--", "TODO: institution", "--email--" );
	return utility::vector1< UnpublishedModuleInfoCOP > { moduleinfo };
}


////////////////////////////////////////////////////////////////////////////////
	/// private methods ///
	///////////////////////


std::ostream &
operator<<( std::ostream & os, --class-- const & mover )
{
	mover.show(os);
	return os;
}


--end_namespace--
