// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/loops_definers/LoopsDatabaseDefiner.cc
/// @brief  A loops definer is creates a serialized loops list
/// @author Matthew O'Meara (mattjomear@gmail.com)

// Unit Headers
#include <protocols/loops/loops_definers/LoopsDatabaseDefiner.hh>

// Package Headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/loops_definers/util.hh>

// Project Headers
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

// Basic Headers
#include <basic/database/sql_utils.hh>
#include <basic/datacache/DataMap.hh>

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <string>
#include <utility/excn/Exceptions.hh>
#include <sstream>


using std::string;
using std::endl;
using std::stringstream;
using utility::tag::TagCOP;
using core::pose::Pose;
using basic::datacache::DataMap;
using utility::sql_database::sessionOP;
using protocols::jd2::JobDistributor;
using basic::database::parse_database_connection;
using cppdb::result;

namespace protocols {
namespace loops {
namespace loops_definers {

LoopsDatabaseDefiner::LoopsDatabaseDefiner() :
	db_session_(),
	database_table_()
{}

LoopsDatabaseDefiner::~LoopsDatabaseDefiner() {}

LoopsDatabaseDefiner::LoopsDatabaseDefiner(LoopsDatabaseDefiner const & src) : LoopsDefiner(src),
	db_session_(src.db_session_),
	database_table_(src.database_table_)
{}


/// @brief Create another loops definer of the type matching the most-derived
/// version of the class.
LoopsDefinerOP
LoopsDatabaseDefiner::clone(
) const {
	return LoopsDefinerOP( new LoopsDatabaseDefiner(*this) );
}

/// @brief Used to parse an xml-like tag to load parameters and properties.
void
LoopsDatabaseDefiner::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap const &,
	Pose const &
) {

	db_session_ = parse_database_connection(tag);

	database_table_ =
		tag->getOption<std::string>("database_table", "loops");

	// the name attribute, as Matt sees it, is not optional. Its absence will
	// cause the following code to throw an exception. However, the name is actually
	// not needed by this class, but rather, by the LoopDefinerLoader. The following
	// code could be reasonably removed.
	if ( !tag->hasOption("name") ) {
		throw utility::excn::EXCN_RosettaScriptsOption(
			"Unable to create unnamed LoopsDefiner (type: " + tag->getName() + ")" );
	}
	string const loops_name(tag->getOption<string>("name"));

}

SerializedLoopList
LoopsDatabaseDefiner::apply(
	Pose const &
) {

	string pose_tag(JobDistributor::get_instance()->current_job()->input_tag());

	stringstream sql_stmt;
	sql_stmt
		<< "SELECT start, stop, cut, skip_rate, extended FROM " << database_table_
		<< " WHERE tag='" << pose_tag << "';";
	result res = (*db_session_) << sql_stmt.str();

	SerializedLoopList loop_list;
	while ( res.next() ) {
		SerializedLoop loop;
		int extended;
		res >> loop.start >> loop.stop >> loop.cut >> loop.skip_rate >> extended;
		loop.extended = extended;

		loop_list.push_back(loop);
	}

	if ( loop_list.size() == 0 ) {
		stringstream error_message;
		error_message
			<< "Unable to locate loops for job distributor input tag '"
			<< pose_tag << "' in database." << endl;
		utility_exit_with_message(error_message.str());
	}

	return loop_list;
}

std::string LoopsDatabaseDefiner::class_name()
{
	return "LoopsDatabase";
}

void LoopsDatabaseDefiner::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::attribute_w_default( "database_table", xs_string,
		"The name of the table in the database from which the loops are to be read", "loops" )
		+ required_name_attribute();
	basic::database::attributes_for_parse_database_connection( attributes, xsd );

	xsd_type_definition_w_attributes( xsd, class_name(), "Load the loop defintions from a table in a database; the table from"
		" which this LoopsDefiner reads may be specified, but if it is not then it will look for a table named 'loops'."
		" Beause of the way the apply function is written, this class is specific for JD2, which is a shame; it will"
		" have to be updated for JD3.", attributes );
}



} //namespace
} //namespace
} //namespace


