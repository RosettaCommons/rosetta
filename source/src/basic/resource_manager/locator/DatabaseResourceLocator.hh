// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/locator/DatabaseResourceLocater.hh
/// @brief  Locate data in a datastore: locator_id (std::string) -> data (std::istream)
/// @author Matthew O'Meara
/// @author Andrew Leaver-Fay
/// @author Brian Weitzner

#ifndef INCLUDED_basic_resource_manager_locator_DatabaseResourceLocater_hh
#define INCLUDED_basic_resource_manager_locator_DatabaseResourceLocater_hh

//unit headers
#include <basic/resource_manager/ResourceLocator.hh>
#include <basic/resource_manager/locator/DatabaseResourceLocator.fwd.hh>

//project headers
#include <utility/tag/Tag.fwd.hh>

//C++ headers
#include <istream>

namespace basic {
namespace resource_manager {
namespace locator {

/// @brief The %DatabaseResourceLocator class is responsible for retreiving data
/// from a Database so that that data can then be used to construct a Resource.
///
/// @details  Upon construction or in its parse_my_tag method, the
/// %DatabaseResourceLocator needs to be given the name of the database-session
/// resource that it will use to communicate with the database and a partially
/// formatted SQL command (a SELECT statement) that will be used to query the
/// database for the resource that it will be pulling from the database.  This
/// SQL command should have a single question mark ("?") for the variable that
/// will be replaced by the locator_id that will be given to the
/// %DatabaseResourceLocator when its locate_resource_stream method is invoked.
/// In that method, the %DatabaseResourceLocator will ask the ResourceManager
/// for the database session object and then bind the input "locator_id" to the
/// previously provided SQL command, and finally will query the database with that
/// statement.  The resulting output is packaged in a string stream and returned
/// to the code that called locate_resource_stream  (e.g. the ResourceManager)
/// and then can be used to construct a resource.
class DatabaseResourceLocator : public basic::resource_manager::ResourceLocator
{
public:
	/// @brief The default constructor sets empty strings for both the database session
	/// and for the SQL command; using this constructor requires initializing those two
	/// pieces of data using the parse-my-tag method.
	DatabaseResourceLocator();

	/// @brief Constructor that takes in both the name of the database session resource
	/// as well as the SQL command
	DatabaseResourceLocator(
		std::string const & database_session_resource_tag,
		std::string const & sql_command);

	DatabaseResourceLocator(
		DatabaseResourceLocator const & src);

	virtual ~DatabaseResourceLocator();

	/// @brief Describe the %DatabaseResourceLocator instance to the given output stream
	virtual
	void
	show(
		std::ostream & out) const;

	/// @brief Return the typename for this class: "DatabaseResourceLocator"
	virtual
	std::string
	type() const;

	/// @brief Create a ResourceStream object from the given resource
	/// source, so that its stream can be passed to the ResourceLoader
	/// using the input "locator_tag" which is bound to the partially formed
	/// SQL select statement that was provided at construction or in
	/// parse_my_tag
	virtual
	ResourceStreamOP
	locate_resource_stream(
		std::string const & locator_tag
	) const;

	/// @brief Initialize the %DatabaseResourceLoader from the input set of tags
	/// which should contain both the name of the database session resource that
	/// will be used to talk to the database, and the partially formed SQL
	/// select statement that will be used to query the database.
	virtual
	void
	parse_my_tag(
		utility::tag::TagCOP tag
	);

private:

	std::string database_session_resource_tag_;
	std::string sql_command_;
	std::string column_separator_;
};

} // namespace locator
} // namespace resource_manager
} // namespace basic


#endif //INCLUDED_basic_resource_manager_locator_DatabaseResourceLocator_hh
