// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file basic/database/insert_statement_generator/DataType.hh
///
/// @brief Column data object for the insert statement generator
/// @author Sam DeLuca

#ifndef INCLUDED_basic_database_insert_statement_generator_RowData_HH
#define INCLUDED_basic_database_insert_statement_generator_RowData_HH

#include <platform/types.hh>
#include <basic/database/insert_statement_generator/RowData.fwd.hh>
#include <boost/uuid/uuid.hpp>
#include <cppdb/frontend.h>

#include <utility/pointer/ReferenceCount.hh>

#include <string>

namespace basic {
namespace database {
namespace insert_statement_generator {

class RowDataBase : public utility::pointer::ReferenceCount {
public:

	RowDataBase(std::string const & column_name);
	virtual ~RowDataBase();

	std::string get_column_name() const;
	virtual void bind_data(platform::Size index, cppdb::statement & statement) = 0;


private:
	std::string column_name_;
};

template<class T>
class RowData : public RowDataBase {
public:
	RowData(
		std::string const & column_name,
		T const & data) : RowDataBase(column_name),data_(data)
	{}
	virtual void bind_data(
		platform::Size index,
		cppdb::statement & statement)
	{
		statement.bind(index,data_);
	}
private:
	T data_;
};


}
}
}


#endif /* ROWDATA_HH_ */
