// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file basic/database/schema_generator/DataType.hh
///
/// @brief DbDataType class for the schema generator framework
/// @author Tim Jacobs

#ifndef INCLUDED_basic_database_schema_generator_DbDataType_HH
#define INCLUDED_basic_database_schema_generator_DbDataType_HH

#include <basic/database/schema_generator/DbDataType.fwd.hh>
#include <string>                                              // for string
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>  // for sessionOP
#include <utility/pointer/ReferenceCount.hh>

namespace basic {
namespace database {
namespace schema_generator {

//Class definitions for datatypes
class DbDataType : public utility::pointer::ReferenceCount {

public:
	DbDataType();
	DbDataType(int size);

	virtual
	std::string print(
		utility::sql_database::sessionOP
	) const = 0;

protected:
	int size_;
	std::string type_string_;
};

//General text data type
class DbText : public  DbDataType {
public:
	DbText();
	DbText(int size);
	std::string print(utility::sql_database::sessionOP) const;
};

class DbTextKey : public  DbDataType {
public:
	DbTextKey();
	std::string print(utility::sql_database::sessionOP) const;
};

//General integer data type
class DbInteger : public DbDataType {
public:
	DbInteger();
	std::string print(utility::sql_database::sessionOP) const;
};

/// DOES NOT WORK WITH CPPDB - USE INTEGER INSTEAD
////boolean data type
//class DbBoolean : public DbDataType {
//public:
// DbBoolean();
//
//};

class DbBigInt : public DbDataType {
public:
	DbBigInt();
	std::string print(utility::sql_database::sessionOP) const;
};

//Double data type
class DbDouble : public DbDataType {
public:
	DbDouble();
	std::string print(utility::sql_database::sessionOP) const;
};

//Real data type
class DbReal : public DbDataType {
public:
	DbReal();
	std::string print(utility::sql_database::sessionOP) const;
};

//Struct Id has its own type due to incompatibilities between backends in regards to unsigned long longs
class DbUUID : public DbDataType {
public:
	DbUUID();
	std::string print(utility::sql_database::sessionOP) const;
};

} // schema_generator
} // namespace database
} // namespace utility

#endif
