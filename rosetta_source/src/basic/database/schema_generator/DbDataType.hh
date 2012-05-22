// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file DataType.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_basic_database_schema_generator_DbDataType_HH
#define INCLUDED_basic_database_schema_generator_DbDataType_HH

//C++
#include <string>

namespace basic{
namespace database{
namespace schema_generator{

//Class definitions for datatypes
class DbDataType {

public:
	DbDataType();
	std::string print() const;

protected:
	std::string database_mode_;
	std::string type_string_;
};

//General text data type
class DbText : public  DbDataType {
public:
	DbText();

	DbText(int size);
};

class DbTextKey : public  DbDataType {
public:
	DbTextKey();
};

//General integer data type
class DbInteger : public DbDataType {
public:
	DbInteger();

	DbInteger(int size);
};

/// DOES NOT WORK WITH CPPDB - USE INTEGER INSTEAD
////boolean data type
//class DbBoolean : public DbDataType {
//public:
//	DbBoolean();
//
//};

class DbBigInt : public DbDataType {
public:
	DbBigInt();

};
	
//Double data type
class DbDouble : public DbDataType {
public:
	DbDouble();
};
	
//Real data type
class DbReal : public DbDataType {
public:
	DbReal();
};

//Struct Id has its own type due to incompatibilities between backends in regards to unsigned long longs
class DbUUID : public DbDataType {
public:
	DbUUID();
};

} // schema_generator
} // namespace database
} // namespace utility

#endif
