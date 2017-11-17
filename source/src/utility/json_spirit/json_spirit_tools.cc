// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/json_spirit/json_spirit_tools.cc
/// @brief  Extra utilitiy functions not provided by json spirit library
/// @author Mike Tyka

#include <utility/json_spirit/json_spirit_value.h>
#include <utility/json_spirit/json_spirit_reader.h>
#include <utility/json_spirit/json_spirit_tools.hh>
#include <utility/excn/Exceptions.hh>

namespace utility {
namespace json_spirit {

bool has_value(const mObject& obj, const std::string& name )
{
	auto i = obj.find(name);

	if ( i != obj.end() ) return true;
	return false;
}

mValue get_value(const mObject& obj, const std::string& name )
{
	auto i = obj.find(name);
	if ( i == obj.end() ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,"Cannot find member '" + name + "'" );
	}
	return i->second;
}

mObject get_mObject(const mObject& obj, const std::string& name )
{
	mValue value = get_value( obj, name );
	if ( value.type() != obj_type ) {  throw CREATE_EXCEPTION(utility::excn::Exception,"JSON error: '" + name + "' is not an object" ); }
	return value.get_obj();
}

mArray get_mArray(const mObject& obj, const std::string& name )
{
	const mValue value = get_value( obj, name );
	if ( value.type() != array_type ) {  throw CREATE_EXCEPTION(utility::excn::Exception,"JSON error: '" + name + "' is not an Array" ); }
	return value.get_array();
}

mArray read_mArray(  const std::string& json_string ){
	mValue value;
	read_or_throw(json_string, value);
	if ( value.type() != array_type ) {  throw CREATE_EXCEPTION(utility::excn::Exception,"JSON read error: expected an Array"); }
	return value.get_array();
}


mObject read_mObject(  const std::string& json_string ){
	mValue value;
	read_or_throw(json_string, value);
	if ( value.type() != obj_type ) {  throw CREATE_EXCEPTION(utility::excn::Exception,"JSON read error: expected an Object"); }
	return value.get_obj();
}

std::string get_string(const mObject& obj, const std::string& name )
{
	mValue value = get_value( obj, name );
	if ( value.type() != str_type ) {  throw CREATE_EXCEPTION(utility::excn::Exception,"JSON error: '" + name + "' is not a string" ); }
	return value.get_str();
}

std::string get_string_or_empty(const mObject& obj, const std::string& name )
{
	const std::string empty_string("");
	if ( !has_value( obj, name ) ) return empty_string;
	mValue value = get_value( obj, name );
	if ( value.type() != str_type ) { return empty_string; }
	return value.get_str();
}

double get_real(const mObject& obj, const std::string& name )
{
	mValue value = get_value( obj, name );
	if ( value.type() == real_type ) return value.get_real();
	if ( value.type() == int_type ) return get_int( obj, name );
	throw CREATE_EXCEPTION(utility::excn::Exception,"JSON error: '" + name + "' is not a number" );
	return 0.0;
}

double get_real_or_zero(const mObject& obj, const std::string& name )
{
	if ( !has_value( obj, name ) ) return 0.0;
	return get_real( obj, name );
}

int get_int(const mObject& obj, const std::string& name )
{
	mValue value = get_value( obj, name );
	if ( value.type() != int_type ) {  throw CREATE_EXCEPTION(utility::excn::Exception,"JSON error: '" + name + "' is not a int" ); }
	return value.get_int();
}

int get_int_or_zero(const mObject& obj, const std::string& name )
{
	if ( !has_value( obj, name ) ) return 0;
	return get_real( obj, name );
}

}
}
