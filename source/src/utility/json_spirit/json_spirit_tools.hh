// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/json_spirit/json_spirit_tools.hh
/// @brief  Extra utilitiy functions not provided by json spirit library
/// @author Mike Tyka

#ifndef JSON_SPIRIT_TOOLS
#define JSON_SPIRIT_TOOLS

#include <utility/json_spirit/json_spirit_value.h>

namespace utility {
namespace json_spirit {


bool has_value(const mObject& obj, const std::string& name );

mValue get_value(const mObject& obj, const std::string& name );

mObject get_mObject(const mObject& obj, const std::string& name );

mArray get_mArray(const mObject& obj, const std::string& name );

mObject read_mObject(  const std::string& json_string );

mArray read_mArray(  const std::string& json_string );

std::string get_string(const mObject& obj, const std::string& name );

std::string get_string_or_empty(const mObject& obj, const std::string& name );

double get_real(const mObject& obj, const std::string& name );

int get_int(const mObject& obj, const std::string& name );

}
}

#endif


