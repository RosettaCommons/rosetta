// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/Emitter.cc
///
/// @brief  Lightweight class to ease writting YAML documents
/// @author Ian W. Davis

#include <basic/Emitter.hh>

// Boost Headers
#include <boost/foreach.hpp>

#include <platform/types.hh>
#include <utility/down_cast.hh>
#include <utility/vector1_bool.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <ostream>
#include <sstream>
#include <basic/Tracer.fwd.hh>
#include <basic/Tracer.hh>

namespace basic {


void Emitter::start_list(bool indent/*=true*/)
{
	assert_in(false, "Tried to write list data inside a map");
	do_indent();
	start_raw(false, indent);
}

void Emitter::start_map(bool indent/*=true*/)
{
	assert_in(false, "Tried to write list data inside a map");
	do_indent();
	start_raw(true, indent);
}

void Emitter::start_list(const char * label, bool indent/*=true*/)
{
	assert_in(true, "Tried to write map data inside a list");
	do_indent();
	write_label(label);
	start_raw(false, indent);
}

void Emitter::start_map(const char * label, bool indent/*=true*/)
{
	assert_in(true, "Tried to write map data inside a list");
	do_indent();
	write_label(label);
	start_raw(true, indent);
}

void Emitter::end_list()
{
	if(assert_in(false, "Tried to end list inside map")) end_raw();
}

void Emitter::end_map()
{
	if(assert_in(true, "Tried to end map inside list")) end_raw();
}

void Emitter::end(size_t desired_depth/*=0*/)
{
	//while( depth() > desired_depth ) end_raw();
	for(size_t i = depth(); i > desired_depth; --i) end_raw();
}

/// @details Anything but the most basic characters needs to be quoted and escaped.
/// For normal YAML, very simple text can be output without quotes, though.
/// @param needs_quotes_out will be set to true if string contains "special" characters.
std::string Emitter::escape_string(std::string const & s, bool & needs_quotes_out)
{
	using std::string;
	// Characters that need no escaping, and no quoting (conservative)
	// Most punctuation characters have special meaning in YAML!
	static string bare_chars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_. ";
	static string hex = "0123456789ABCDEF";
	// If all chars in s are also in bare_chars, no quotes are needed in YAML
	needs_quotes_out = !(s.find_first_not_of(bare_chars) == string::npos);
	// Special cases: starts/ends with whitespace, zero length
	if( s.empty() || s[0] == ' ' || s[s.size()-1] == ' ' ) needs_quotes_out = true;
	if( !needs_quotes_out ) return s;

	std::ostringstream o;

	BOOST_FOREACH(char ch, s){
		if( ' ' <= ch && ch <= '~' ) o << ch; // the range of printable ASCII characters
		else if( ch == '"'  ) o << "\\\"";
		else if( ch == '\\' ) o << "\\\\";
		else if( ch == '\n' ) o << "\\n";
		else if( ch == '\r' ) o << "\\r";
		else if( ch == '\t' ) o << "\\t";
		else if( ch == '\f' ) o << "\\f";
		else if( ch == '\b' ) o << "\\b";
		else { // Unicode escape
			o << "\\u";
			o << hex[ ((ch>>24)&0x000F) ];
			o << hex[ ((ch>>16)&0x000F) ];
			o << hex[ ((ch>> 8)&0x000F) ];
			o << hex[ ((ch    )&0x000F) ];
		}
	}
	return o.str();
}



} // basic

