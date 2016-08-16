// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/UTools.hh
/// @brief  Unit test Tracer/Diff system
/// @author Sergey Lyskov


#ifndef INCLUDED_UTools_HH
#define INCLUDED_UTools_HH

#include <string>
#include <vector>
#include <ostream>
#include <sstream>
#include <cmath>
#include <cstdlib> //required by GCC 4.3.2

namespace test {
namespace utools {

enum ItemType {
	_IT_None_,  // type was not set yet.
	_IT_String_,
	_IT_Number_,
	_IT_MarkupAbsTolerance_,
	_IT_MarkupRelTolerance_
};

struct Item {
	Item(ItemType t=_IT_None_, double n=0.0, std::string const & s="") : type(t), num(n), str(s) {}

	ItemType type;
	double num; /// number and markup value stored here
	std::string str;
};


inline std::ostream& operator <<(std::ostream &os, std::vector<Item> const & v)
{
	os << "[";
	for(unsigned int i=0; i<v.size(); i++) {
		os << "(" << v[i].type << " ";

		if( v[i].type == _IT_String_ ) os << v[i].str;
		if( v[i].type == _IT_Number_ ) os << v[i].num;
		if( v[i].type == _IT_MarkupAbsTolerance_ ) os << v[i].num;
		if( v[i].type == _IT_MarkupRelTolerance_ ) os << v[i].num;

		os << "), ";
	}
	os << "]";
	return os;
}


inline bool isFloatNumber(std::string const &s, unsigned int start_pos, double &res, unsigned int &end_pos, bool isMarkup=false)
{
	if( start_pos >= s.size() ) return false; // if( s.size() == 0 ) return false;

	if( isMarkup ) {
		if( s[start_pos] == ' ' ) return false; // we don't want it to skip spaces
	}
	else {
		if( s[start_pos] != ' ' ) return false; // we always want at least one space before number
		if( start_pos + 1 < s.size() ) {
			if( s[start_pos+1] == ' ' ) return false; // only one space before number is allowed
		}
	}

	char const * end = &s[start_pos];
	res = strtod(&s[start_pos], (char **)&end);

	end_pos = end - (&s[0]);

	if( !isMarkup ) {
		if( end_pos < s.size() ) {  // check if number end up with space (only if its not the end of a line)
			if( s[end_pos] != ' ' ) return false;
		}
	}

	if( end_pos > start_pos ) return true;
	else return false;
}


inline bool isMarkup(std::string const &s, unsigned int start_pos, double &res, unsigned int &end_pos, std::string markup)
{
	//std::string markup="set_abs_tolerance(";
	if( start_pos >= s.size() ) return false; // if( s.size() == 0 ) return false;
	unsigned int r = s.find(markup, start_pos);
	if( r != start_pos ) return false;

	unsigned int end;
	if( !isFloatNumber(s, r+markup.size(), res, end, true) ) return false;
	if( end >= s.size() ) return false;
	if( s[end] != ')' ) return false;
	end_pos = end + 1;
	return true;
}

inline bool smartAdd(bool f, double const &res, Item &I, ItemType type, std::vector<Item> &V)
{
	if(f) {
		if( I.type == _IT_String_ ) {
			//I.num = 0;
			V.push_back(I);
		}
		I = Item(type, res);
		//I.type = type;
		//I.num = res;
		V.push_back(I);
		//I = Item();
		I.type = _IT_None_;
	}
	return f;
}

/// Parse string in to vector of items
inline std::vector<Item> ParseString(std::string const &s)
{
	std::vector<Item> V;
	Item r;

	for(unsigned int i=0; i<s.size(); ) {
		double res;
		unsigned int end_pos(0); // initialized to zero to keep compiler happy

		if( smartAdd(isFloatNumber(s, i, res, end_pos), res, r, _IT_Number_, V) )
			{ i=end_pos;  continue; }
		if( smartAdd(isMarkup(s, i, res, end_pos, "set_abs_tolerance("), res, r, _IT_MarkupAbsTolerance_, V) )
			{ i=end_pos;  continue; }
		if( smartAdd(isMarkup(s, i, res, end_pos, "set_rel_tolerance("), res, r, _IT_MarkupRelTolerance_, V) )
			{  i=end_pos;  continue; }

		if( r.type == _IT_String_ ) {
			//r.num = 0;
			r.str.push_back( s[i] );
		}
		else {
			r.type = _IT_String_;
			r.str = "";
			r.str.push_back( s[i] );
		}
		i++;
		continue;
	}
	if( r.type == _IT_String_ ) {
		V.push_back(r);
	}
	return V;
}

/// Smart compare strings, return true if they eq (with tolerance).
inline bool isEq(std::string const & s1, std::string const & s2,
		  double &abs_tolerance, double &rel_tolerance, std::string &error_message)
{
	error_message = "";

	std::vector<Item> v1 = ParseString(" " + s1 + " ");
	std::vector<Item> v2 = ParseString(" " + s2 + " ");

	// Removing markup from vector 2.
	for(int i=v2.size()-1; i>=0; i--) {
		if( v2[i].type == _IT_None_ ) {
			error_message = "Wrong Item type: _IT_None_! Probably a bug in UTools.hh...";
			return false;
		}

		if( v2[i].type == _IT_MarkupAbsTolerance_ || v2[i].type == _IT_MarkupRelTolerance_ ) {
			v2.erase( v2.begin() + i );
			continue;
		}
	}


	unsigned int i=0, j=0;
	for(; i<v1.size(); i++) {
		Item &I1(v1[i]);
		Item &I2(v2[j]);

		if( I1.type == _IT_None_ ) {
			error_message = "Wrong Item type: _IT_None_! Probably a bug in UTools.hh...";
			return false;
		}

		if( I1.type == _IT_MarkupAbsTolerance_ ) {
			abs_tolerance = I1.num;
			continue;
		}

		if( I1.type == _IT_MarkupRelTolerance_ ) {
			rel_tolerance = I1.num;
			continue;
		}

		if( j>= v2.size() ) break;

		if( I1.type != I2.type ) return false;

		if( I1.type == _IT_String_ ) {
			if( I1.str != I2.str ) {
				error_message = I1.str + " != " + I2.str;
				return false;
			}
			j++;  continue;
		}

		if( I1.type == _IT_Number_ ) {
			double p = (I1.num + I2.num)/2. * rel_tolerance;
			p = std::max(abs_tolerance, p);

			if( fabs(I1.num - I2.num) > p ) {
				std::ostringstream s;
				s <<  I1.num << "!=" << I2.num << " [abs_tolerance=" << abs_tolerance << ", rel_tolerance=" << rel_tolerance << "]";
				error_message = s.str();
				return false;
			}

			j++;  continue;
		}

		error_message = "Unknow Item type! Probably a bug in UTools.hh...";
		return false;
	}

	if( i<v1.size() || j<v2.size() ) return false; // Some elements left in vectors - false!

	return true;
}

} // namespace utools
} // namespace test


#endif // INCLUDED_test_UTools_hh
