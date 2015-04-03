// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/JobOptions.cc
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <basic/resource_manager/JobOptions.hh>

// Platform Headers
#include <platform/types.hh>
#include <utility/vector1.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.hh>


//C++ Headers
#include <string>
#include <map>
#include <iomanip>

namespace basic {
namespace resource_manager {

/// @details Auto-generated virtual destructor
JobOptions::~JobOptions() {}

using platform::Real;
using utility::options::BooleanOptionKey;
using utility::options::BooleanVectorOptionKey;
using utility::options::FileOptionKey;
using utility::options::FileVectorOptionKey;
using utility::options::IntegerOptionKey;
using utility::options::IntegerVectorOptionKey;
using utility::options::PathOptionKey;
using utility::options::PathVectorOptionKey;
using utility::options::RealOptionKey;
using utility::options::RealVectorOptionKey;
using utility::options::StringOptionKey;
using utility::options::StringVectorOptionKey;
using utility::excn::EXCN_Msg_Exception;
using utility::file::FileName;
using utility::file::PathName;
using utility::vector1;
using std::map;
using std::string;
using std::endl;
using std::setw;

void
JobOptions::show(
	std::ostream & out
) const {
	using namespace utility::options;
	out
		<< std::setiosflags(std::ios::left) << setw(16) << "OptionKeyType"
		<< std::setiosflags(std::ios::left) << setw(16) << "OptionKey"
		<< "OptionValue" << endl;
	for(
		map<BooleanOptionKey, bool>::const_iterator
			o = boolean_options_.begin(), oe = boolean_options_.end(); o != oe; ++o){
		out
			<< std::setiosflags(std::ios::left) << setw(16) << "Boolean"
			<< std::setiosflags(std::ios::left) << setw(16) << o->first.id()
			<< o->second << endl;
	}

	for(
		map<BooleanVectorOptionKey, vector1<bool> >::const_iterator
			o = boolean_vector_options_.begin(), oe = boolean_vector_options_.end();
		o != oe; ++o){
		out
			<< std::setiosflags(std::ios::left) << setw(16) << "BooleanVector"
			<< std::setiosflags(std::ios::left) << setw(16) << o->first.id();
		bool first(true);
		for(
			vector1<bool>::const_iterator
				k = o->second.begin(), ke = o->second.end(); k != ke; ++k){
			if(!first){
				out << " ";
			} else {
				first = false;
			}

			out	<< (*k ? "true" : "false");
		}
		out << endl;
	}
}

std::ostream &
operator<<(
	std::ostream & out,
	const JobOptions & job_options
) {
	job_options.show(out);
	return out;
}

void
JobOptions::add_option(
	BooleanOptionKey key,
	bool val
) {
	boolean_options_[key] = val;
}

bool
JobOptions::has_option(
	BooleanOptionKey key
) const {
	return boolean_options_.find(key) != boolean_options_.end();
}

bool
JobOptions::get_option(
	BooleanOptionKey key
) const {
	map< BooleanOptionKey, bool >::const_iterator val(
		boolean_options_.find(key));
	if(val == boolean_options_.end()){
		throw EXCN_Msg_Exception("BooleanOptionKey '" + key.identifier() + "' not found in JobOptions");
	} else {
		return val->second;
	}
}


void
JobOptions::add_option(
	BooleanVectorOptionKey key,
	vector1< bool > const & val
) {
	boolean_vector_options_[key] = val;
}

bool
JobOptions::has_option(
	BooleanVectorOptionKey key
) const {
	return boolean_vector_options_.find(key) != boolean_vector_options_.end();
}

vector1< bool > const &
JobOptions::get_option(
	BooleanVectorOptionKey key
) const {
	map< BooleanVectorOptionKey, vector1< bool > >::const_iterator val(
		boolean_vector_options_.find(key));
	if(val == boolean_vector_options_.end()){
		throw EXCN_Msg_Exception("BooleanVectorOptionKey '" + key.identifier() + "' not found in JobOptions");
	} else {
		return val->second;
	}
}


void
JobOptions::add_option(
	FileOptionKey key,
	FileName const & val
) {
	file_options_[key] = val;
}

bool
JobOptions::has_option(
	FileOptionKey key
) const {
	return file_options_.find(key) != file_options_.end();
}

FileName const &
JobOptions::get_option(
	FileOptionKey key
) const {
	map< FileOptionKey, FileName >::const_iterator val(
		file_options_.find(key));
	if(val == file_options_.end()){
		throw EXCN_Msg_Exception("FileOptionKey '" + key.identifier() + "' not found in JobOptions");
	} else {
		return val->second;
	}
}


void
JobOptions::add_option(
	FileVectorOptionKey key,
	vector1< FileName > const & val
) {
	file_vector_options_[key] = val;
}

bool
JobOptions::has_option(
	FileVectorOptionKey key
) const {
	return file_vector_options_.find(key) != file_vector_options_.end();
}


vector1< FileName > const &
JobOptions::get_option(
	FileVectorOptionKey key
) const {
	map< FileVectorOptionKey, vector1< FileName > >::const_iterator val(
		file_vector_options_.find(key));
	if(val == file_vector_options_.end()){
		throw EXCN_Msg_Exception("BooleanOptionKey '" + key.identifier() + "' not found in JobOptions");
	} else {
		return val->second;
	}
}


void
JobOptions::add_option(
	IntegerOptionKey key,
	int val
) {
	integer_options_[key] = val;
}

bool
JobOptions::has_option(
	IntegerOptionKey key
) const {
	return integer_options_.find(key) != integer_options_.end();
}

int
JobOptions::get_option(
	IntegerOptionKey key
) const {
	map< IntegerOptionKey, int >::const_iterator val(
		integer_options_.find(key));
	if(val == integer_options_.end()){
		throw EXCN_Msg_Exception("IntegerOptionKey '" + key.identifier() + "' not found in JobOptions");
	} else {
		return val->second;
	}
}


void
JobOptions::add_option(
	IntegerVectorOptionKey key,
	vector1< int > const & val
) {
	integer_vector_options_[key] = val;
}

bool
JobOptions::has_option(
	IntegerVectorOptionKey key
) const {
	return integer_vector_options_.find(key) != integer_vector_options_.end();
}

vector1< int > const &
JobOptions::get_option(
	IntegerVectorOptionKey key
) const {
	map< IntegerVectorOptionKey, vector1< int > >::const_iterator val(
		integer_vector_options_.find(key));
	if(val == integer_vector_options_.end()){
		throw EXCN_Msg_Exception("IntegerVectorOptionKey '" + key.identifier() + "' not found in JobOptions");
	} else {
		return val->second;
	}
}


void
JobOptions::add_option(
	PathOptionKey key,
	PathName const & val
) {
	path_options_[key] = val;
}

bool
JobOptions::has_option(
	PathOptionKey key
) const {
	return path_options_.find(key) != path_options_.end();
}

PathName const &
JobOptions::get_option(
	PathOptionKey key
) const {
	map< PathOptionKey, PathName >::const_iterator val(
		path_options_.find(key));
	if(val == path_options_.end()){
		throw EXCN_Msg_Exception("PathOptionKey '" + key.identifier() + "' not found in JobOptions");
	} else {
		return val->second;
	}
}


void
JobOptions::add_option(
	PathVectorOptionKey key,
	vector1< PathName > const & val
) {
	path_vector_options_[key] = val;
}

bool
JobOptions::has_option(
	PathVectorOptionKey key
) const {
	return path_vector_options_.find(key) != path_vector_options_.end();
}

vector1< PathName > const &
JobOptions::get_option(
	PathVectorOptionKey key
) const {
	map< PathVectorOptionKey, vector1< PathName > >::const_iterator val(
		path_vector_options_.find(key));
	if(val == path_vector_options_.end()){
		throw EXCN_Msg_Exception("PathVectorOptionKey '" + key.identifier() + "' not found in JobOptions");
	} else {
		return val->second;
	}
}


void
JobOptions::add_option(
	RealOptionKey key,
	Real val
) {
	real_options_[key] = val;
}

bool
JobOptions::has_option(
	RealOptionKey key
) const {
	return real_options_.find(key) != real_options_.end();
}

Real
JobOptions::get_option(
	RealOptionKey key
) const {
	map< RealOptionKey, Real >::const_iterator val(
		real_options_.find(key));
	if(val == real_options_.end()){
		throw EXCN_Msg_Exception("RealOptionKey '" + key.identifier() + "' not found in JobOptions");
	} else {
		return val->second;
	}
}


void
JobOptions::add_option(
	RealVectorOptionKey key,
	vector1< Real > const & val
) {
	real_vector_options_[key] = val;
}

bool
JobOptions::has_option(
	RealVectorOptionKey key
) const {
	return real_vector_options_.find(key) != real_vector_options_.end();
}

vector1< Real > const &
JobOptions::get_option(
	RealVectorOptionKey key
) const {
	map< RealVectorOptionKey, vector1< Real> >::const_iterator val(
		real_vector_options_.find(key));
	if(val == real_vector_options_.end()){
		throw EXCN_Msg_Exception("RealVectorOptionKey '" + key.identifier() + "' not found in JobOptions");
	} else {
		return val->second;
	}
}


void
JobOptions::add_option(
	StringOptionKey key,
	string val
) {
	string_options_[key] = val;
}

bool
JobOptions::has_option(
	StringOptionKey key
) const {
	return string_options_.find(key) != string_options_.end();
}

string const &
JobOptions::get_option(
	StringOptionKey key
) const {
	map< StringOptionKey, string >::const_iterator val(
		string_options_.find(key));
	if(val == string_options_.end()){
		throw EXCN_Msg_Exception("StringOptionKey '" + key.identifier() + "' not found in JobOptions");
	} else {
		return val->second;
	}
}


void
JobOptions::add_option(
	StringVectorOptionKey key,
	vector1< string > const & val
) {
	string_vector_options_[key] = val;
}

bool
JobOptions::has_option(
	StringVectorOptionKey key
) const {
	return string_vector_options_.find(key) != string_vector_options_.end();
}

vector1< string > const &
JobOptions::get_option(
	StringVectorOptionKey key
) const {
	map< StringVectorOptionKey, vector1< string > >::const_iterator val(
		string_vector_options_.find(key));
	if(val == string_vector_options_.end()){
		throw EXCN_Msg_Exception("StringVectorOptionKey '" + key.identifier() + "' not found in JobOptions");
	} else {
		return val->second;
	}
}

} // namespace
} // namespace
