// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   binder/config.hpp
/// @brief  Support for Binder Config file
/// @author Sergey Lyskov

#include <config.hpp>

#include <util.hpp>


#include <llvm/Support/raw_ostream.h>

#include <stdexcept>
#include <fstream>

using namespace llvm;

using std::string;

namespace binder {

/// Split string in two by ether space or tab character
std::pair<string, string> split_in_two(string const &s, string const & error_string)
{
	size_t space = s.find(' ');
	size_t tab = s.find('\t');

	size_t split = std::min(space, tab);

	if( split == string::npos ) throw std::runtime_error(error_string);
	else return std::make_pair( s.substr(0, split), s.substr(split+1) );
}


Config &Config::get()
{
	static Config * instance = nullptr;

	if( !instance ) instance = new Config();

	return *instance;
}

/// Read config setting from file
void Config::read(string const &file_name)
{
	string const _namespace_     {"namespace"};
	string const _function_      {"function"};
	string const _class_         {"class"};
	string const _include_       {"include"};
	string const _binder_        {"binder"};
	string const _add_on_binder_ {"add_on_binder"};

	string const _default_pointer_return_value_policy_          {"default_pointer_return_value_policy"};
	string const _default_lvalue_reference_return_value_policy_	{"default_lvalue_reference_return_value_policy"};
	string const _default_rvalue_reference_return_value_policy_ {"default_rvalue_reference_return_value_policy"};

	std::ifstream f(file_name);
	string line;

	while( std::getline(f, line) ) {
		if( line.size() ) {
			if( line[0] == '#' ) continue;

			if( line[0] == '+'  or  line[0] == '-' ) {
				size_t space = line.find(' ');
				if(  space == string::npos ) throw std::runtime_error("Invalid line in config file! Each line must have token separated with space from object name. For example: '+function aaa::bb::my_function'. Line: " + line);
				else {
					bool bind = line[0] == '+' ? true : false;
					string token = line.substr(1, space-1);
					string name = line.substr(space+1);
					string name_without_spaces = name;
					name_without_spaces.erase(std::remove(name_without_spaces.begin(), name_without_spaces.end(), ' '), name_without_spaces.end());

					//outs() << token << " " << name << "\n";

					if( token == _namespace_ ) {

						if(bind) namespaces_to_bind.push_back(name_without_spaces); else namespaces_to_skip.push_back(name_without_spaces);

					} else if( token == _class_ ) {

						if(bind) classes_to_bind.push_back(name_without_spaces); else classes_to_skip.push_back(name_without_spaces);

					} else if( token == _function_ ) {

						if(bind) functions_to_bind.push_back(name_without_spaces); else functions_to_skip.push_back(name_without_spaces);

					} else if( token == _include_ ) {

						if(bind) includes_to_add.push_back(name_without_spaces);
						else includes_to_skip.push_back(name_without_spaces);

					} else if( token == _binder_ ) {

						if(bind) {
							auto binder_function = split_in_two(name, "Invalid line for binder specification! Must be: name_of_type + <space or tab> + name_of_binder. Got: " + line);
							binders_[binder_function.first] = binder_function.second;
						}

					} else if( token == _add_on_binder_ ) {

						if(bind) {
							auto binder_function = split_in_two(name, "Invalid line for add_on_binder specification! Must be: name_of_type + <space or tab> + name_of_binder. Got: " + line);
							add_on_binders_[binder_function.first] = binder_function.second;
						}
					}
					else if( token == _default_pointer_return_value_policy_ )          default_pointer_return_value_policy_ = name_without_spaces;
					else if( token == _default_lvalue_reference_return_value_policy_ ) default_lvalue_reference_return_value_policy_ = name_without_spaces;
					else if( token == _default_rvalue_reference_return_value_policy_ ) default_rvalue_reference_return_value_policy_ = name_without_spaces;
					else throw std::runtime_error("Invalid token in config file! Each token must be ether: namespace, class or function! For example: '+function aaa::bb::my_function'. Token: '" + token + "' Line: '" + line + '\'');
				}
			}
			else throw std::runtime_error("Invalid token at the begining of line in config file! Each line should begin with ether '+' or '-' or '#'! Line: " + line);
		}
	}
}


bool Config::is_namespace_binding_requested(string const &namespace_) const
{
	bool to_bind_flag=false, to_skip_flag=false;
	string to_bind, to_skip;

	for(auto &n : namespaces_to_bind) {
		if( n.size() == 0 ) to_bind_flag=true;
		if( begins_with(namespace_, n) ) {
			if( n.size() > to_bind.size() ) { to_bind = n; to_bind_flag=true; }
		}
	}

	for(auto &s : namespaces_to_skip) {
		if( begins_with(namespace_, s) ) {
			if( s.size() > to_skip.size() ) { to_skip = s; to_skip_flag=true; }
		}
	}

	if( to_bind.size() > to_skip.size() ) return true;
	if( to_bind.size() < to_skip.size() ) return false;

	if( to_bind_flag and to_skip_flag ) throw std::runtime_error("Could not determent if namespace '" + namespace_ + "' should be binded or not... please check if options --bind and --skip conflicting!!!");

	if( to_bind_flag ) return true;

	return false;
}

bool Config::is_namespace_skipping_requested(string const &namespace_) const
{
	bool to_bind_flag=false, to_skip_flag=false;
	string to_bind, to_skip;

	for(auto &n : namespaces_to_bind) {
		if( begins_with(namespace_, n) ) {
			if( n.size() > to_bind.size() ) { to_bind = n; to_bind_flag=true; }
		}
	}

	for(auto &s : namespaces_to_skip) {
		if( begins_with(namespace_, s) ) {
			if( s.size() > to_skip.size() ) { to_skip = s; to_skip_flag=true; }
		}
	}

	if( to_bind.size() > to_skip.size() ) return false;
	if( to_bind.size() < to_skip.size() ) return true;

	if( to_bind_flag and to_skip_flag ) throw std::runtime_error("Could not determent if namespace '" + namespace_ + "' should be binded or not... please check if options --bind and --skip conflicting!!!");

	if( to_skip_flag ) return true;

	return false;
}


bool Config::is_function_binding_requested(string const &function_) const
{
	string function {function_};
	function.erase(std::remove(function.begin(), function.end(), ' '), function.end());

	auto bind = std::find(functions_to_bind.begin(), functions_to_bind.end(), function);

	if( bind != functions_to_bind.end() ) return true;

	return false;
}

bool Config::is_function_skipping_requested(string const &function_) const
{
	string function {function_};
	function.erase(std::remove(function.begin(), function.end(), ' '), function.end());

	auto bind = std::find(functions_to_skip.begin(), functions_to_skip.end(), function);

	if( bind != functions_to_skip.end() ) {
		//outs() << "Skipping: " << function << "\n";
		return true;
	}

	return false;
}



bool Config::is_class_binding_requested(string const &class__) const
{
	string class_ {class__};
	class_.erase(std::remove(class_.begin(), class_.end(), ' '), class_.end());

	auto bind = std::find(classes_to_bind.begin(), classes_to_bind.end(), class_);

	if( bind != classes_to_bind.end() ) return true;

	return false;
}


bool Config::is_class_skipping_requested(string const &class__) const
{
	string class_ {class__};
	class_.erase(std::remove(class_.begin(), class_.end(), ' '), class_.end());

	auto bind = std::find(classes_to_skip.begin(), classes_to_skip.end(), class_);

	if( bind != classes_to_skip.end() ) {
		//outs() << "Skipping: " << class_ << "\n";
		return true;
	}

	return false;
}


bool Config::is_include_skipping_requested(string const &include) const
{
	// auto skip = std::find(includes_to_skip.begin(), includes_to_skip.end(), include);
	// if( skip != includes_to_skip.end() ) return true;
	// else return false;
	for(auto & i : includes_to_skip)
		if( begins_with(include, i) ) return true;

	return false;
}


string Config::includes_code() const
{
	string c;
	for(auto & i: includes_to_add) c += "#include " + i + "\n";
	return c.size() ? c+'\n' : c;
}

} // namespace binder
