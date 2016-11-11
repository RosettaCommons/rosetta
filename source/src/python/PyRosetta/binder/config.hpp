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


#ifndef _INCLUDED_config_hpp_
#define _INCLUDED_config_hpp_

#include <string>
#include <vector>
#include <map>

namespace binder {

class Config
{
	typedef std::string string;

	Config() {}

	Config(string const &root_module_, std::vector<string> namespaces_to_bind_, std::vector<string> namespaces_to_skip_, string const &prefix_, std::size_t maximum_file_length_) :
		root_module(root_module_), namespaces_to_bind(namespaces_to_bind_), namespaces_to_skip(namespaces_to_skip_), prefix(prefix_), maximum_file_length(maximum_file_length_) {}

private:
	std::map<string, string> binders_, add_on_binders_;

	string default_static_pointer_return_value_policy_          = "pybind11::return_value_policy::automatic";
	string default_static_lvalue_reference_return_value_policy_ = "pybind11::return_value_policy::automatic";
	string default_static_rvalue_reference_return_value_policy_ = "pybind11::return_value_policy::automatic";

	string default_member_pointer_return_value_policy_          = "pybind11::return_value_policy::automatic";
	string default_member_lvalue_reference_return_value_policy_ = "pybind11::return_value_policy::automatic";
	string default_member_rvalue_reference_return_value_policy_ = "pybind11::return_value_policy::automatic";

public:
	static Config &get();

	/// Read config setting from the file
	void read(string const &file_name);

	string root_module;

	std::vector<string> namespaces_to_bind, classes_to_bind, functions_to_bind,
						namespaces_to_skip, classes_to_skip, functions_to_skip,
						includes_to_add, includes_to_skip;

	std::map<string, string> const &binders() const { return binders_; }
	std::map<string, string> const &add_on_binders() const { return add_on_binders_; }

	string const &default_static_pointer_return_value_policy()          { return default_static_pointer_return_value_policy_; }
	string const &default_static_lvalue_reference_return_value_policy() { return default_static_lvalue_reference_return_value_policy_; }
	string const &default_static_rvalue_reference_return_value_policy() { return default_static_rvalue_reference_return_value_policy_; }

	string const &default_member_pointer_return_value_policy()          { return default_member_pointer_return_value_policy_; }
	string const &default_member_lvalue_reference_return_value_policy() { return default_member_lvalue_reference_return_value_policy_; }
	string const &default_member_rvalue_reference_return_value_policy() { return default_member_rvalue_reference_return_value_policy_; }

	string prefix;

	std::size_t maximum_file_length;

	/// check if user requested binding for given declaration
	bool is_namespace_binding_requested(string const &namespace_) const;
	bool is_namespace_skipping_requested(string const &namespace_) const;

	bool is_function_binding_requested(string const &function) const;
	bool is_function_skipping_requested(string const &function) const;

	bool is_class_binding_requested(string const &class_) const;
	bool is_class_skipping_requested(string const &class_) const;

	bool is_include_skipping_requested(string const &include) const;

	string includes_code() const;
};


} // namespace binder

#endif // _INCLUDED_config_hpp_
