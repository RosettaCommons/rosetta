// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   binder/context.cpp
/// @brief  Data structures to represent root context and modules
/// @author Sergey Lyskov


#include <context.hpp>
#include <binder.hpp>

#include <type.hpp>
#include <class.hpp>
#include <util.hpp>
#include <fmt/format.h>

#include <clang/AST/ASTContext.h>

#include <sstream>
#include <fstream>
#include <set>

using namespace llvm;
using namespace clang;
using std::string;
using std::vector;
using std::unordered_map;
using std::make_pair;

using namespace fmt::literals;
using fmt::format;

namespace binder {

//const std::string _module_variable_name_{"M"};



// check if declaration is already in stack with level at lease as 'level' or lower and add it if it is not - return true if declaration was added
bool IncludeSet::add_decl(clang::NamedDecl const *D, int level)
{
	// auto r = stack_.insert( make_pair(D, level) );
	// if( !r.second  and  (*r.first).second <= level ) return false;

	if( stack_.count(D)  and  stack_[D] <= level ) return false;

	stack_[D] = level;
	return true;
}


// remove all includes and clear up the stack
void IncludeSet::clear()
{
	includes_.clear();
	stack_.clear();
}


/// return true if object declared in system header
bool Binder::is_in_system_header()
{
	NamedDecl * decl( named_decl() );
	ASTContext & ast_context( decl->getASTContext() );
	SourceManager & sm( ast_context.getSourceManager() );

	return FullSourceLoc( decl->getLocation(), sm ).isInSystemHeader();
}


// return true if code was already generate for this object
bool Binder::is_binded() const
{
	return code().size()  or is_python_builtin( named_decl() );
}



llvm::raw_ostream & operator << (llvm::raw_ostream & os, Binder const &b)
{
	clang::NamedDecl *decl = b.named_decl();

	string name  = decl->getNameAsString();
	string qualified_name = decl->getQualifiedNameAsString();
	string path = decl->getQualifiedNameAsString().substr(0, qualified_name.size() - name.size() );

	return os << "B{name=" << name << ", path=" << path << "\n";  //<< ", include= " code=\n" << b("module") << "\n}
}

void Context::add(BinderOP &b)
{
	if( ids.count( b->id() ) ) {
		//errs() << "Skipping adding duplicate binder for: " + b->id() + "...\n";
		return;
	}

	binders.push_back(b);
	ids.insert( b->id() );

	if( TypeDecl * type_decl = dyn_cast<TypeDecl>( b->named_decl() ) ) types[ typename_from_type_decl(type_decl) ] = b;
}


/// check if forward declaration for CXXRecordDecl needed
bool Context::is_forward_needed(clang::CXXRecordDecl const *C)
{
	return !binded.count( class_qualified_name(C) )  and  !is_python_builtin(C);
}


/// add given class to 'aleady binded' set
void Context::add_to_binded(CXXRecordDecl const *C)
{
	binded.insert( class_qualified_name(C) );
}


/// examine binded objects and recursivly create all nested namespaces
std::set<string> Context::create_all_nested_namespaces()
{
	vector<string> namespaces;//{""};

	for(auto & b : binders) {
		if( b->code().size() ) {
			string ns = namespace_from_named_decl( b->named_decl() );

			while( ns.size() ) {
				namespaces.push_back(ns);
				ns = base_namespace(ns);
			}
		}
	}

	std::set<string> s( namespaces.begin(), namespaces.end() );

	return s;
}


// generate code for include directives and cleanup the includes vector
string generate_include_directives(IncludeSet const &include_set)
{
	string r;
	for(auto &i : std::set<string>(include_set.includes().begin(), include_set.includes().end() ) )
		if( !Config::get().is_include_skipping_requested(i) ) r += "#include " + i + '\n';

	//includes.resize(0);
	return r;
}


std::string Context::module_variable_name(std::string const& namespace_)
{
	return "M(\"" + namespace_ + "\")";
}


// find binder related to given object name and bind it
void Context::request_bindings(std::string const &type)
{
	if( types.count(type)  and  !types[type]->is_binded()  and  types[type]->bindable()  and  !types[type]->skipping_requested()  and  !is_python_builtin( types[type]->named_decl() ) ) {
		//if( begins_with(type, "std::__1::basic_ostream") ) outs() << "Context::bind: requesting bindings for: " << type << "...\n";
		//if( begins_with(type, "std::basic_ostream") ) outs() << "Context::bind: requesting bindings for: " << type << "...\n";
		types[type]->request_bindings();
	}
	//else outs() << "Context::bind: Could not find generator for type:" << type << "\n";
}


/// walk over all binders and bind one that belong to requested namespaces
void Context::bind(Config const &config)
{
	for(auto & sp : binders) {
		Binder & b( *sp );
		if( !b.is_in_system_header() and  b.bindable() ) b.request_bindings_and_skipping(config);
	}

	bool flag = true;
	int pass = 1;
	while(flag) {
		flag = false;

		outs() << "Generate bindings, pass №" << pass << "...\n";

		for(auto & sp : binders) {
			Binder & b( *sp );
			if( !b.is_binded()  and  b.bindable() and  b.binding_requested() ) {
				//outs() << "Binding: " << b.id() /*named_decl()->getQualifiedNameAsString()*/ << "\n";
				b.bind(*this);
				flag=true;
			}
		}

		++pass;
	}
}


/// Generate file name where binding for given generator should be stored
string file_name_prefix_for_binder(BinderOP &b)
{
	clang::NamedDecl *decl = b->named_decl();

	string include = relevant_include(decl);

	if( include.size() <= 2 ) { include = "<unknown/unknown.hh>";  outs() << "Warning: file_name_for_decl could not determent file name for decl: " + string(*b) + ", result is too short!\n"; } //throw std::runtime_error( "file_name_for_decl failed!!! include name for decl: " + string(*b) + " is too short!");
	include = include.substr(1, include.size()-2);

	//outs() << "qn:" << decl->getQualifiedNameAsString() << " namespace_from_named_decl(decl):" << namespace_from_named_decl(decl) << "\n";

	if( namespace_from_named_decl(decl) == "std"  or  begins_with(namespace_from_named_decl(decl), "std::" ) ) include = "std/" +  ( begins_with(include, "bits/") ? include.substr(5) : include );

	//return replace( replace( replace( replace(include, ".hh", ""), ".hpp", ""), ".h", ""), ".", "_");
	replace(include, ".hh", ""); replace(include, ".hpp", ""); replace(include, ".h", ""); replace(include, ".", "_");
	return include;
}


const char * main_module_header = R"_(#include <map>
#include <memory>
#include <stdexcept>
#include <functional>

#include <pybind11/pybind11.h>

typedef std::function< pybind11::module & (std::string const &) > ModuleGetter;

{0}

PYBIND11_PLUGIN({1}) {{
	std::map <std::string, std::shared_ptr<pybind11::module> > modules;
	ModuleGetter M = [&](std::string const &namespace_) -> pybind11::module & {{
		auto it = modules.find(namespace_);
		if( it == modules.end() ) throw std::runtime_error("Attempt to access pybind11::module for namespace " + namespace_ + " before it was created!!!");
		return * it->second;
	}};

	modules[""] = std::make_shared<pybind11::module>("{1}", "{1} module");

	std::vector< std::pair<std::string, std::string> > sub_modules {{
{2}	}};
	for(auto &p : sub_modules ) modules[p.first.size() ? p.first+"::"+p.second : p.second] = std::make_shared<pybind11::module>( modules[p.first]->def_submodule(p.second.c_str(), ("Bindings for " + p.first + "::" + p.second + " namespace").c_str() ) );

{3}
	return modules[""]->ptr();
}}
)_";

const char * module_header = "\n#include <pybind11/pybind11.h>\n\n{}#ifndef BINDER_PYBIND11_TYPE_CASTER\n\t#define BINDER_PYBIND11_TYPE_CASTER\n\tPYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);\n#endif\n\n";

const char * module_function_suffix = "(std::function< pybind11::module &(std::string const &namespace_) > &M)";

void Context::generate(Config const &config)
{
	bind(config);

	vector<string> sources;
	vector<string> binding_function_names;

	std::map<string, int> file_names;

	string file_name = config.prefix + config.root_module + ".cpp";
	std::ofstream root_module_file_handle(file_name);
	sources.push_back(file_name);

	outs() << "Writing code...\n";
	for(uint i=0; i<binders.size(); ++i) {
		if( /*binders[i]->is_binded()  and*/  binders[i]->code().size() ) {
			string np = file_name_prefix_for_binder(binders[i]);
			string file_name = np + ( file_names[np] ? "_"+std::to_string(file_names[np]) : "" );
			++file_names[np];

			string function_name = "bind_" + replace_(file_name, "/", "_");
			file_name += ".cpp"; //".b.cpp";

			sources.push_back(file_name);
			binding_function_names.push_back(function_name);
			//outs() << string(*binders[i]) << " → " << file_name << "\n";

			string code, prefix_code;
			string namespace_ = namespace_from_named_decl( binders[i]->named_decl() );

			//vector<string> includes;
			//std::set<NamedDecl const *> stack;
			IncludeSet includes;

			for(; code.size()<config.maximum_file_length  and  i<binders.size()  and  namespace_==namespace_from_named_decl( binders[i]->named_decl() ); ++i) {
				//outs() << "Binding: " << string(*binders[i]) << "\n";
				if( ClassBinder * CB = dynamic_cast<ClassBinder*>( binders[i].get() ) ) {
					std::vector<clang::CXXRecordDecl const *> const dependencies = CB->dependencies();
					for(auto & c : dependencies ) {
						if( is_forward_needed(c) ) {
							code += bind_forward_declaration(c, *this);
							add_to_binded(c);
							outs() << "Adding forward binding for " << class_qualified_name(c) << "\n";
						}
					}
					add_to_binded( dynamic_cast<CXXRecordDecl*>( CB->named_decl() ) );

					prefix_code += CB->prefix_code();
				}
				code += binders[i]->code();
				binders[i]->add_relevant_includes(includes);
			}
			if( i < binders.size() ) --i;

			code = generate_include_directives(includes) + format(module_header, config.includes_code()) + prefix_code + "void " + function_name + module_function_suffix + "\n{\n" + code + "}\n";

			if( O_single_file ) root_module_file_handle << "// File: " << file_name << '\n' << code << "\n\n";
			else update_source_file(config.prefix, file_name, code);
		}
	}
	outs() << "Writing code... Done!\n";

	string namespace_pairs;
	std::set<string> namespaces = create_all_nested_namespaces();
	for(auto & n : namespaces) {
		if( n.size() ) namespace_pairs += "\t\t{{\"{}\", \"{}\"}},\n"_format(base_namespace(n), last_namespace(n));
	}

	string binding_function_decls, binding_function_calls;
	for(auto &f : binding_function_names) {
		binding_function_decls += "void " + f + module_function_suffix + ";\n";
		binding_function_calls += "\t" + f + "(M);\n";
	}

	std::stringstream s;
	s << format(main_module_header, binding_function_decls, config.root_module, namespace_pairs, binding_function_calls);

	root_module_file_handle << s.str();

	if( O_single_file ) {
		root_module_file_handle << "\n// Source list file: " << config.prefix + config.root_module + ".sources\n";
		for(auto &s : sources) root_module_file_handle << "// "<< s << "\n";
	} else {
		std::ofstream f(config.prefix + config.root_module + ".sources");
		for(auto &s : sources) f << s << "\n";
	}
}

} // namespace binder
