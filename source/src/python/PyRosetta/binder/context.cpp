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

#include <type.hpp>
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

using namespace fmt::literals;
using fmt::format;

namespace binder {

//const std::string _module_variable_name_{"M"};

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


	//string ns = namespace_from_named_decl( b->named_decl() );
	//modules[ns].push_back(b);

	//outs() << "Adding to: " << ns << "\n";
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
	//namespaces.assign( s.begin(), s.end() );

	//std::sort( std::begin(namespaces), std::end(namespaces) );

	return s;
}

// /// create vector of all namespaces and sort it
// vector<string> Context::sorted_namespaces()
// {
// 	vector<string> n;
// 	for(auto & p : modules) n.push_back(p.first);

// 	std::sort( std::begin(n), std::end(n) );

// 	return n;
// }


// generate code for include directives and cleanup the includes vector
string generate_include_directives(vector<string> &includes)
{
	string r;
	for(auto &i : std::set<string>(includes.begin(), includes.end() ) ) r += "#include " + i + '\n';
	includes.resize(0);
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
		//outs() << "Context::bind: requesting bindings for: " << type << "...\n";
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

	if( include.size() <= 2 ) throw std::runtime_error( "file_name_for_decl failed!!! include name for decl: " + string(*b) + " is too short!");
	include = include.substr(1, include.size()-2);

	//outs() << "qn:" << decl->getQualifiedNameAsString() << " namespace_from_named_decl(decl):" << namespace_from_named_decl(decl) << "\n";

	if( namespace_from_named_decl(decl) == "std" ) include = "std/" +  ( begins_with(include, "bits/") ? include.substr(5) : include );

	return replace( replace( replace( replace(include, ".hh", ""), ".hpp", ""), ".h", ""), ".", "_");
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

const char * module_header = "\n#include <pybind11/pybind11.h>\n#include <pybind11/stl.h>\n\n//PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);\n";

const char * module_function_suffix = "(std::function< pybind11::module &(std::string const &namespace_) > &M)";

void Context::generate(Config const &config)
{
	bind(config);

	vector<string> sources;
	vector<string> binding_function_names;

	std::map<string, int> file_names;

	outs() << "Writing code...\n";
	for(uint i=0; i<binders.size(); ++i) {
		if( /*binders[i]->is_binded()  and*/  binders[i]->code().size() ) {
			string np = file_name_prefix_for_binder(binders[i]);
			string file_name = np + ( file_names[np] ? "_"+std::to_string(file_names[np]) : "" );
			++file_names[np];

			string function_name = "bind_" + replace(file_name, "/", "_");
			file_name += ".cpp"; //".b.cpp";

			sources.push_back(file_name);
			binding_function_names.push_back(function_name);
			//outs() << string(*binders[i]) << " → " << file_name << "\n";

			string code;
			string namespace_ = namespace_from_named_decl( binders[i]->named_decl() );
			vector<string> includes;
			std::set<NamedDecl const *> stack;

			for(; code.size()<config.maximum_file_length  and  i<binders.size()  and  namespace_==namespace_from_named_decl( binders[i]->named_decl() ); ++i) {
				//outs() << "Binding: " << string(*binders[i]) << "\n";
				code += binders[i]->code();
				binders[i]->add_relevant_includes(includes, stack);
			}
			if( i < binders.size() ) --i;

			code = generate_include_directives(includes) + module_header + "void " + function_name + module_function_suffix + "\n{\n" + code + "}\n";

			update_source_file(config.prefix, file_name, code);
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


	//s << "#include <pybind11/pybind11.h>\nPYBIND11_PLUGIN(example) {\n\tpybind11::module m(\"example\", \"example module\");\n";

	// string module_var_defs, module_function_defs, module_function_calls;

	// for(auto & n : namespaces) {
	// 	if( n.size() ) {
	// 		module_var_defs += "\tpybind11::module {} = {}.def_submodule(\"{}\", \"Bindings for {} namespace\");\n"_format( module_variable(n), module_variable( base_namespace(n) ), last_namespace(n), n);
	// 	}

	// 	module_function_defs += module_binder_function(n, true) + ";\n";
	// 	module_function_calls += module_binder_function(n, false, "", module_variable(n)) + ";\n";
	// }

	//s << format(main_module_header, module_function_defs, root_module, _root_module_variable_name_, root_module, root_module + " module");

	// s << module_var_defs << "\n" << module_function_calls;
	// s << "\treturn m.ptr();\n}\n";

	// //for(auto & p : modules) outs() << "Ns: " << p.first << "\n";

	// // vector<BinderOP> & binders( modules[""] );
	// // //for(auto &b : binders) s << indent(b("m"), "\t");
	// // for(auto &b : binders) s << (*b)("m", "\t");
	// // s << "\treturn m.ptr();\n}\n";

	// for(auto &ns : modules) {
	// 	vector<string> files = bind_namespaces(ns.first, maximum_file_length);

	// 	//for(auto &f : files) s << f;
	// 	for(uint i=0; i<files.size(); ++i) {
	// 		string file_name = prefix + ( ns.first.size() ? replace(ns.first, "::", "-") : root_module ) + "-" + std::to_string(i) + ".cpp";

	// 		std::ofstream(file_name) << files[i];
	// 		sources.push_back(file_name);
	// 	}
	// }

	//errs() << s.str() << "\n";
	string file_name = config.prefix + config.root_module + ".cpp";
	std::ofstream(file_name) << s.str();
	sources.push_back(file_name);

	std::ofstream f(config.prefix + config.root_module + ".sources");
	for(auto &s : sources) f << s << "\n";
}

} // namespace binder
