// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   binder/function.cpp
/// @brief  Binding generation for static and member functions
/// @author Sergey Lyskov

#include <function.hpp>

#include <class.hpp>
#include <type.hpp>
#include <util.hpp>
#include <fmt/format.h>

#include <clang/AST/DeclCXX.h>
#include <clang/AST/ASTContext.h>

#include <clang/AST/ExprCXX.h>


#include <vector>


using namespace llvm;
using namespace clang;

using std::string;
using std::pair;
using std::tuple;
using std::vector;
using std::unordered_map;

using namespace fmt::literals;

namespace binder {

static std::map<string, string > const cpp_python_operator_map {
	{"operator+", "__add__"},
	{"operator-", "__sub__"},

	// {"operator+=", "__iadd__"},
	// {"operator-=", "__isub__"},
	{"operator*=", "__imul__"},
	{"operator/=", "__idiv__"},

	{"operator()", "__call__"},
	{"operator==", "__eq__"},
	{"operator!=", "__ne__"},
	{"operator[]", "__getitem__"},
	{"operator=",  "assign"},
	{"operator++", "plus_plus"},
	{"operator--", "minus_minus"},
};

// Generate function argument list separate by comma: int, bool, std::string
string function_arguments(clang::FunctionDecl const *record)
{
	string r;

	for(uint i=0; i<record->getNumParams(); ++i) {
		r += record->getParamDecl(i)->getOriginalType().getCanonicalType().getAsString();
		if( i+1 != record->getNumParams() ) r += ", ";
	}

	fix_boolean_types(r);

	return r;
}


// Generate function argument list separate by comma
// name_arguments - if arguments should be named: a1, a2, ...
// n - number of arguments to generate. If n > num_of_function_parameters - generate only list with num_of_function_parameters
pair<string, string> function_arguments_for_lambda(clang::FunctionDecl const *record, uint n)
{
	string r, a;

	for(uint i=0; i<record->getNumParams()  and  i<n; ++i) {
		QualType qt = record->getParamDecl(i)->getOriginalType().getCanonicalType();
		r += qt.getAsString() + ' ';
		if( !qt->isReferenceType()  and  !qt->isPointerType() ) r += !qt.isConstQualified() ? " const &" : " &";
		r += "a" + std::to_string(i);
		a += "a" + std::to_string(i);
		if( i+1 != record->getNumParams()  and  i+1 != n ) { r += ", ";  a += ", "; }
	}

	fix_boolean_types(r);

	return std::make_pair(r, a);
}


// Generate three version of function argument list: (with types separate by comma, only arguments names, only argument names with by-reference arguments converted to pointers by adding '&'
// example: ("string const & a0, int *a1, float a2", "a0, a1, a2", "&a0, a1, a2")
tuple<string, string, string> function_arguments_for_py_overload(clang::FunctionDecl const *record)
{
	string r, a, p;

	for(uint i=0; i<record->getNumParams(); ++i) {
		QualType qt = record->getParamDecl(i)->getOriginalType().getCanonicalType();
		r += qt.getAsString() + ' ' + "a" + std::to_string(i);
		a += "a" + std::to_string(i);
		p += string(qt->isLValueReferenceType() ? "&" : "" ) + "a" + std::to_string(i);
		if( i+1 != record->getNumParams() ) { r += ", ";  a += ", ";  p += ", ";  }
	}

	fix_boolean_types(r);

	return std::make_tuple(r, a, p);
}


// generate class template specialization for ClassTemplateSpecializationDecl or empty string otherwise
string template_specialization(FunctionDecl const *F)
{
	string templ;

	if( F->getTemplatedKind() == FunctionDecl::TK_MemberSpecialization  or   F->getTemplatedKind() == FunctionDecl::TK_FunctionTemplateSpecialization ) {
		if( TemplateArgumentList const *ta = F->getTemplateSpecializationArgs() ) {
			templ += "<";
			for(uint i=0; i < ta->size(); ++i) {
				//outs() << "function template argument: " << template_argument_to_string( ta->get(i) ) << "\n";
				templ += template_argument_to_string( ta->get(i) ) + ",";

				//if( t->getTemplateArgs()[i].ArgKind() == TemplateArgument::ArgKind::Integral ) outs() << " template arg:" << t->getTemplateArgs()[i].<< "\n";
				//outs() << expresion_to_string( t->getTemplateArgs()[i].getAsExpr() ) << "\n";
			}
			templ.back() = '>';

			fix_boolean_types(templ);
		}
	}

	return templ;
}


// generate string represetiong class name that could be used in python
string python_function_name(FunctionDecl const *F)
{
	if( F->isOverloadedOperator() ) return cpp_python_operator_map.at( F->getNameAsString() );
	else return mangle_type_name( F->getNameAsString() + template_specialization(F) );
}

// Generate function pointer type string for given function: void (*)(int, doule)_ or  void (ClassName::*)(int, doule)_ for memeber function
string function_pointer_type(FunctionDecl const *F)
{
	string r;
	string prefix, maybe_const;
	if( auto m = dyn_cast<CXXMethodDecl>(F) ) {
		prefix = m->isStatic() ? "" : class_qualified_name( cast<CXXRecordDecl>( F->getParent() ) ) + "::";
	    maybe_const = m->isConst() ? " const" : "";
	}

	r += F->getReturnType().getCanonicalType().getAsString();  r+= " ({}*)("_format(prefix);

	r += function_arguments(F);

	r += ")" + maybe_const;

	fix_boolean_types(r);

	return r;
}


// generate qualified function name that could be used in bindings code indcluding template specialization if any
string function_qualified_name(FunctionDecl const *F)
{
	string maybe_const;
	if( auto m = dyn_cast<CXXMethodDecl>(F) ) maybe_const = m->isConst() ? " const" : "";

	string r = F->getReturnType().getCanonicalType().getAsString() + " "+ standard_name( F->getQualifiedNameAsString() + template_specialization(F) ) + "(" + function_arguments(F) + ")" + maybe_const;
	fix_boolean_types(r);
	return r;
}


// generate vector<QualType> with all types that function uses including: return type, types of function arguments and template arguments
vector<QualType> get_type_dependencies(FunctionDecl const *F)
{
	vector<QualType> r;

	r.push_back( F->getReturnType() ); //.getDesugaredType(F->getASTContext()) );
	for(uint i=0; i<F->getNumParams(); ++i) r.push_back(F->getParamDecl(i)->getOriginalType()/*.getDesugaredType(F->getASTContext())*/ );

	// r.push_back( F->getReturnType().getDesugaredType(F->getASTContext()) );
	// for(uint i=0; i<F->getNumParams(); ++i) r.push_back(F->getParamDecl(i)->getOriginalType().getDesugaredType(F->getASTContext()) );

	if( F->getTemplatedKind() == FunctionDecl::TK_MemberSpecialization  or   F->getTemplatedKind() == FunctionDecl::TK_FunctionTemplateSpecialization ) {
		if( TemplateArgumentList const *tal = F->getTemplateSpecializationArgs() ) {
			for(uint i=0; i < tal->size(); ++i) {
				TemplateArgument const &ta( tal->get(i) );
				if( ta.getKind() == TemplateArgument::Type ) r.push_back( ta.getAsType() );
			}
		}
	}

	return r;
}


/// check if user requested binding for the given declaration
bool is_binding_requested(FunctionDecl const *F, Config const &config)
{
	bool bind = config.is_function_binding_requested( F->getQualifiedNameAsString() ) or  config.is_function_binding_requested( function_qualified_name(F) )  or  config.is_namespace_binding_requested( namespace_from_named_decl(F) );

	for(auto & t : get_type_dependencies(F) ) bind |= binder::is_binding_requested(t, config);

	return bind;
}


/// check if user requested skipping for the given declaration
bool is_skipping_requested(FunctionDecl const *F, Config const &config)
{
	string name = standard_name( F->getQualifiedNameAsString() );
	bool skip = config.is_function_skipping_requested(name) or config.is_function_skipping_requested( function_qualified_name(F) ) or config.is_namespace_skipping_requested( namespace_from_named_decl(F) );

    // moved to config -> name.erase(std::remove(name.begin(), name.end(), ' '), name.end());
	skip |= config.is_function_skipping_requested(name);

	// calculating skipping for template classes without template specialization specified as: myclass::member_function_to_skip
	//outs() << "Checking skipping for function: " << name << "... ";
	if( CXXMethodDecl const *M = dyn_cast<CXXMethodDecl>(F) ) {
		CXXRecordDecl const *C = M->getParent();
		if( dyn_cast<ClassTemplateSpecializationDecl>(C) ) {
			//outs() << C->getQualifiedNameAsString() << "::" << F->getNameAsString() << "\n";
			skip |= config.is_function_skipping_requested( standard_name( C->getQualifiedNameAsString() + "::" + F->getNameAsString() ) );
		}
	}
	//outs() << "OK\n";

	for(auto & t : get_type_dependencies(F) ) skip |= is_skipping_requested(t, config);

	return skip;
}


// Generate binding for given function: .def("foo", (std::string (aaaa::A::*)(int) ) &aaaa::A::foo, "doc")
string bind_function(FunctionDecl const *F, uint args_to_bind, bool request_bindings_f, Context &context)
{
	string function_name = python_function_name(F);
	string function_qualified_name { F->getQualifiedNameAsString() };

	CXXMethodDecl const * m = dyn_cast<CXXMethodDecl>(F);
	string maybe_static = m and m->isStatic() ? "_static" : "";


	string function;
	if( args_to_bind == F->getNumParams() ) {
		function = "({}) &{}{}"_format(function_pointer_type(F), function_qualified_name, template_specialization(F));
	}
	else {
		pair<string, string> args = function_arguments_for_lambda(F, args_to_bind);
		//string args; for(uint i=0; i<args_to_bind; ++i) args += "a" + std::to_string(i) + ( i+1 == args_to_bind ? "" : ", " );

		string return_type = F->getReturnType().getCanonicalType().getAsString();  fix_boolean_types(return_type);

		if( m and !m->isStatic() ) {
			string object = class_qualified_name( m->getParent() ) + (m->isConst() ? " const" : "") + " &o" + ( args_to_bind ? ", " : "" );
			function = "[]({}{}) -> {} {{ return o.{}({}); }}"_format(object, args.first, return_type, F->getNameAsString(), args.second);
		}
		else {
			function = "[]({}) -> {} {{ return {}({}); }}"_format(args.first, return_type, function_qualified_name, args.second);
		}
	}

	string maybe_return_policy = "";
	if     ( F->getReturnType()->isPointerType() )         maybe_return_policy = ", " + Config::get().default_pointer_return_value_policy();
	else if( F->getReturnType()->isLValueReferenceType() ) maybe_return_policy = ", " + Config::get().default_lvalue_reference_return_value_policy();
	else if( F->getReturnType()->isRValueReferenceType() ) maybe_return_policy = ", " + Config::get().default_rvalue_reference_return_value_policy();

	//string r = R"(.def{}("{}", ({}) &{}{}, "doc")"_format(maybe_static, function_name, function_pointer_type(F), function_qualified_name, template_specialization(F));
	string r = R"(.def{}("{}", {}, "doc"{})"_format(maybe_static, function_name, function, maybe_return_policy);

	if(request_bindings_f) request_bindings(F->getReturnType().getCanonicalType(), context);

	for(uint i=0; i<F->getNumParams()  and  i < args_to_bind; ++i) {
		r += ", pybind11::arg(\"{}\")"_format( string( F->getParamDecl(i)->getName() ) );

		if(request_bindings_f) request_bindings( F->getParamDecl(i)->getOriginalType(), context);
	}

	// for(auto p = F->param_begin(); p != F->param_end(); ++p) {
	// 	string default_argument;
	// 	if( !(*p)->hasUninstantiatedDefaultArg() ) default_argument = expresion_to_string( (*p)->getDefaultArg() );
	// 	bool is_function_call = ( default_argument.find("(") != std::string::npos  and  default_argument.find(")") != std::string::npos )  or  default_argument.find("new ") != std::string::npos;  // filter 'function call' default arguments
	// 	string arg_type = (*p)->getOriginalType().getCanonicalType().getAsString();  fix_boolean_types(arg_type);
	// 	bool good_default = default_argument.find("std::cout") == std::string::npos  and  default_argument.find("std::cerr") == std::string::npos;
	// 	if( (*p)->hasDefaultArg()  and  !(*p)->hasUninstantiatedDefaultArg()  and  !is_function_call  and  good_default  and  false) {
	// 		r += ", pybind11::arg_t<{}>(\"{}\", {}, \"{}\")"_format(arg_type, string( (*p)->getName() ), default_argument, replace(default_argument, "\"", "\\\"") );
	// 	}
	// 	else r += ", pybind11::arg(\"{}\")"_format( string( (*p)->getName() ) );
	// 	if(request_bindings_f) request_bindings((*p)->getOriginalType(), context);
	// }

	r += ");";

	return r;
}

// Generate binding for given function. If function have default arguments generate set of bindings by creating separate bindings for each argument with default.
string bind_function(string const & module, FunctionDecl const *F, Context &context)
{
	string code;

	uint args_to_bind = 0;
	for(; args_to_bind < F->getNumParams(); ++args_to_bind) {
		if( F->getParamDecl(args_to_bind)->hasDefaultArg() ) break;
	}

	for(; args_to_bind <= F->getNumParams(); ++args_to_bind) code += module + bind_function(F, args_to_bind, args_to_bind == F->getNumParams(), context) + '\n';

	return code;
	// string function_name = python_function_name(F);
	// string function_qualified_name { F->getQualifiedNameAsString() };

	// CXXMethodDecl * m = dyn_cast<CXXMethodDecl>(F);
	// string maybe_static = m and m->isStatic() ? "_static" : "";

	// string r = R"(.def{}("{}", ({}) &{}{}, "doc")"_format(maybe_static, function_name, function_pointer_type(F), function_qualified_name, template_specialization(F));

	// request_bindings(F->getReturnType(), context);

	// for(auto p = F->param_begin(); p != F->param_end(); ++p) {

	// 	string default_argument;
	// 	if( !(*p)->hasUninstantiatedDefaultArg() ) default_argument = expresion_to_string( (*p)->getDefaultArg() );

	// 	bool is_function_call = ( default_argument.find("(") != std::string::npos  and  default_argument.find(")") != std::string::npos )  or  default_argument.find("new ") != std::string::npos;  // filter 'function call' default arguments

	// 	string arg_type = (*p)->getOriginalType().getCanonicalType().getAsString();  fix_boolean_types(arg_type);

	// 	bool good_default = default_argument.find("std::cout") == std::string::npos  and  default_argument.find("std::cerr") == std::string::npos;

	// 	if( (*p)->hasDefaultArg()  and  !(*p)->hasUninstantiatedDefaultArg()  and  !is_function_call  and  good_default  and  false) {
	// 		r += ", pybind11::arg_t<{}>(\"{}\", {}, \"{}\")"_format(arg_type, string( (*p)->getName() ), default_argument, replace(default_argument, "\"", "\\\"") );
	// 	}
	// 	else r += ", pybind11::arg(\"{}\")"_format( string( (*p)->getName() ) );

	// 	request_bindings((*p)->getOriginalType(), context);
	// }

	// r += ')';

	// return r;
}


/// extract include needed for this generator and add it to includes vector
void add_relevant_includes(FunctionDecl const *F, IncludeSet &includes, int level /*, bool for_template_arg_only*/)
{
	//if( stack.count(F) ) return; else stack.insert(F);
	if( !includes.add_decl(F, level) ) return;

	// if( begins_with(F->getQualifiedNameAsString(), "boost::get_property_value") ) outs() << "add_relevant_includes(function): " << F->getQualifiedNameAsString() << " templ:" << template_specialization(F) << "\n";

	add_relevant_include_for_decl(F, includes);

	for(auto & t : get_type_dependencies(F) ) binder::add_relevant_includes(t, includes, level);
}


/// Generate string id that uniquly identify C++ binding object. For functions this is function prototype and for classes forward declaration.
string FunctionBinder::id() const
{
	return function_qualified_name(F);
}


/// check if generator can create binding
bool is_bindable(FunctionDecl const *F)
{
	//bool r = true;
	bool r = !F->isDeleted();

	if( F->isOverloadedOperator() ) {
		//outs() << "Operator: " << F->getNameAsString() << '\n';
		if( !isa<CXXMethodDecl>(F)  or  !cpp_python_operator_map.count( F->getNameAsString() ) ) return false;
	}

	r &= F->getTemplatedKind() != FunctionDecl::TK_FunctionTemplate  /*and  !F->isOverloadedOperator()*/  and  !isa<CXXConversionDecl>(F)  and  !F->isDeleted();

	QualType rt( F->getReturnType() );

	r &= is_bindable(rt);

	for(auto p = F->param_begin(); p != F->param_end(); ++p) r &= is_bindable( (*p)->getOriginalType().getCanonicalType() );
	//outs() << "is_bindable: " << F->getQualifiedNameAsString() << " " << r << "\n";

	return r;
}



bool FunctionBinder::bindable() const
{
	return binder::is_bindable(F);
}


/// check if user requested binding for the given declaration
void FunctionBinder::request_bindings_and_skipping(Config const &config)
{
	if( is_skipping_requested(F, config) ) Binder::request_skipping();
	else if( is_binding_requested(F, config) ) Binder::request_bindings();
}


/// extract include needed for this generator and add it to includes vector
void FunctionBinder::add_relevant_includes(IncludeSet &includes) const
{
	binder::add_relevant_includes(F, includes, 0);
}


/// generate binding code for this object and all its dependencies
void FunctionBinder::bind(Context &context)
{
	if( is_binded() ) return;

	string const module_variable_name = context.module_variable_name( namespace_from_named_decl(F) );
	string const include = relevant_include(F);
	//string const namespace_ = namespace_from_named_decl(F);

	code()  = "\t// " + F->getQualifiedNameAsString() + "(" + function_arguments(F) + ") file:" + (include.size() ? include.substr(1, include.size()-2) : "") + " line:" + line_number(F) + "\n";
	//if( namespace_.size() ) code() += "\t\tusing namespace " + namespace_ + ";\n";
	code() += bind_function("\t"+ module_variable_name, F, context);
	code() += "\n";
}

} // namespace binder
