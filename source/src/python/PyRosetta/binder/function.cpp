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
using std::vector;
using std::unordered_map;

using namespace fmt::literals;

namespace binder {


// Generate function argument list separate by comma: int, bool, std::sting
string function_arguments(clang::FunctionDecl const *record)
{
	string r;

	for(uint i=0; i<record->getNumParams(); ++i) {
		r += record->getParamDecl(i)->getOriginalType().getCanonicalType().getAsString();
		if( i != record->getNumParams()-1 ) r += ", ";
	}

	fix_boolean_types(r);

	return r;
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

	string r = F->getReturnType().getCanonicalType().getAsString() + " "+ F->getQualifiedNameAsString() + template_specialization(F) + "(" + function_arguments(F) + ")" + maybe_const;
	fix_boolean_types(r);
	return r;
}


// generate vector<QualType> with all types that function uses including: return type, types of function arguments and template arguments
vector<QualType> get_type_dependencies(FunctionDecl const *F)
{
	vector<QualType> r;

	r.push_back( F->getReturnType() ); //.getDesugaredType(F->getASTContext()) );

	for(uint i=0; i<F->getNumParams(); ++i) r.push_back(F->getParamDecl(i)->getOriginalType() );

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


// /// check if user requested binding for the given declaration
// bool is_binding_requested(FunctionDecl const *F, Config const &config)
// {
// 	bool bind = config.is_function_binding_requested( F->getQualifiedNameAsString() )  or  config.is_namespace_binding_requested( namespace_from_named_decl(F) );
// 	if(bind) {
// 		bind &= !binder::is_skipping_requested(F->getReturnType().getDesugaredType(F->getASTContext()), config);
// 		for(uint i=0; i<F->getNumParams(); ++i) {
// 			QualType qt = F->getParamDecl(i)->getOriginalType();
// 			bind &= !binder::is_skipping_requested( qt.getDesugaredType(F->getASTContext()), config);
// 		}
// 		if( F->getTemplatedKind() == FunctionDecl::TK_MemberSpecialization  or   F->getTemplatedKind() == FunctionDecl::TK_FunctionTemplateSpecialization ) {
// 			if( TemplateArgumentList const *tal = F->getTemplateSpecializationArgs() ) {
// 				for(uint i=0; i < tal->size(); ++i) {
// 					TemplateArgument const &ta( tal->get(i) );
// 					if( ta.getKind() == TemplateArgument::Type ) {
// 						bind &= !binder::is_skipping_requested( ta.getAsType(), config);
// 					}
// 				}
// 			}
// 		}
// 	}
// 	return bind;
// }


/// check if user requested skipping for the given declaration
bool is_skipping_requested(FunctionDecl const *F, Config const &config)
{
	bool skip = config.is_function_skipping_requested( F->getQualifiedNameAsString() ) or config.is_class_skipping_requested( function_qualified_name(F) ) or config.is_namespace_skipping_requested( namespace_from_named_decl(F) );

	for(auto & t : get_type_dependencies(F) ) skip |= is_skipping_requested(t, config);

	return skip;
}

// /// check if user requested skipping for the given declaration
// bool is_skipping_requested(FunctionDecl const *F, Config const &config)
// {
// 	bool skip = config.is_function_skipping_requested( F->getQualifiedNameAsString() ) or config.is_class_skipping_requested( function_qualified_name(F) ) or config.is_namespace_skipping_requested( namespace_from_named_decl(F) );
// 	if(!skip) {
// 		skip |= binder::is_skipping_requested(F->getReturnType().getDesugaredType(F->getASTContext()), config);
// 		for(uint i=0; i<F->getNumParams(); ++i) {
// 			QualType qt = F->getParamDecl(i)->getOriginalType();
// 			skip |= binder::is_skipping_requested( qt.getDesugaredType(F->getASTContext()), config);
// 		}
// 		if( F->getTemplatedKind() == FunctionDecl::TK_MemberSpecialization  or   F->getTemplatedKind() == FunctionDecl::TK_FunctionTemplateSpecialization ) {
// 			if( TemplateArgumentList const *tal = F->getTemplateSpecializationArgs() ) {
// 				for(uint i=0; i < tal->size(); ++i) {
// 					TemplateArgument const &ta( tal->get(i) );
// 					if( ta.getKind() == TemplateArgument::Type ) {
// 						skip |= binder::is_skipping_requested( ta.getAsType(), config);
// 					}
// 				}
// 			}
// 		}
// 	}
// 	return skip;
// }


// Generate binding for given function: .def("foo", (std::string (aaaa::A::*)(int) ) &aaaa::A::foo, "doc")
string bind_function(FunctionDecl *F, Context &context)
{
	string function_name { F->getNameAsString() };
	string function_qualified_name { F->getQualifiedNameAsString() };

	CXXMethodDecl * m = dyn_cast<CXXMethodDecl>(F);
	string maybe_static = m and m->isStatic() ? "_static" : "";

	string r = R"(.def{}("{}", ({}) &{}{}, "doc")"_format(maybe_static, function_name, function_pointer_type(F), function_qualified_name, template_specialization(F));

	request_bindings(F->getReturnType(), context);

	for(auto p = F->param_begin(); p != F->param_end(); ++p) {
		// .def("myFunction", py::arg("arg") = SomeType(123));
		//string defalt_argument = (*p)->hasDefaultArg() ? " = ({})({})"_format( (*p)->getOriginalType().getCanonicalType().getAsString(), expresion_to_string( (*p)->getDefaultArg() ) ) : "";
		//string defalt_argument = (*p)->hasDefaultArg() ? expresion_to_string( (*p)->getDefaultArg() ) : "";
		//r += ", pybind11::arg(\"{}\"){}"_format( string( (*p)->getName() ), defalt_argument );

		// .def("myFunction", py::arg_t<int>("arg", 123, "123"));
		if( (*p)->hasDefaultArg()  and  !(*p)->hasUninstantiatedDefaultArg() ) {
			string arg_type = (*p)->getOriginalType().getCanonicalType().getAsString();  fix_boolean_types(arg_type);
			r += ", pybind11::arg_t<{}>(\"{}\", {}, \"{}\")"_format(arg_type,
																	string( (*p)->getName() ),
																	expresion_to_string( (*p)->getDefaultArg() ),
																	replace(expresion_to_string( (*p)->getDefaultArg() ), "\"", "\\\"") );

			// clang::Expr *e = (*p)->getDefaultArg();
			// clang::Expr::EvalResult result;
			// outs() << "ex: " << expresion_to_string(e) << " EvaluateAsRValue: " << e->EvaluateAsRValue(result, F->getASTContext());
			// outs() << " res: " << result.Val.getAsString(F->getASTContext(), (*p)->getOriginalType()) << "\n";
		}
		else r += ", pybind11::arg(\"{}\")"_format( string( (*p)->getName() ) );

		//add_relevant_include(*p, includes);

		//outs() << (*p)->getDefaultArg()->getAsString();
		//(*p)->dump();
		// if( (*p)->hasDefaultArg() ) {
		// 	outs() << "  CXXDefaultArgExpr: " << (*p)->getName() << " = " << expresion_to_string( (*p)->getDefaultArg() ) << "\n";

		// 	//SourceRange s = (*p)->getDefaultArgRange();
		// 	//s.getBegin().dump( F->getASTContext().getSourceManager() );


		// 	//if( auto e = dyn_cast<clang::CXXDefaultArgExpr>( (*p)->getDefaultArg() ) ) {
		// 		//e->dump();
		// 	// }

		// 	// Expr::EvalResult result;
		// 	// if( (*p)->getDefaultArg()->EvaluateAsRValue(result, F->getASTContext() ) ) {
		// 	// 	outs() << "  Default for: " << (*p)->getName() << " = " << result.Val.getAsString(F->getASTContext(), (*p)->getOriginalType() ) << "\n";
		// 	// }
		// }
		request_bindings((*p)->getOriginalType(), context);
	}

	r += ')';

	return r;
}


/// extract include needed for this generator and add it to includes vector
void add_relevant_includes(FunctionDecl const *F, std::vector<std::string> &includes, std::set<NamedDecl const *> &stack /*, bool for_template_arg_only*/)
{
	if( stack.count(F) ) return; else stack.insert(F);

	add_relevant_include_for_decl(F, includes);

	for(auto & t : get_type_dependencies(F) ) binder::add_relevant_includes(t, includes, stack);
}


/// extract include needed for this generator and add it to includes vector
// void add_relevant_includes(FunctionDecl const *F, std::vector<std::string> &includes, std::set<NamedDecl const *> &stack /*, bool for_template_arg_only*/)
// {
// 	if( stack.count(F) ) return;
// 	else stack.insert(F);
// 	add_relevant_include_for_decl(F, includes);
// 	//if( !for_template_arg_only ) {
// 	binder::add_relevant_includes(F->getReturnType().getDesugaredType(F->getASTContext()), includes, stack);
// 	for(uint i=0; i<F->getNumParams(); ++i) {
// 		QualType qt = F->getParamDecl(i)->getOriginalType();
// 		//binder::add_relevant_includes(qt, includes);
// 		binder::add_relevant_includes( qt.getDesugaredType(F->getASTContext()), includes, stack);
// 	}
// 	if( F->getTemplatedKind() == FunctionDecl::TK_MemberSpecialization  or   F->getTemplatedKind() == FunctionDecl::TK_FunctionTemplateSpecialization ) {
// 		if( TemplateArgumentList const *tal = F->getTemplateSpecializationArgs() ) {
// 			for(uint i=0; i < tal->size(); ++i) {
// 				TemplateArgument const &ta( tal->get(i) );
// 				if( ta.getKind() == TemplateArgument::Type ) {
// 					binder::add_relevant_includes( ta.getAsType(), includes, stack);
// 				}
// 			}
// 		}
// 	}
// }


/// Generate string id that uniquly identify C++ binding object. For functions this is function prototype and for classes forward declaration.
string FunctionBinder::id() const
{
	return function_qualified_name(F);
}


/// check if generator can create binding
bool is_bindable(FunctionDecl const *F)
{
	bool r = true;

	// todo: bindging for operators and type conversion
	r &= F->getTemplatedKind() != FunctionDecl::TK_FunctionTemplate  and  !F->isOverloadedOperator()  and  !isa<CXXConversionDecl>(F)  and  !F->isDeleted();

	QualType rt( F->getReturnType() );

	r &= is_bindable(rt);

	for(auto p = F->param_begin(); p != F->param_end(); ++p) r &= is_bindable( (*p)->getOriginalType() );

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
void FunctionBinder::add_relevant_includes(std::vector<std::string> &includes, std::set<clang::NamedDecl const *> &stack) const
{
	binder::add_relevant_includes(F, includes, stack);
}


/// generate binding code for this object and all its dependencies
void FunctionBinder::bind(Context &context)
{
	if( is_binded() ) return;

	string const module_variable_name = context.module_variable_name( namespace_from_named_decl(F) );
	string const include = relevant_include(F);

	code()  = "\t{ // " + F->getQualifiedNameAsString() + "(" + function_arguments(F) + ") file:" + include.substr(1, include.size()-2) + " line:" + line_number(F) + "\n";
	code() += "\t\tusing namespace " + namespace_from_named_decl(F) + ";\n";
	code() += "\t\t"+ module_variable_name + bind_function(F, context) + ";\n";
	code() += "\t}\n\n";
}

} // namespace binder
