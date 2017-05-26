// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// Copyright (c) 2016 Sergey Lyskov <sergey.lyskov@jhu.edu>
//
// All rights reserved. Use of this source code is governed by a
// MIT license that can be found in the LICENSE file.

/// @file   binder/class.cpp
/// @brief  Binding generation for C++ struct and class objects
/// @author Sergey Lyskov

#include <class.hpp>
#include <function.hpp>
#include <enum.hpp>
#include <type.hpp>
#include <util.hpp>

#include <fmt/format.h>

#include <clang/AST/DeclTemplate.h>
//#include <clang/AST/TemplateBase.h>


using namespace llvm;
using namespace clang;

using std::string;
using std::pair;
using std::tuple;
using std::make_pair;
using std::vector;
using std::set;

//using std::unordered_map;

using namespace fmt::literals;

namespace binder {

// generate class template specialization for ClassTemplateSpecializationDecl or empty string otherwise
string template_specialization(clang::CXXRecordDecl const *C)
{
	string templ;

	if( auto t = dyn_cast<ClassTemplateSpecializationDecl>(C) ) {
		templ += "<";

		for(uint i=0; i < t->getTemplateArgs().size(); ++i) {
			//if( t->getTemplateArgs()[i].isInstantiationDependent() ) break;  // avoid explicitly specifying SFINAE related arguments

			//outs() << " template argument: " << template_argument_to_string(t->getTemplateArgs()[i]) << "\n";
			templ += template_argument_to_string(t->getTemplateArgs()[i]) + ",";

			//if( t->getTemplateArgs()[i].ArgKind() == TemplateArgument::ArgKind::Integral ) outs() << " template arg:" << t->getTemplateArgs()[i].<< "\n";
			//outs() << expresion_to_string( t->getTemplateArgs()[i].getAsExpr() ) << "\n";
		}
		templ.back() = '>';
	}

	fix_boolean_types(templ);
	return templ;
}

// generate class name that could be used in bindings code indcluding template specialization if any
string class_name(CXXRecordDecl const *C)
{
	string res = standard_name( C->getNameAsString() + template_specialization(C) );

	if(  namespace_from_named_decl(C) == "std" ) res = simplify_std_class_name(res);

	return res;

	//return standard_name(C->getNameAsString() + template_specialization(C));
}


// generate string represetiong class name that could be used in python
string python_class_name(CXXRecordDecl const *C)
{
	string name = class_name(C);
	return mangle_type_name(name);
}


// generate qualified class name that could be used in bindings code indcluding template specialization if any
string class_qualified_name(CXXRecordDecl const *C)
{
	return standard_name(C->getQualifiedNameAsString() + template_specialization(C));
}

// generate vector<QualType> with all types that class uses if it tempalated
vector<QualType> get_type_dependencies(CXXRecordDecl const *C /*, bool include_members*/)
{
	vector<QualType> r;

	if( auto t = dyn_cast<ClassTemplateSpecializationDecl>(C) ) {
		for(uint i=0; i < t->getTemplateArgs().size(); ++i) {
			if( t->getTemplateArgs()[i].getKind() == TemplateArgument::Type ) {
				r.push_back( t->getTemplateArgs()[i].getAsType() /*.getDesugaredType(C->getASTContext())*/ );
			}
		}
	}

	return r;
}


// Return true if class have direct or inderect std::enable_shared_from_this as base class
bool is_inherited_from_enable_shared_from_this(CXXRecordDecl const *C)
{
	if( C->getQualifiedNameAsString() == "std::enable_shared_from_this" ) return true;
	if( C->isCompleteDefinition() ) {
		for(auto b = C->bases_begin(); b!=C->bases_end(); ++b) {
			if( auto r = dyn_cast<RecordType>(b->getType().getCanonicalType().getTypePtr() ) ) {
				CXXRecordDecl *rd = cast<CXXRecordDecl>(r->getDecl());
			if( rd /*and rd->isCompleteDefinition()*/  and  is_inherited_from_enable_shared_from_this(rd) ) return true;
			}
		}
	}
	return false;
}


bool is_field_assignable(FieldDecl const *f)
{
	if( RecordType const* r = dyn_cast<RecordType>( f->getType() ) ) {
		if( CXXRecordDecl *C = cast<CXXRecordDecl>(r->getDecl() ) ) { // checking if this type has deleted operator=
			for(auto m = C->method_begin(); m != C->method_end(); ++m) {
				//if( m->getAccess() == AS_public  and  m->isCopyAssignmentOperator()  /*and  !m->doesThisDeclarationHaveABody()*/ and  m->isDeleted() ) return false;
				if(  m->isCopyAssignmentOperator()  and  ( m->getAccess() != AS_public  or  m->isDeleted() ) ) return false;
			}
		}
	}

	return true;
}

/// check if generator can create binding
bool is_bindable(FieldDecl *f)
{
	if( f->getType()->isAnyPointerType() or f->getType()->isReferenceType()  or  f->getType()->isArrayType() ) return false;

	if( !is_field_assignable(f) ) return false;

	return true;
}


// Generate bindings for class data member
string bind_data_member(FieldDecl const *d, string const &class_qualified_name)
{
	if( d->getType().isConstQualified() ) return ".def_readonly(\"{}\", &{}::{})"_format(d->getNameAsString(), class_qualified_name, d->getNameAsString());
	else return ".def_readwrite(\"{}\", &{}::{})"_format(d->getNameAsString(), class_qualified_name, d->getNameAsString());
}



/// check if generator can create binding
bool is_bindable(clang::CXXRecordDecl const *C)
{
	bool r = true;

	// outs() << "is_bindable(CXXRecordDecl): " << C->getQualifiedNameAsString() << template_specialization(C)
	// 	   << " C->hasDefinition():" << C->hasDefinition()
	// 	   << " C->isCompleteDefinition():" << C->isCompleteDefinition()
	// 	   // << " C->isThisDeclarationADefinition():" << C->isThisDeclarationADefinition()
	// 	   // << " C->getDefinition():" << C->getDefinition()
	// 	   << " C->isDependentType():" << C->isDependentType()
	// 	   <<"\n";
	string qualified_name = C->getQualifiedNameAsString();
	// if( qualified_name != "std::pair"  and  qualified_name != "std::tuple" ) {
	// 	if( C->isDependentType() ) return false;
	// 	if( !C->isCompleteDefinition() ) return false;
	// 	if( C->getAccess() == AS_protected  or  C->getAccess() == AS_private ) return false;
	// }

	if( qualified_name == "(anonymous)" ) return false;
	if( C->isDependentType() ) return false;
	if( C->getAccess() == AS_protected  or  C->getAccess() == AS_private ) return false;

	if( !C->isCompleteDefinition() ) {
		if( auto ts = dyn_cast<ClassTemplateSpecializationDecl>(C) ) {
			if( ts->getPointOfInstantiation()/* SourceLocation */.isInvalid()  and  not is_python_builtin(C) ) {
				//outs() << "is_bindable( " << class_qualified_name(C) << " ): no point of instantiation  found, skipping...\n";
				return false;
			}
		}
		else return false;
	}

	if( auto t = dyn_cast<ClassTemplateSpecializationDecl>(C) ) {
		for(uint i=0; i < t->getTemplateArgs().size(); ++i) {

			if( t->getTemplateArgs()[i].getKind() == TemplateArgument::Type ) {
				if( !is_bindable( t->getTemplateArgs()[i].getAsType() ) ) return false;
			}

			if( t->getTemplateArgs()[i].getKind() == TemplateArgument::Declaration )  {
				if( ValueDecl *v = t->getTemplateArgs()[i].getAsDecl() ) {
					if( v->getAccess() == AS_protected   or  v->getAccess() == AS_private ) {
						outs() << "Private template VALUE arg: " << v->getNameAsString() << "\n";
						return false;
					}
				}
			}

		}
	}

	return r;
}

/// check if user requested binding for the given declaration
bool is_binding_requested(clang::CXXRecordDecl const *C, Config const &config)
{
	if( dyn_cast<ClassTemplateSpecializationDecl>(C) ) return false;
	bool bind = config.is_class_binding_requested( standard_name( C->getQualifiedNameAsString() ) ) or config.is_class_binding_requested( class_qualified_name(C) ) or config.is_namespace_binding_requested( namespace_from_named_decl(C) );
	for(auto & t : get_type_dependencies(C) ) bind &= !is_skipping_requested(t, config);
	return bind;
}

// check if user requested skipping for the given declaration
bool is_skipping_requested(clang::CXXRecordDecl const *C, Config const &config)
{
	bool skip = config.is_class_skipping_requested( standard_name( C->getQualifiedNameAsString() ) ) or config.is_class_skipping_requested( class_qualified_name(C) ) or config.is_namespace_skipping_requested( namespace_from_named_decl(C) );

	for(auto & t : get_type_dependencies(C) ) skip |= is_skipping_requested(t, config);

	return skip;
}


// Check if given method is const overload of some other method in class
bool is_const_overload(CXXMethodDecl *mc)
{
	string name = function_qualified_name(mc);

	string _const = " const";

	if( ends_with(name, _const) ) {
		name.resize( name.size() - _const.size() );

		CXXRecordDecl const *C = mc->getParent();

		for(auto m = C->method_begin(); m != C->method_end(); ++m) {
			if( m->getAccess() == mc->getAccess()  and  function_qualified_name(*m) == name ) return true;
		}
	}

	return false;;
}

// extract include needed for declaration and add it to includes
void add_relevant_includes(clang::CXXRecordDecl const *C, IncludeSet &includes, int level)
{
	if( !includes.add_decl(C, level) ) return;

	add_relevant_include_for_decl(C, includes);

	if( auto t = dyn_cast<ClassTemplateSpecializationDecl>(C) ) {

		for(uint i=0; i < t->getTemplateArgs().size(); ++i) {
			if( t->getTemplateArgs()[i].getKind() == TemplateArgument::Type ) {
				add_relevant_includes( t->getTemplateArgs()[i].getAsType().getDesugaredType(C->getASTContext()) , includes, level+1);
			}
			if( t->getTemplateArgs()[i].getKind() == TemplateArgument::Template ) {
				add_relevant_include_for_decl( t->getTemplateArgs()[i].getAsTemplate().getAsTemplateDecl()->getTemplatedDecl(), includes);
			}

			if( t->getTemplateArgs()[i].getKind() == TemplateArgument::Declaration ) {
				ValueDecl *v = t->getTemplateArgs()[i].getAsDecl();
				add_relevant_include_for_decl(v, includes);
			}
			if( t->getTemplateArgs()[i].getKind() == TemplateArgument::Integral ) {
				add_relevant_includes(t->getTemplateArgs()[i].getIntegralType(), includes, level+1);
			}
		}
	}

	// member function templates: we have to treat them separatly because instantiation types might not be included in class header
	for(auto d = C->decls_begin(); d != C->decls_end(); ++d) {
		if(FunctionTemplateDecl *ft = dyn_cast<FunctionTemplateDecl>(*d) ) {
			for(auto s = ft->spec_begin(); s != ft->spec_end(); ++s) {
				if(CXXMethodDecl *m = dyn_cast<CXXMethodDecl>(*s) ) {
					if( m->getAccess() == AS_public
						and  is_bindable(m)  //and  !is_skipping_requested(FunctionDecl const *F, Config const &config)
						and  !is_skipping_requested(m, Config::get())
						and  !isa<CXXConstructorDecl>(m)  and   !isa<CXXDestructorDecl>(m)
						and  !is_const_overload(m) ) {
						//m->dump();
						add_relevant_includes(m, includes, level+1);
					}
				}
			}

			// ArrayRef< TemplateArgument > args = ft->getInjectedTemplateArgs();
			// for(auto a = args.begin(); a != args.end(); ++a) {
			// 	outs() << "function template argument: " << template_argument_to_string( *a ) << " type:" << a->getKind() << "\n";
			// 	if( a->isDependent() ) {
			// 		outs() << "isDependent!" << "\n";
			// 	}
			// 	if( a->isInstantiationDependent() ) {
			// 		outs() << "isInstantiationDependent!" << "\n";
			// 	}
			// }

			// if( FunctionDecl const *td = ft->getTemplatedDecl() ) {
			// 	td->dump();
			// 	if( TemplateArgumentList const *ta = td->getTemplateSpecializationArgs() ) {
			// 		for(uint i=0; i < ta->size(); ++i) {
			// 			outs() << "function template argument: " << template_argument_to_string( ta->get(i) ) << "\n";
			// 			if( ta->get(i).isDependent() ) {
			// 				outs() << "isDependent!" << "\n";
			// 			}
			// 			if( ta->get(i).isInstantiationDependent() ) {
			// 				outs() << "isInstantiationDependent!" << "\n";
			// 			}
			// 		}
			// 	}
			// }

		}
	}

	//outs() << "isCompleteDefinition:" << C->isCompleteDefinition() << " id: " << C->getQualifiedNameAsString() << "\n";
	if( level < 2 ) {
		for(auto m = C->method_begin(); m != C->method_end(); ++m) {
			if( m->getAccess() == AS_public  and  is_bindable(*m)  /*and  !isa<CXXConstructorDecl>(*m)*/  and   !isa<CXXDestructorDecl>(*m) ) {
				add_relevant_includes(*m, includes, level+1);
			}
		}
	}
}


/// DEPRECATED Check if all bases have public default constructors
// bool is_default_default_constructor_available(CXXRecordDecl const *C)
// {
// 	for(auto b = C->bases_begin(); b!=C->bases_end(); ++b) {
// 		if( auto rt = dyn_cast<RecordType>(b->getType().getCanonicalType().getTypePtr() ) ) {
// 			if(CXXRecordDecl *R = cast<CXXRecordDecl>(rt->getDecl()) ) {
// 				if( !R->hasDefaultConstructor() ) return false;
// 				bool default_constructor_processed = false;
// 				for(auto t = R->ctor_begin(); t != R->ctor_end(); ++t) {
// 					if( t->isDefaultConstructor() ) {
// 						if( t->getAccess() == AS_private  or   !t->isUserProvided()  or  t->isDeleted() ) return false;
// 					}
// 					default_constructor_processed = true;
// 				}
// 				if( R->ctor_begin() != R->ctor_end()   and   !default_constructor_processed ) return false;
// 				if( !is_default_default_constructor_available(R) ) return false;
// 			}
// 		}
// 	}
// 	return true;
// }

// Check if all bases have public default constructors
bool base_default_default_constructor_available(CXXRecordDecl const *C)
{
	for(auto b = C->bases_begin(); b!=C->bases_end(); ++b) {
		if( auto rt = dyn_cast<RecordType>(b->getType().getCanonicalType().getTypePtr() ) ) {
			if(CXXRecordDecl *R = cast<CXXRecordDecl>(rt->getDecl()) ) {

				if( !R->hasDefaultConstructor() ) return false;

				for(auto t = R->ctor_begin(); t != R->ctor_end(); ++t) {
					if( t->isDefaultConstructor() ) {
						if( t->getAccess() == AS_private  /*or  !t->isUserProvided()*/  or  t->isDeleted() ) return false;
					}
				}

				if( !base_default_default_constructor_available(R) ) return false;
			}
		}
	}

	return true;
}



/// Generate string id that uniquly identify C++ binding object. For functions this is function prototype and for classes forward declaration.
string ClassBinder::id() const
{
	return class_qualified_name(C);
}


/// check if generator can create binding
bool ClassBinder::bindable() const
{
	return is_bindable(C) and !skipping_requested();
}


/// check if user requested binding for the given declaration
void ClassBinder::request_bindings_and_skipping(Config const &config)
{
	if( is_skipping_requested(C, config) ) Binder::request_skipping();
	else if( is_binding_requested(C, config) ) Binder::request_bindings();
}


/// extract include needed for this generator and add it to includes vector
void ClassBinder::add_relevant_includes(IncludeSet &includes) const
{
	for(auto & m : prefix_includes ) binder::add_relevant_includes(m, includes, 0);
	binder::add_relevant_includes(C, includes, 0);

	includes.add_include("<sstream> // __str__");
}

string binding_public_data_members(CXXRecordDecl const *C)
{
	string c;

	// binding protected data member that was made public in child class by 'using' declaration
	for(auto d = C->decls_begin(); d != C->decls_end(); ++d) {
		if(UsingDecl *u = dyn_cast<UsingDecl>(*d) ) {
			if( u->getAccess() == AS_public ) {
				for(auto s = u->shadow_begin(); s != u->shadow_end(); ++s) {
					if(UsingShadowDecl *us = dyn_cast<UsingShadowDecl>(*s) ) {
						if( FieldDecl *f = dyn_cast<FieldDecl>( us->getTargetDecl() ) ) {
							if( is_bindable(f) ) c += "\tcl" + bind_data_member(f, class_qualified_name(C)) + ";\n";
						}
					}
				}
			}
		}
	}

	for(auto d = C->decls_begin(); d != C->decls_end(); ++d) {
		if(FieldDecl *f = dyn_cast<FieldDecl>(*d) ) {
			if( f->getAccess() == AS_public  and  is_bindable(f) ) c += "\tcl" + bind_data_member(f, class_qualified_name(C)) + ";\n";
		}
	}
	return c;
}

// generate call-back structure name for given class
inline string callback_structure_name(CXXRecordDecl const *C)
{
	return "PyCallBack_" + python_class_name(C);
}


// Check if binding this class require creation of call-back structure to allow overriding virtual functions in Python
bool is_callback_structure_needed(CXXRecordDecl const *C)
{
	//C->dump();
	if( C->hasAttr<FinalAttr>() ) return false;

	for(auto m = C->method_begin(); m != C->method_end(); ++m) {
		if( m->getAccess() != AS_private  and  is_bindable(*m)   and  m->isVirtual()  and  !isa<CXXDestructorDecl>(*m) ) return true;
	}

	for(auto b = C->bases_begin(); b!=C->bases_end(); ++b) {
		if( b->getAccessSpecifier() != AS_private ) {
			if( auto rt = dyn_cast<RecordType>(b->getType().getCanonicalType().getTypePtr() ) ) {
				if(CXXRecordDecl *R = cast<CXXRecordDecl>(rt->getDecl()) ) {
					if( is_callback_structure_needed(R) ) return true;
				}
			}
		}
	}

	return false;
}

// Check if call-back structure that we can create will be constructible
bool is_callback_structure_constructible(CXXRecordDecl const *C)
{
	if( C->isAbstract() ) {
		for(auto m = C->method_begin(); m != C->method_end(); ++m) {
			if( m->isPure()  and  !isa<CXXConstructorDecl>(*m)  and  ( m->getAccess() == AS_private  or  !is_bindable(*m)  or  is_skipping_requested(*m, Config::get()) ) ) return false;
		}

		for(auto b = C->bases_begin(); b!=C->bases_end(); ++b) {
			if( auto rt = dyn_cast<RecordType>(b->getType().getCanonicalType().getTypePtr() ) ) {
				if(CXXRecordDecl *R = cast<CXXRecordDecl>(rt->getDecl()) ) {
				if( b->getAccessSpecifier() != AS_private ) {
					if( !is_callback_structure_constructible(R) ) return false;
				}
				else if( R->isAbstract() ) return false;
				}
			}
		}
	}

	return true;
}


// call_back_function_body_template is almost like PYBIND11_OVERLOAD_INT but specify pybind11::return_value_policy::reference
// #define PYBIND11_OVERLOAD_INT(ret_type, cname, name, ...) { \
//         pybind11::gil_scoped_acquire gil; \
//         pybind11::function overload = pybind11::get_overload(static_cast<const cname *>(this), name); \
//         if (overload) \
//             return overload(__VA_ARGS__).template cast<ret_type>();  }
//
// 		pybind11::gil_scoped_acquire gil; \
// 		pybind11::function overload = pybind11::get_overload(static_cast<const cname *>(this), name); \
// 		if (overload) { \
// 			auto o = overload(__VA_ARGS__); \
// 			if (pybind11::detail::cast_is_temporary_value_reference<ret_type>::value) { \
// 				static pybind11::detail::overload_caster_t<ret_type> caster; \
// 				return pybind11::detail::cast_ref<ret_type>(std::move(o), caster); \
// 			} \
// 			else return pybind11::detail::cast_safe<ret_type>(std::move(o)); \
// 		} \
// 	}


const char * call_back_function_body_template = R"_(
pybind11::gil_scoped_acquire gil;
pybind11::function overload = pybind11::get_overload(static_cast<const {0} *>(this), "{1}");
if (overload) {{
	auto o = overload.operator()<pybind11::return_value_policy::reference>({2});
	if (pybind11::detail::cast_is_temporary_value_reference<{3}>::value) {{
		static pybind11::detail::overload_caster_t<{3}> caster;
		return pybind11::detail::cast_ref<{3}>(std::move(o), caster);
	}}
	else return pybind11::detail::cast_safe<{3}>(std::move(o));
}}
)_";


// generate call-back overloads for all public virtual functions in C including it bases
string bind_member_functions_for_call_back(CXXRecordDecl const *C, string const & class_name, /*string const & base_type_alias,*/ set<string> &binded, int &ret_id, std::vector<clang::FunctionDecl const *> &prefix_includes/*, std::set<clang::NamedDecl const *> &prefix_includes_stack*/)
{
	string c;

	for(auto m = C->method_begin(); m != C->method_end(); ++m) {
		if( (m->getAccess() != AS_private)  and  is_bindable(*m)  and  is_overloadable(*m)
			and  !is_skipping_requested(*m, Config::get())
			and  !isa<CXXConstructorDecl>(*m)  and   !isa<CXXDestructorDecl>(*m)  and  m->isVirtual() and  !is_const_overload(*m)
			) {


			//if( m->hasAttr<NoExceptAttr>() ) {
			// 	(*m)->dump();
			// }

			string return_type = m->getReturnType().getCanonicalType().getAsString();  fix_boolean_types(return_type);
			tuple<string, string, string> args = function_arguments_for_py_overload(*m);

			string key = /*return_type + ' ' +*/ m->getNameAsString() + '(' + std::get<0>(args) + (m->isConst() ? ") const" : ")");

			if( !binded.count(key) ) {
				binded.insert(key);

				// m->hasAttr<OverrideAttr>, attribute list is automatically generated, see <release>/tools/clang/include/clang/Basic/AttrList.inc for attribute list
				if( m->hasAttr<FinalAttr>() ) continue; // we can not test this condition in a big if above because we need 'final' function to be marked as 'binded' so they will be be binded by some of the base classes

				if( return_type.find(',') != std::string::npos ) {
					string return_type_alias = "_binder_ret_" + std::to_string(ret_id); ++ret_id;
					c += "\tusing {} = {};\n"_format(return_type_alias, return_type);
					return_type = std::move(return_type_alias);
				}

				string python_name = python_function_name(*m);

				string exception_specification;
				if(FunctionProtoType const *fpt = dyn_cast<FunctionProtoType>( m->getType().getTypePtr() ) ) {
					if( fpt->getExceptionSpecType() & clang::ExceptionSpecificationType::EST_BasicNoexcept ) exception_specification = "noexcept ";
					if( fpt->getExceptionSpecType() & (clang::ExceptionSpecificationType::EST_DynamicNone | clang::ExceptionSpecificationType::EST_Dynamic | clang::ExceptionSpecificationType::EST_MSAny) ) exception_specification = "throw() ";
				}

				c += "\t{} {}({}){} {}override {{ "_format(return_type, m->getNameAsString(), std::get<0>(args), m->isConst() ? " const" : "", exception_specification);

				c += indent( fmt::format(call_back_function_body_template, class_name, /*class_qualified_name(C), */python_name, std::get<1>(args), return_type), "\t\t");
				if( m->isPure() ) c+= "\t\tpybind11::pybind11_fail(\"Tried to call pure virtual function \\\"{}::{}\\\"\");\n"_format(C->getNameAsString(), python_name);
				else c+= "\t\treturn {}::{}({});\n"_format(C->getNameAsString(), m->getNameAsString(), std::get<1>(args));
				c += "\t}\n";

				// c += string(m->isPure() ? "PYBIND11_OVERLOAD_PURE_NAME" : "PYBIND11_OVERLOAD_NAME") + '(';
				// c += return_type + ", " + base_type_alias + ", \"" + python_function_name(*m) + "\", " + C->getNameAsString() + "::" + m->getNameAsString();
				// c += ( args.second.size() ? ", " + args.second : "");
				// c += "); }\n";

				// Note: does not work when argument is std::shared_ptr passed by-reference...
				// string python_name = python_function_name(*m);
				// c += "PYBIND11_OVERLOAD_INT({}, \"{}\"{}); "_format(return_type, python_name, std::get<2>(args).size() ? ", " + std::get<2>(args) : "");
				// if( m->isPure() ) c+= "pybind11::pybind11_fail(\"Tried to call pure virtual function \\\"{}::{}\\\"\");"_format(C->getNameAsString(), python_name);
				// else c+= "return {}::{}({});"_format(C->getNameAsString(), m->getNameAsString(), std::get<1>(args));
				// c += " }\n";

				prefix_includes.push_back(*m);
				//add_relevant_includes(*m, prefix_includes, prefix_includes_stack, 0);
			}
		}
	}

	for(auto b = C->bases_begin(); b!=C->bases_end(); ++b) {
		if( b->getAccessSpecifier() != AS_private ) {
			if( auto rt = dyn_cast<RecordType>(b->getType().getCanonicalType().getTypePtr() ) ) {
				if(CXXRecordDecl *R = cast<CXXRecordDecl>(rt->getDecl()) ) {
					c += bind_member_functions_for_call_back(R, class_name, /*base_type_alias,*/ binded, ret_id, prefix_includes);
					//add_relevant_includes(R, prefix_includes, prefix_includes_stack, 0);
				}
			}
		}
	}

	return c;
}

// Genarate code for defining 'call-back struct' that act as 'trampoline' and allows overlading virtual functions in Python
void ClassBinder::generate_prefix_code()
{
	if( !is_callback_structure_needed(C) ) return;

	prefix_code_  = generate_comment_for_declaration(C);
	prefix_code_ += "struct {0} : public {1} {{\n\tusing {1}::{2};\n\n"_format(callback_structure_name(C), class_qualified_name(C), C->getNameAsString());

	// string base_type_alias = "_binder_base_";
	// prefix_code_ += "\tusing {} = {};\n\n"_format(base_type_alias, class_qualified_name(C));

	set<string> binded;  int ret_id = 0;
	prefix_code_ += bind_member_functions_for_call_back(C, class_qualified_name(C), /*base_type_alias,*/ binded, ret_id, prefix_includes);
	prefix_code_ += "};\n\n";
}

string binding_public_member_functions(CXXRecordDecl const *C, bool callback_structure, bool callback_structure_constructible, Context &context)
{
	string c;
	// binding protected member functions that was made public in child class by 'using' declaration
	for(auto d = C->decls_begin(); d != C->decls_end(); ++d) {
		if(UsingDecl *u = dyn_cast<UsingDecl>(*d) ) {
			if( u->getAccess() == AS_public ) {
				for(auto s = u->shadow_begin(); s != u->shadow_end(); ++s) {
					if(UsingShadowDecl *us = dyn_cast<UsingShadowDecl>(*s) ) {
						if( CXXMethodDecl *m = dyn_cast<CXXMethodDecl>( us->getTargetDecl() ) ) {
							if( is_bindable(m)  and  !is_skipping_requested(m, Config::get())  and  !isa<CXXConstructorDecl>(m)  and   !isa<CXXDestructorDecl>(m)  and  !is_const_overload(m) ) {
								// Create a new CXXRecordDecl and insert base method into inherited class so bind_function correctly resolve parent namespace for function as 'child::' instead of 'base::'
								// CXXRecordDecl NC(*C);
								// CXXMethodDecl *nm = CXXMethodDecl::Create(m->getParentASTContext(), &NC, m->getLocStart(), m->getNameInfo(), m->getType(), m->getTypeSourceInfo(),
								// 										  m->getStorageClass(), m->isInlineSpecified(), m->isConstexpr() , m->getLocStart());
								// it looks like LLVM will delete this object when parent CXXRecordDecl is destroyed so commenting out for now... // delete nm;
								c += bind_function("\tcl", m, context, C);
							}
						}
					}
				}
			}
		}

		if(FunctionTemplateDecl *ft = dyn_cast<FunctionTemplateDecl>(*d) ) {
			for(auto s = ft->spec_begin(); s != ft->spec_end(); ++s) {
				if(CXXMethodDecl *m = dyn_cast<CXXMethodDecl>(*s) ) {
					if( m->getAccess() == AS_public
						and  is_bindable(m)  //and  !is_skipping_requested(FunctionDecl const *F, Config const &config)
						and  !is_skipping_requested(m, Config::get())
						and  !isa<CXXConstructorDecl>(m)  and   !isa<CXXDestructorDecl>(m)
						and  !is_const_overload(m) ) {
						//m->dump();

						c += bind_function("\tcl", m, context);
					}
				}
			}
		}
	}

	for(auto m = C->method_begin(); m != C->method_end(); ++m) {
		if( m->getAccess() == AS_public
			and  is_bindable(*m)  //and  !is_skipping_requested(FunctionDecl const *F, Config const &config)
			and  !is_skipping_requested(*m, Config::get())
			and  !isa<CXXConstructorDecl>(*m)  and   !isa<CXXDestructorDecl>(*m)
			and  !is_const_overload(*m) ) {
			//(*m)->dump();

			c += bind_function("\tcl", *m, context);
		}
	}

	return c;
}


string binding_template_bases(CXXRecordDecl const *C, bool callback_structure, bool callback_structure_constructible, Context &context)
{
	string c;

	if( dyn_cast<ClassTemplateSpecializationDecl>(C) ) { // for template classes explicitly bind data members and member functions from public base classes
		for(auto b = C->bases_begin(); b!=C->bases_end(); ++b) {
			if( b->getAccessSpecifier() == AS_public) {
				if( auto rt = dyn_cast<RecordType>(b->getType().getCanonicalType().getTypePtr() ) ) {
					if(CXXRecordDecl *R = cast<CXXRecordDecl>(rt->getDecl()) ) {
						if( !is_skipping_requested(R, Config::get()) ) {
							c += binding_public_data_members(R);
							c += binding_public_member_functions(R, callback_structure, callback_structure_constructible, context);
							c += binding_template_bases(R, callback_structure, callback_structure_constructible, context);
						}
					}
				}
			}
		}
	}

	return c;
}

/// Create forward-binding for given class which consist of only class type without any member, function or constructors
string bind_forward_declaration(CXXRecordDecl const *C, Context &context)
{
	string const qualified_name{ class_qualified_name(C) };
	string const module_variable_name = context.module_variable_name( namespace_from_named_decl(C) );
	//string const decl_namespace = namespace_from_named_decl(C);

	string const include = relevant_include(C);

	string c = "\t// Forward declaration for: " + qualified_name + " file:" + (include.size() ? include.substr(1, include.size()-2) : "") + " line:" + line_number(C) + "\n";

	string maybe_holder_type =  ", std::shared_ptr<{}>"_format(qualified_name);
	if( is_inherited_from_enable_shared_from_this(C) ) maybe_holder_type = ", std::shared_ptr<{}>"_format(qualified_name);
	else if( CXXDestructorDecl * d = C->getDestructor() ) {
		if( d->getAccess() != AS_public ) maybe_holder_type = ", " + qualified_name + '*';
	}

	c += '\t' + R"(pybind11::class_<{}{}>({}, "{}");)"_format(qualified_name, maybe_holder_type, module_variable_name, python_class_name(C)) + "\n\n";

	return c;
}


/// check if any of the base classes is wrappable and if generate a string describing them: , pybind11::base<BaseClass>()
string ClassBinder::maybe_base_classes(Context &context)
{
	string r;

	static std::vector<string> const skip_list = {"std::enable_shared_from_this", "std::string", "std::basic_string", "std::pair", "std::tuple"};

	for(auto b = C->bases_begin(); b!=C->bases_end(); ++b) {
		if( b->getAccessSpecifier() == AS_public) {
			if( auto rt = dyn_cast<RecordType>(b->getType().getCanonicalType().getTypePtr() ) ) {
				if(CXXRecordDecl *R = cast<CXXRecordDecl>(rt->getDecl()) ) {
					auto e = std::find(skip_list.begin(), skip_list.end(), standard_name( R->getQualifiedNameAsString() ));

					if( e == skip_list.end()  and  is_bindable(R)  and  !is_skipping_requested(R, Config::get()) ) {
						//r = ", pybind11::base<{}>()"_format( class_qualified_name(R) );
						r += ", {}"_format( class_qualified_name(R) );

						binder::request_bindings(b->getType().getCanonicalType(), context);

						dependencies_.push_back(R);
					}
				}
			}
		}
	}

	return r;
}

/// generate binding code for this object by using external user-provided binder
void ClassBinder::bind_with(string const &binder, Context &context)
{
	string c = '\t' + generate_comment_for_declaration(C);

	string const module_variable_name =  context.module_variable_name( namespace_from_named_decl(C) );

	c += '\t' + binder + standard_name( template_specialization(C) ) + '(' + module_variable_name;

	if( auto t = dyn_cast<ClassTemplateSpecializationDecl>(C) ) {
		for(uint i=0; i < t->getTemplateArgs().size(); ++i) {

			string templ = mangle_type_name( standard_name( template_argument_to_string(t->getTemplateArgs()[i]) ), true);
			fix_boolean_types(templ);
			c += ", \"" + templ + '"';
		}
	}

	c += ");\n\n";

	code() = c;
}


char const * constructor_template = "\tcl.def(\"__init__\", []({2} *self_{0}) {{ new (self_) {2}({1}); }}, \"doc\");";
char const * constructor_if_template = "\tcl.def(\"__init__\", [cl_type](pybind11::handle self_{0}) {{ if (self_.get_type() == cl_type) new (self_.cast<{2} *>()) {2}({1}); else new (self_.cast<{3} *>()) {3}({1}); }}, \"doc\");";

// Generate binding for given function: .def("foo", (std::string (aaaa::A::*)(int) ) &aaaa::A::foo, "doc")
// constructor_types is pair<Base, Alias> - if one of these is absent empty string is expected
string bind_constructor(CXXConstructorDecl const *T, pair<string, string> const &constructor_types, uint args_to_bind, bool request_bindings_f, Context &context)
{
	//string function_name = python_function_name(F);
	//string function_qualified_name { F->getQualifiedNameAsString() };

	string c;
	if( args_to_bind == T->getNumParams() ) {
		c = "\tcl.def(pybind11::init<{}>()"_format( function_arguments(T) );

		for(uint i=0; i<T->getNumParams()  and  i < args_to_bind; ++i) {
			c += ", pybind11::arg(\"{}\")"_format( string( T->getParamDecl(i)->getName() ) );

			if(request_bindings_f) request_bindings( T->getParamDecl(i)->getOriginalType(), context);
		}
		c += ");\n";
	}
	else {
		pair<string, string> args = function_arguments_for_lambda(T, args_to_bind);

		string params = args_to_bind ? ", " + args.first : "";

		if( constructor_types.first.size()  and  constructor_types.second.size()  ) c = fmt::format(constructor_if_template, params, args.second, constructor_types.first, constructor_types.second);
		else if( constructor_types.first.size() ) c = fmt::format(constructor_template, params, args.second, constructor_types.first);
		else c = fmt::format(constructor_template, params, args.second, constructor_types.second);
	}

	return c;
}

// Generate binding for given constructor. If constructor have default arguments generate set of bindings by creating separate bindings for each argument with default.
string bind_constructor(CXXConstructorDecl const *T, pair<string, string> const &constructor_types, Context &context)
{
	string code;

	uint args_to_bind = 0;
	for(; args_to_bind < T->getNumParams(); ++args_to_bind) {
		if( T->getParamDecl(args_to_bind)->hasDefaultArg() ) break;
	}

	for(; args_to_bind <= T->getNumParams(); ++args_to_bind) code += bind_constructor(T, constructor_types, args_to_bind, args_to_bind == T->getNumParams(), context) + '\n';

	return code;
}


/// generate (if any) bindings for Python __str__ by using appropriate global operator<<
std::string ClassBinder::bind_repr(Context &context)
{
	string c;

	if( FunctionDecl const * F = context.global_insertion_operator(C) ) {
		string qualified_name = class_qualified_name(C);
		//outs() << "Found insertion operator for: " << class_qualified_name(C) << "\n";

		c += "\n\tcl.def(\"__str__\", []({} const &o) -> std::string {{ std::ostringstream s; s << o; return s.str(); }} );\n"_format(qualified_name);

		prefix_includes.push_back(F);
	}

	return c;
}


/// generate binding code for this object and all its dependencies
void ClassBinder::bind(Context &context)
{
	if( is_binded() ) return;

	string const qualified_name_without_template = standard_name( C->getQualifiedNameAsString() );
	std::map<string, string> const &external_binders = Config::get().binders();
	if( external_binders.count(qualified_name_without_template) ) {
		bind_with( external_binders.at(qualified_name_without_template), context );
		return;
	}

	assert( bindable() && "Attempt to bind non-bindable CXXRecord!");

	bool callback_structure = is_callback_structure_needed(C);
	bool callback_structure_constructible = callback_structure and is_callback_structure_constructible(C);
	if(callback_structure  and  callback_structure_constructible) generate_prefix_code();

	string const qualified_name{ class_qualified_name(C) };
	string const module_variable_name = context.module_variable_name( namespace_from_named_decl(C) );
	//string const decl_namespace = namespace_from_named_decl(C);

	string c = "{ " + generate_comment_for_declaration(C);

	string const binding_qualified_name = callback_structure_constructible ? callback_structure_name(C) : qualified_name;

	string maybe_holder_type =  ", std::shared_ptr<{}>"_format(qualified_name); // for now enable std::shared_ptr by default
	if( is_inherited_from_enable_shared_from_this(C) ) maybe_holder_type = ", std::shared_ptr<{}>"_format(qualified_name);
	else if( CXXDestructorDecl * d = C->getDestructor() ) {
		if( d->getAccess() != AS_public ) maybe_holder_type = ", " + qualified_name + '*';
	}

	string maybe_trampoline = callback_structure_constructible ? ", " + binding_qualified_name : "";

	c += '\t' + R"(pybind11::class_<{}{}{}{}> cl({}, "{}", "{}");)"_format(qualified_name, maybe_holder_type, maybe_trampoline, maybe_base_classes(context), module_variable_name, python_class_name(C), generate_documentation_string_for_declaration(C)) + '\n';
	c += "\tpybind11::handle cl_type = cl;\n\n";

	//if( C->isAbstract()  and  callback_structure) c += "\tcl.def(pybind11::init<>());\n";

	if( !C->isAbstract()  or  callback_structure_constructible) {
		bool default_constructor_processed = false;

		string constructors;
		for(auto t = C->ctor_begin(); t != C->ctor_end(); ++t) {
			if( t->getAccess() == AS_public  and  !t->isMoveConstructor()  and  is_bindable(*t)  /*and  t->doesThisDeclarationHaveABody()*/ ) {
				if( t->isCopyConstructor()  and  callback_structure_constructible) constructors += "\tcl.def(pybind11::init<{} const &>());\n"_format(binding_qualified_name);
				else constructors += bind_constructor(*t, std::make_pair(!C->isAbstract() ? qualified_name : "",
																		 callback_structure_constructible ? callback_structure_name(C) : ""), context);
			}
			if( t->isDefaultConstructor() ) default_constructor_processed = true;
		}
		//c += "\t// hasDefaultConstructor={} needsImplicitDefaultConstructor={} base_default_default_constructor_available={}\n"_format(C->hasDefaultConstructor(), C->needsImplicitDefaultConstructor(), base_default_default_constructor_available(C));

		if( !default_constructor_processed  and  C->needsImplicitDefaultConstructor()  and  base_default_default_constructor_available(C)
			// ( C->ctor_begin() == C->ctor_end() and  is_default_default_constructor_available(C)  and  !default_constructor_processed)
			// or (C->hasDefaultConstructor()  and  !default_constructor_processed )

			// if( C->hasDefaultConstructor()  or ( !default_constructor_processed and  is_default_default_constructor_available(C)) /*and  !C->needsImplicitDefaultConstructor() and !C->hasNonTrivialDefaultConstructor()*/

			/*and  !C->needsImplicitDefaultConstructor() and !C->hasNonTrivialDefaultConstructor()*/
			) {  // No constructors defined, adding default constructor

			c += "\tcl.def(pybind11::init<>());\n";  // making sure that default is appering first
		}
		c += constructors;
	}

	// binding public enums
	for(auto d = C->decls_begin(); d != C->decls_end(); ++d) {
		if(EnumDecl *e = dyn_cast<EnumDecl>(*d) ) {
			if( e->getAccess() == AS_public ) {
				//outs() << "Enum: " << e->getQualifiedNameAsString() << "\n";
				c += bind_enum("cl", e);
			}
		}
	}

	c += binding_public_data_members(C);
	c += binding_public_member_functions(C, callback_structure,  callback_structure_constructible, context);

	c += binding_template_bases(C, callback_structure,  callback_structure_constructible, context);

	c += bind_repr(context);

	std::map<string, string> const &external_add_on_binders = Config::get().add_on_binders();
	if( external_add_on_binders.count(qualified_name_without_template) ) c += "\n\t{}(cl);\n"_format(external_add_on_binders.at(qualified_name_without_template));

	c += "}\n";

	code() = indent(c, "\t");
}

} // namespace binder
