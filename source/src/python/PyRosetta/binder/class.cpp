// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
	return C->getNameAsString() + template_specialization(C);
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
	return C->getQualifiedNameAsString() + template_specialization(C);
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
	//outs() << "is_inherited_from_enable_shared_from_this: " << C->getQualifiedNameAsString() << " " << C->isCompleteDefinition() << "\n";
	if( C->getQualifiedNameAsString() == "std::enable_shared_from_this" ) return true;
	if( C->isCompleteDefinition() ) {
		for(auto b = C->bases_begin(); b!=C->bases_end(); ++b) {
			//if( b->getAccessSpecifier() == AS_public)
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
				if( m->getAccess() == AS_public  and  m->isCopyAssignmentOperator()  and  !m->doesThisDeclarationHaveABody()  ) return false;
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


			//errs() << "is_bindable(CXXRecordDecl): " << C->getQualifiedNameAsString() << "  -  " << const_cast<ClassTemplateSpecializationDecl*>(ts)->getMostRecentDecl()->isCompleteDefinition() << "\n";
		}
		else return false;
	}

	//outs() << "is_bindable(CXXRecordDecl): " << C->getQualifiedNameAsString() << template_specialization(C) << " " << r << "\n";

	if( auto t = dyn_cast<ClassTemplateSpecializationDecl>(C) ) {
		//C->dump();
		for(uint i=0; i < t->getTemplateArgs().size(); ++i) {

			if( t->getTemplateArgs()[i].getKind() == TemplateArgument::Type ) {
				if( !is_bindable( t->getTemplateArgs()[i].getAsType() ) ) return false;
			}

			//if( t->getTemplateArgs()[i].getKind() == TemplateArgument::Integral ) return true;

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
	bool bind = config.is_class_binding_requested( C->getQualifiedNameAsString() ) or config.is_class_binding_requested( class_qualified_name(C) ) or config.is_namespace_binding_requested( namespace_from_named_decl(C) );
	for(auto & t : get_type_dependencies(C) ) bind &= !is_skipping_requested(t, config);
	return bind;
}

// check if user requested skipping for the given declaration
bool is_skipping_requested(clang::CXXRecordDecl const *C, Config const &config)
{
	bool skip = config.is_class_skipping_requested( C->getQualifiedNameAsString() ) or config.is_class_skipping_requested( class_qualified_name(C) ) or config.is_namespace_skipping_requested( namespace_from_named_decl(C) );

	for(auto & t : get_type_dependencies(C) ) skip |= is_skipping_requested(t, config);

	return skip;
}



// extract include needed for declaration and add it to includes
void add_relevant_includes(clang::CXXRecordDecl const *C, vector<string> &includes, set<NamedDecl const *> &stack, int level)
{
	if( stack.count(C)  /*or  stack.size() > 16*/ ) return; else stack.insert(C);

	//outs() << stack.size() << " ";
	// if( begins_with(C->getQualifiedNameAsString(), "boost::property")
	// 	) outs() << "add_relevant_includes(class): " << C->getQualifiedNameAsString() << template_specialization(C) << "\n";

	add_relevant_include_for_decl(C, includes);

	if( auto t = dyn_cast<ClassTemplateSpecializationDecl>(C) ) {

		for(uint i=0; i < t->getTemplateArgs().size(); ++i) {
			if( t->getTemplateArgs()[i].getKind() == TemplateArgument::Type ) {
				add_relevant_includes( t->getTemplateArgs()[i].getAsType().getDesugaredType(C->getASTContext()) , includes, stack, level+1);
			}
			if( t->getTemplateArgs()[i].getKind() == TemplateArgument::Template ) {
				add_relevant_include_for_decl( t->getTemplateArgs()[i].getAsTemplate().getAsTemplateDecl()->getTemplatedDecl(), includes);
			}

			if( t->getTemplateArgs()[i].getKind() == TemplateArgument::Declaration ) {
				ValueDecl *v = t->getTemplateArgs()[i].getAsDecl();
				add_relevant_include_for_decl(v, includes);
			}
			if( t->getTemplateArgs()[i].getKind() == TemplateArgument::Integral ) {
				add_relevant_includes(t->getTemplateArgs()[i].getIntegralType(), includes, stack, level+1);
			}
		}
	}

	//outs() << "isCompleteDefinition:" << C->isCompleteDefinition() << " id: " << C->getQualifiedNameAsString() << "\n";
	if( level < 2 ) {
		for(auto m = C->method_begin(); m != C->method_end(); ++m) {
			if( m->getAccess() == AS_public  and  is_bindable(*m)  /*and  !isa<CXXConstructorDecl>(*m)*/  and   !isa<CXXDestructorDecl>(*m) ) {
				add_relevant_includes(*m, includes, stack, level+1);
			}
		}
	}
}

/// Check if all bases have public default constructors
bool is_default_default_constructor_available(CXXRecordDecl const *C)
{
	for(auto b = C->bases_begin(); b!=C->bases_end(); ++b) {
		if( auto rt = dyn_cast<RecordType>(b->getType().getCanonicalType().getTypePtr() ) ) {
			if(CXXRecordDecl *R = cast<CXXRecordDecl>(rt->getDecl()) ) {
				if( !R->hasDefaultConstructor() ) return false;
				for(auto t = R->ctor_begin(); t != R->ctor_end(); ++t) {
					if( t->isDefaultConstructor() ) {
						if( t->getAccess() == AS_private  or   !t->isUserProvided() ) return false;
					}
				}
				if( !is_default_default_constructor_available(R) ) return false;
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
void  ClassBinder::request_bindings_and_skipping(Config const &config)
{
	if( is_skipping_requested(C, config) ) Binder::request_skipping();
	else if( is_binding_requested(C, config) ) Binder::request_bindings();
}


/// extract include needed for this generator and add it to includes vector
void ClassBinder::add_relevant_includes(std::vector<std::string> &includes, std::set<clang::NamedDecl const *> &stack) const
{
	binder::add_relevant_includes(C, includes, stack, 0);
}

string binding_public_data_members(CXXRecordDecl const *C)
{
	string c;
	for(auto d = C->decls_begin(); d != C->decls_end(); ++d) {
		if(FieldDecl *f = dyn_cast<FieldDecl>(*d) ) {
			if( f->getAccess() == AS_public  and  is_bindable(f) ) c += "\tcl" + bind_data_member(f, class_qualified_name(C)) + ";\n";
		}
	}
	return c;
}

string binding_public_member_functions(CXXRecordDecl const *C, /*std::set<clang::NamedDecl const *> const &members_to_skip, */Context &context)
{
	string c;
	for(auto m = C->method_begin(); m != C->method_end(); ++m) {
		if( m->getAccess() == AS_public  and  is_bindable(*m)  //and  !is_skipping_requested(FunctionDecl const *F, Config const &config)
			//and  !members_to_skip.count(*m)
			and  !is_skipping_requested(*m, Config::get())
			and  !isa<CXXConstructorDecl>(*m)  and   !isa<CXXDestructorDecl>(*m) ) {
			//(*m)->dump();
			c += bind_function("\tcl", *m, context);
		}
	}

	return c;
}

string binding_template_bases(CXXRecordDecl const *C, Context &context)
{
	// TODO: add somekind of redundancy detection here

	string c;

	if( dyn_cast<ClassTemplateSpecializationDecl>(C) ) { // for template classes explicitly bind data members and member functions from public base classes
		for(auto b = C->bases_begin(); b!=C->bases_end(); ++b) {
			if( b->getAccessSpecifier() == AS_public) {
				if( auto rt = dyn_cast<RecordType>(b->getType().getCanonicalType().getTypePtr() ) ) {
					if(CXXRecordDecl *R = cast<CXXRecordDecl>(rt->getDecl()) ) {
						c += binding_public_data_members(R);
						c += binding_public_member_functions(R, /*members_to_skip,*/ context);
						c += binding_template_bases(R, context);
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
					auto e = std::find(skip_list.begin(), skip_list.end(), R->getQualifiedNameAsString());

					if( e == skip_list.end()  and  is_bindable(R)  and  !is_skipping_requested(R, Config::get()) ) {
						r = ", pybind11::base<{}>()"_format( class_qualified_name(R) );
						binder::request_bindings(b->getType().getCanonicalType(), context);
						dependencies_.push_back(R);
						break; // right now pybind11 support only one base class
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
	string const qualified_name{ class_qualified_name(C) };
	string const include = relevant_include(C);

	string c = "// " + qualified_name + " file:" + (include.size() ? include.substr(1, include.size()-2) : "") + " line:" + line_number(C) + "\n";

	string const module_variable_name =  context.module_variable_name( namespace_from_named_decl(C) );

	c += binder + template_specialization(C) + '(' + module_variable_name;

	if( auto t = dyn_cast<ClassTemplateSpecializationDecl>(C) ) {
		for(uint i=0; i < t->getTemplateArgs().size(); ++i) {

			string templ = mangle_type_name( template_argument_to_string(t->getTemplateArgs()[i]), true);
			fix_boolean_types(templ);
			c += ", \"" + templ + '"';
		}
	}

	c += ");\n\n";

	code() = indent(c, "\t");
}

/// generate binding code for this object and all its dependencies
void ClassBinder::bind(Context &context)
{
	if( is_binded() ) return;

	string const qualified_name_without_template = C->getQualifiedNameAsString();
	std::map<string, string> const &external_binders = Config::get().binders();
	if( external_binders.count(qualified_name_without_template) ) {
		bind_with( external_binders.at(qualified_name_without_template), context );
		return;
	}

	assert( bindable() && "Attempt to bind non-bindable CXXRecord!");

	string const qualified_name{ class_qualified_name(C) };
	string const module_variable_name = context.module_variable_name( namespace_from_named_decl(C) );
	//string const decl_namespace = namespace_from_named_decl(C);

	string const include = relevant_include(C);

	string c = "{ // " + qualified_name + " file:" + (include.size() ? include.substr(1, include.size()-2) : "") + " line:" + line_number(C) + "\n";
	//if( decl_namespace.size() ) c += "\tusing namespace " + decl_namespace + ";\n";

	string maybe_holder_type =  ", std::shared_ptr<{}>"_format(qualified_name);
	if( is_inherited_from_enable_shared_from_this(C) ) maybe_holder_type = ", std::shared_ptr<{}>"_format(qualified_name);
	else if( CXXDestructorDecl * d = C->getDestructor() ) {
		if( d->getAccess() != AS_public ) maybe_holder_type = ", " + qualified_name + '*';
	}

	c += '\t' + R"(pybind11::class_<{}{}> cl({}, "{}"{});)"_format(qualified_name, maybe_holder_type, module_variable_name, python_class_name(C), maybe_base_classes(context)) + '\n';

	if( !C->isAbstract() ) {
		//c += "// hasDefaultConstructor:{}\n"_format(C->hasDefaultConstructor() );
		bool default_constructor_processed = false;

		string constructors;
		for(auto t = C->ctor_begin(); t != C->ctor_end(); ++t) {
			if( t->getAccess() == AS_public  and  !t->isMoveConstructor()  and  is_bindable(*t)  /*and  t->doesThisDeclarationHaveABody()*/ ) {
				constructors += "\tcl.def(pybind11::init<{}>());\n"_format( function_arguments(*t) );

			}
			if( t->isDefaultConstructor() ) default_constructor_processed = true;
		}

		if( /*C->ctor_begin() == C->ctor_end()  and*/  C->hasDefaultConstructor()  and  !default_constructor_processed  and  is_default_default_constructor_available(C) /*and  !C->needsImplicitDefaultConstructor() and !C->hasNonTrivialDefaultConstructor()*/ ) {  // No constructors defined, adding default constructor
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
	c += binding_public_member_functions(C, /*members_to_skip,*/ context);

	c += binding_template_bases(C, context);

	//outs() << "typename_from_type_decl: " << typename_from_type_decl(C) << "\n";

	c += "}\n";

	code() = indent(c, "\t");
}




} // namespace binder
