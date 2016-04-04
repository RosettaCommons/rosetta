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


// generate vector<NamedDecl const *> with all declarations related to class
// vector<NamedDecl const *> get_decl_dependencies(CXXRecordDecl const *C)
// {
// 	vector<NamedDecl const *> r;

// 	if( auto t = dyn_cast<ClassTemplateSpecializationDecl>(C) ) {

// 		for(uint i=0; i < t->getTemplateArgs().size(); ++i) {
// 			if( t->getTemplateArgs()[i].getKind() == TemplateArgument::Template ) {
// 				r.push_back( t->getTemplateArgs()[i].getAsTemplate().getAsTemplateDecl()->getTemplatedDecl() );
// 			}
// 		}
// 	}

// 	for(auto m = C->method_begin(); m != C->method_end(); ++m) {
// 		if( m->getAccess() == AS_public  and  is_bindable(*m)  /*and  !isa<CXXConstructorDecl>(*m)*/  and   !isa<CXXDestructorDecl>(*m) ) {
// 			r.push_back(*m);
// 		}
// 	}

// 	return r;
// }


// Return true if class have direct or inderect std::enable_shared_from_this as base class
bool is_inherited_from_enable_shared_from_this(CXXRecordDecl const *C)
{
	// non recursive 1-level deep version
	// if( C->isCompleteDefinition() ) {
	// 	for(auto b = C->bases_begin(); b!=C->bases_end(); ++b) {
	// 		//if( b->getAccessSpecifier() == AS_public)
	// 		if( auto r = dyn_cast<RecordType>(b->getType().getCanonicalType().getTypePtr() ) ) {
	// 			CXXRecordDecl *rd = cast<CXXRecordDecl>(r->getDecl());
	// 		if( rd  and  rd->getQualifiedNameAsString() == "std::enable_shared_from_this") return true;
	// 		}
	// 	}
	// }
	// return false;

	// recursive version, but it appears that we only need to loop up for direct bases only so using simplified version above for now
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


/// generate list of class/enum names on which this CXXRecordDecl depend to get binded ommiting build-in types
// vector<string> calculate_dependency(CXXRecordDecl *C)
// {
// }


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

	if( C->isDependentType() ) return false;
	if( !C->isCompleteDefinition()  and  !dyn_cast<ClassTemplateSpecializationDecl>(C) ) return false;
	if( C->getAccess() == AS_protected  or  C->getAccess() == AS_private ) return false;


	//r &= C->isCompleteDefinition() /* and C->getDefinition() */  /*and  C->hasDefinition()*/;

	//outs() << "is_bindable(CXXRecordDecl): " << C->getQualifiedNameAsString() << template_specialization(C) << " " << r << "\n";

	if( auto t = dyn_cast<ClassTemplateSpecializationDecl>(C) ) {
		//C->dump();
		for(uint i=0; i < t->getTemplateArgs().size(); ++i) {
			// if( template_specialization(C) == "<1,utility::options::OptionCollection::OptionTypes,std::allocator<utility::options::OptionCollection::OptionTypes>>" ) {
			// 	if( t->getTemplateArgs()[i].getKind() == TemplateArgument::Type ) {
			// 		Type const *tp = t->getTemplateArgs()[i].getAsType().getTypePtrOrNull();
			// 		if(tp) tp->dump();
			// 		if( TagDecl *td = tp->getAsTagDecl() ) {
			// 			if( td->getAccess() != AS_public ) {
			// 				outs() << "Access NOT public!\n";
			// 			//outs() << "Private template TYPE arg: " << td->getNameAsString() << "\n";
			// 				return false;
			// 			}
			// 		} else {
			// 			outs() << "Access public!\n";
			// 		}
			// 	}
			// 	//t->getTemplateArgs()[i].dump( outs() );
			// }

			if( t->getTemplateArgs()[i].getKind() == TemplateArgument::Type ) {
				if( !is_bindable( t->getTemplateArgs()[i].getAsType() ) ) return false;
				// Type const *tp = t->getTemplateArgs()[i].getAsType().getTypePtrOrNull();
				// if( TagDecl *td = tp->getAsTagDecl() ) {
				// 	if( td->getAccess() == AS_protected  or  td->getAccess() == AS_private  ) {
				// 		//outs() << "Private template TYPE arg: " << td->getNameAsString() << "\n";
				// 		return false;
				// 	}
				// }


				// if( tp  and  (tp->isRecordType() or tp->isEnumeralType())  and  !tp->isBuiltinType()  and  !begins_wtih(template_argument_to_string(t->getTemplateArgs()[i]), "std::") ) {
				// 	TagDecl *td = tp->getAsTagDecl();

				// 	//if(td)
				// 	if(FieldDecl *fd = dyn_cast<FieldDecl>(td) ) {
				// 		outs() << "FieldDecl!!! " << fd-> getNameAsString() << "\n";
				// 		if( fd->getAccess() != AS_public ) return false;
				// 		//if( f->getAccess() == AS_public  and  is_bindable(f) ) c+= '\t' + bind_data_member(f, class_qualified_name(C)) + '\n';
				// 	}
				// 	//CXXRecordDecl *rd =	tp->getAsCXXRecordDecl();
				// }
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
	// todo: bindging for abstract classes
	//if(r) r &= !C->isAbstract();  // need an 'if' here or clang assert got triggered on classed with incomplete definitions

	// if(r) {
	// 	outs() << C->getQualifiedNameAsString() << " isCXXClassMember:" << C->isCXXClassMember() << " isCompleteDefinition:" << C->isCompleteDefinition() //<< " isBeingDefined:" << C->isBeingDefined()
	// 		   << " isDependentType:" << C->isDependentType() << " isCXXInstanceMember:" << C->isCXXInstanceMember() << " isExternallyVisible:" << C->isExternallyVisible()
	// 		   << " getNumTemplateParameterLists:" << C->getNumTemplateParameterLists()  << " isa<ClassTemplateSpecializationDecl>:" << isa<ClassTemplateSpecializationDecl>(C) << "\n";
	// 	C->dump();
	// 	//auto l = C->getTemplateParameterList();
	// 	//for(auto p = l.begin(); p!=l.end(); ++p) outs() << (*p)->getNameAsString() << "\n";
	// }

	return r;
}

/// check if user requested binding for the given declaration
bool is_binding_requested(clang::CXXRecordDecl const *C, Config const &config)
{
	if( dyn_cast<ClassTemplateSpecializationDecl>(C) ) return false;
	bool bind = config.is_class_binding_requested( C->getQualifiedNameAsString() ) or config.is_class_binding_requested( class_qualified_name(C) ) or config.is_namespace_binding_requested( namespace_from_named_decl(C) );
	for(auto & t : get_type_dependencies(C) ) bind &= !is_skipping_requested(t, config);
	return bind;

	// crashing on boost boost::mpl::is_na
	// bool bind = config.is_namespace_binding_requested( namespace_from_named_decl(C) );
	// for(auto & t : get_type_dependencies(C) ) bind |= is_binding_requested(t, config);
	// return bind;

	// bool bind = config.is_namespace_binding_requested( namespace_from_named_decl(C) );
	// //outs() << "Skipping: " << b.named_decl()->getQualifiedNameAsString() << "\n";
	// if(bind) {
	// 	if( auto t = dyn_cast<ClassTemplateSpecializationDecl>(C) ) {
	// 		for(uint i=0; i < t->getTemplateArgs().size(); ++i) {
	// 			if( t->getTemplateArgs()[i].getKind() == TemplateArgument::Type ) {
	// 				Type const *tp = t->getTemplateArgs()[i].getAsType().getTypePtrOrNull();
	// 				if( tp  and  (tp->isRecordType() or tp->isEnumeralType()) and  !tp->isBuiltinType() ) {
	// 					if(CXXRecordDecl *rd = tp->getAsCXXRecordDecl() ) bind &= is_binding_requested(rd, config);
	// 				}
	// 			}
	// 		}
	// 	}
	// }
	// return bind;
}

// check if user requested skipping for the given declaration
bool is_skipping_requested(clang::CXXRecordDecl const *C, Config const &config)
{
	bool skip = config.is_class_skipping_requested( C->getQualifiedNameAsString() ) or config.is_class_skipping_requested( class_qualified_name(C) ) or config.is_namespace_skipping_requested( namespace_from_named_decl(C) );

	for(auto & t : get_type_dependencies(C) ) skip |= is_skipping_requested(t, config);

	return skip;

	// bool skip = config.is_class_skipping_requested( C->getQualifiedNameAsString() ) or config.is_class_skipping_requested( class_qualified_name(C) ) or config.is_namespace_skipping_requested( namespace_from_named_decl(C) );
	// if(!skip) {
	// 	if( auto t = dyn_cast<ClassTemplateSpecializationDecl>(C) ) {
	// 		for(uint i=0; i < t->getTemplateArgs().size(); ++i) {
	// 			if( t->getTemplateArgs()[i].getKind() == TemplateArgument::Type ) {
	// 				skip |= is_skipping_requested(t->getTemplateArgs()[i].getAsType(), config);
	// 				// Type const *tp = t->getTemplateArgs()[i].getAsType().getTypePtrOrNull();
	// 				// if( tp  and  (tp->isRecordType() or tp->isEnumeralType()) and  !tp->isBuiltinType() ) {
	// 				// 	if(CXXRecordDecl *rd = tp->getAsCXXRecordDecl() ) skip |= is_skipping_requested(rd, config);
	// 				// }
	// 			}
	// 		}
	// 	}
	// }
	// //outs() << "is_skipping_requested: " << C->getQualifiedNameAsString() << " - " << skip << "\n";
	// return skip;
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

				// Type const *tp = t->getTemplateArgs()[i].getAsType().getTypePtrOrNull();
				// if( tp  and  (tp->isRecordType() or tp->isEnumeralType()) and  !tp->isBuiltinType() ) {
				// 	CXXRecordDecl *rd = tp->getAsCXXRecordDecl();
				// 	//TagDecl *td = tp->getAsTagDecl();
				// 	add_relevant_includes(rd, includes);
				// }
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
	// if( C->isCXXInstanceMember() ) {
	// 	if( auto r = dyn_cast<CXXRecordDecl>(C->getDeclContext()) ) add_relevant_includes(r, includes, stack);
	// 	if( auto e = dyn_cast<EnumDecl>(C->getDeclContext()) ) add_relevant_includes(e, includes, stack);
	// }
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
	for(auto m = C->method_begin(); m != C->method_end(); ++m) {
		if( m->getAccess() == AS_public  and  is_bindable(*m)  /*and  !isa<CXXConstructorDecl>(*m)*/  and   !isa<CXXDestructorDecl>(*m) ) {
			if( is_skipping_requested(*m, config) ) members_to_skip.insert(*m);
		}
	}

	if( is_skipping_requested(C, config) ) Binder::request_skipping();
	else if( is_binding_requested(C, config) ) Binder::request_bindings();
}


/// extract include needed for this generator and add it to includes vector
void ClassBinder::add_relevant_includes(std::vector<std::string> &includes, std::set<clang::NamedDecl const *> &stack) const
{
	binder::add_relevant_includes(C, includes, stack, 0);
}


/// generate binding code for this object and all its dependencies
void ClassBinder::bind(Context &context)
{
	if( is_binded() ) return;

	assert( bindable() && "Attempt to bind non-bindable CXXRecord!");

	string const indentation="\t";
	string const module_variable_name =  context.module_variable_name( namespace_from_named_decl(C) );

	string const qualified_name{ class_qualified_name(C) };
	string const include = relevant_include(C);

	string c = "{ // " + qualified_name + " file:" + include.substr(1, include.size()-2) + " line:" + line_number(C) + "\n";
	c += "\tusing namespace " + namespace_from_named_decl(C) + ";\n\t";

	// class_<A>(module_a, "A") or class_<A, std::shared_ptr<A>>(module_a, "A")
	//string maybe_holder_type;
	// if( is_inherited_from_enable_shared_from_this(C) ) maybe_holder_type = ", std::shared_ptr<{}>"_format(qualified_name);
	// else if( CXXDestructorDecl * d = C->getDestructor() ) {
	// 	if( d->getAccess() != AS_public ) maybe_holder_type = ", " + qualified_name + '*';
	// }

	string maybe_holder_type =  ", std::shared_ptr<{}>"_format(qualified_name);
	if( is_inherited_from_enable_shared_from_this(C) ) maybe_holder_type = ", std::shared_ptr<{}>"_format(qualified_name);
	else if( CXXDestructorDecl * d = C->getDestructor() ) {
		if( d->getAccess() != AS_public ) maybe_holder_type = ", " + qualified_name + '*';
	}

	c += R"(pybind11::class_<{}{}> cl({}, "{}");)"_format(qualified_name, maybe_holder_type, module_variable_name, class_name(C)) + '\n';

	if( !C->isAbstract() ) {
		if( C->ctor_begin() == C->ctor_end() ) {  // No constructors defined, adding default constructor
			c+= "\tcl.def(pybind11::init<>());\n\n";
		}
		else {
			bool added=false;
			for(auto t = C->ctor_begin(); t != C->ctor_end(); ++t) {
				if( t->getAccess() == AS_public  and  !t->isMoveConstructor()  and  is_bindable(*t) ) { added=true;  c+= "\tcl.def(pybind11::init<{}>());\n"_format( function_arguments(*t) ); }
			}
			if(added) c += '\n';
		}
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

	// binding public data members
	for(auto d = C->decls_begin(); d != C->decls_end(); ++d) {
		if(FieldDecl *f = dyn_cast<FieldDecl>(*d) ) {
			if( f->getAccess() == AS_public  and  is_bindable(f) ) c+= "\tcl" + bind_data_member(f, class_qualified_name(C)) + ";\n";
		}
	}

	// binding public member functions
	for(auto m = C->method_begin(); m != C->method_end(); ++m) {
		if( m->getAccess() == AS_public  and  is_bindable(*m)  //and  !is_skipping_requested(FunctionDecl const *F, Config const &config)
			and  !members_to_skip.count(*m)
			and  !isa<CXXConstructorDecl>(*m)  and   !isa<CXXDestructorDecl>(*m) ) {
			//(*m)->dump();
			c += "\tcl" + bind_function(*m, context) + ";\n";
		}
	}

	//outs() << "typename_from_type_decl: " << typename_from_type_decl(C) << "\n";

	c += "}\n\n";

	code() = indent(c, indentation);
}




} // namespace binder
