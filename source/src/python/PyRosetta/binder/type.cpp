// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   binder/type.cpp
/// @brief  Various functionality for handline clang::QualType's
/// @author Sergey Lyskov

#include <type.hpp>

#include <class.hpp>
#include <enum.hpp>
#include <util.hpp>

#include <clang/AST/ASTContext.h>
#include <clang/AST/ExprCXX.h>

//#include <experimental/filesystem>
#include <cstdlib>
#include <fstream>

using namespace llvm;
using namespace clang;
using std::string;
using std::vector;


namespace binder {


/// check if user requested binding for the given QualType
bool is_binding_requested(clang::QualType const &qt, Config const &config)
{
	if( PointerType const *pt = dyn_cast<PointerType>( qt.getTypePtr() ) ) return is_binding_requested(pt->getPointeeType(), config);

	if( ReferenceType const *rt = dyn_cast<ReferenceType>( qt.getTypePtr() ) ) return is_binding_requested(rt->getPointeeType(), config);

	if( CXXRecordDecl *r = qt->getAsCXXRecordDecl() ) return is_binding_requested(r, config);

	return false;
}


// check if user requested skipping for the given declaration
bool is_skipping_requested(QualType const &qt, Config const &config)
{
	if( PointerType const *pt = dyn_cast<PointerType>( qt.getTypePtr() ) ) return is_skipping_requested(pt->getPointeeType(), config);

	if( ReferenceType const *rt = dyn_cast<ReferenceType>( qt.getTypePtr() ) ) return is_skipping_requested(rt->getPointeeType(), config);

	if( CXXRecordDecl *r = qt->getAsCXXRecordDecl() ) return is_skipping_requested(r, config);

	return false;
}


// extract include path needed for declaration itself (without template dependency if any), return empty string if no include could be found
string relevant_include(NamedDecl const *decl)
{
	ASTContext & ast_context( decl->getASTContext() );
	SourceManager & sm( ast_context.getSourceManager() );

	// outs() << "  Source location: " << F->getLocation().printToString(sm) << "\n";
	//outs() << "  Source location: " << sm.getFilename( decl->getLocation() )  << "\n";

	FileID fid = sm.getFileID( decl->getLocation() );
	SourceLocation include = sm.getIncludeLoc(fid);

	//outs() << "  Include location: " << include.printToString( ast_context.getSourceManager() ) << "\n";
	//outs() << "  printToString: " << include.printToString( ast_context.getSourceManager() ) << "\n";
	//include.dump(ast_context.getSourceManager());
	//outs() << "  Line: " << sm.getCharacterData(include) << "\n";
	// outs() << "  Source location: " << sm.getFileEntryForID( FullSourceLoc(F->getLocation(), sm ).getFileID() )->getName() << "\n";

	string include_string;
	if( include.isValid() ) {
		char const *data = sm.getCharacterData(include);

		if( strlen(data) > 2 ) {  // should be at least 3 chars: open/close symbol + file name
			char terminator = *data == '"' ? '"' : '>';

			include_string.push_back(*data); ++data;
			for(; *data and  *data != terminator; ++data ) include_string.push_back(*data);
			if(*data == terminator) include_string.push_back(*data);

			//outs() << "  Include for " << decl->getNameAsString() << ": " << include_string << "\n";
		}
		if( include_string.size()  and  include_string[0]=='"' ) include_string.resize(0); // avoid adding include in quotes because compiler will not be able to find them
	}

	return include_string;
}


// extract include path needed for declaration itself (without template dependency if any), do nothing if include could not be found (ie for build-in's)
void add_relevant_include_for_decl(NamedDecl const *decl, vector<string> &includes/*, std::set<clang::NamedDecl const *> &stack*/)
{
	using std::make_pair;
	static vector< std::pair<string, string> > const name_map = {
		make_pair("std::allocator",        "<memory>"),
		make_pair("std::basic_string",     "<string>"),
		make_pair("std::char_traits",      "<string>"),
		make_pair("std::iterator",         "<iterator>"),
		make_pair("std::list",             "<list>"),
		make_pair("std::locale",           "<locale>"),
		make_pair("std::map",              "<map>"),
		make_pair("std::pair",             "<utility>"),
		make_pair("std::reverse_iterator", "<iterator>"),
		make_pair("std::set",              "<set>"),
		make_pair("std::vector",           "<vector>"),

		make_pair("__gnu_cxx::__normal_iterator", "<iterator>"),
	};

	static vector< std::pair<string, string> > const include_map = {
		make_pair("<bits/postypes.h>", "<ios>"),
		make_pair("<bits/shared_ptr.h>", "<memory>"),
		make_pair("<bits/stl_function.h>", "<functional>"),
		make_pair("<bits/stl_tree.h>", "<map>"),
		make_pair("<bits/ios_base.h>", "<ios>"),
		make_pair("<bits/istream.tcc>", "<istream>"),
		make_pair("<bits/ostream.tcc>", "<ostream>"),
		make_pair("<bits/streambuf.tcc>", "<streambuf>"),
		make_pair("<bits/stl_bvector.h>", "<vector>"),
		make_pair("<bits/algorithmfwd.h>", "<algorithm>"),
		make_pair("<bits/stl_algo.h>", "<algorithm>"),
		make_pair("<bits/fstream.tcc>", "<fstream>"),
		make_pair("<bits/sstream.tcc>", "<sstream>"),
		make_pair("<bits/stl_list.h>", "<list>"),
		make_pair("<bits/stl_deque.h>", "<deque>"),
	};


	string name = decl->getQualifiedNameAsString();

	for(auto & p : name_map) {
		if( begins_with(name, p.first) ) {
			includes.push_back(p.second);
			//#ifndef NDEBUG
			//includes.back() += " // " + name;
			//#endif
			return;
		}
	}


	string include = relevant_include(decl);

	for(auto & i : include_map) {
		if( include == i.first ) { include = i.second; break; }
	}

	//if( ( namespace_from_named_decl(decl) == "std" or namespace_from_named_decl(decl) == "__gnu_cxx" )
	//	and  begins_with(include, "<bits/") ) return;
	//if( begins_with(name, "std::")  and  begins_with(include, "<bits/") ) return;

	if( include.size() ) {
		//#ifndef NDEBUG
		//include += " // " + name;
		//#endif
		includes.push_back(include);
	}
	// else {
	// 	for(;;) {
	// 		if(DeclContext const * d = decl->getDeclContext() ) {
	// 			if( auto r = dyn_cast<CXXRecordDecl>(d) ) {
	// 				outs() << "include from decl context for: " << decl->getQualifiedNameAsString() << "\n";
	// 				add_relevant_includes(r, includes, stack);
	// 			}
	// 		}
	// 		else break;
	// 	}
	// }
}


/// extract include needed for this generator and add it to includes vector
void add_relevant_includes(QualType const &qt, /*const ASTContext &context,*/ std::vector<std::string> &includes, std::set<clang::NamedDecl const *> &stack, int level)
{
	//QualType qt = qt.getDesugaredType(context);
	//outs() << "add_relevant_includes(qt): " << qt.getAsString() << "\n";
	if( PointerType const *pt = dyn_cast<PointerType>( qt.getTypePtr() ) ) add_relevant_includes(pt->getPointeeType(), includes, stack, level);
	if( ReferenceType const *rt = dyn_cast<ReferenceType>( qt.getTypePtr() ) ) add_relevant_includes(rt->getPointeeType(), includes, stack, level);
	if( CXXRecordDecl *r = qt->getAsCXXRecordDecl() ) add_relevant_includes(r, includes, stack, level);
	if( EnumDecl *e = dyn_cast_or_null<EnumDecl>( qt->getAsTagDecl() ) ) add_relevant_includes(e, includes, stack, level);
}


// check if given QualType is bindable
bool is_bindable(QualType const &qt)
{
	// if( qt->isVoidPointerType() ) return false; // for now refuse to bind void* types - note the recent pybind11 can handle bindings for void*

	bool r = true;

	r &= !qt->isFunctionPointerType()  and  !qt->isRValueReferenceType()  and  !qt->isInstantiationDependentType()  and  !qt->isArrayType();  //and  !qt->isConstantArrayType()  and  !qt->isIncompleteArrayType()  and  !qt->isVariableArrayType()  and  !qt->isDependentSizedArrayType()

	if( PointerType const *pt = dyn_cast<PointerType>( qt.getTypePtr() ) ) {
		if( pt->getPointeeType()->isPointerType() ) return false;  // refuse to bind 'value**...' types
		if( pt->getPointeeType()->isArithmeticType() ) return false;  // refuse to bind 'int*, doublle*...' types
		if( pt->getPointeeType()->isArrayType() or pt->getPointeeType()->isConstantArrayType() ) return false;  // refuse to bind 'T* v[]...' types

		//qt->dump();
		r &= is_bindable( pt->getPointeeType()/*.getCanonicalType()*/ );
	}

	if( ReferenceType const *rt = dyn_cast<ReferenceType>( qt.getTypePtr() ) ) {
		QualType pqt = rt->getPointeeType();
		//outs() << "#### " << pqt.getAsString() << "\n";

		// special handling for std::pair&  and  std::tuple&  whitch pybind11 can't pass by refference
		string pqt_name = pqt.getAsString();
		if( begins_with(pqt_name, "std::pair")  or  begins_with(pqt_name, "std::tuple") ) return false;  // but we allow bindings for 'const std::tuple' and 'const std::pair'

		//rt->dump();
		//outs() << "Ref " << qt.getAsString() << " -> " << is_bindable( rt->getPointeeType().getCanonicalType() ) << "\n";
		r &= is_bindable( pqt/*.getCanonicalType()*/ );
	}

	if( Type const *tp = qt/*.getCanonicalType()*/.getTypePtrOrNull() ) {
		if( CXXRecordDecl *rd = tp->getAsCXXRecordDecl() ) {
			//outs() << "is_bindable qt CXXRecordDecl:" << rd->getQualifiedNameAsString() << " " << is_bindable(rd) << "\n";
			r &= is_bindable(rd);
		}
		if( TagDecl *td = tp->getAsTagDecl() ) {
			if( td->getAccess() == AS_protected  or  td->getAccess() == AS_private  ) return false;
		}
	}


	// if( auto rt = dyn_cast<RecordType>( qt.getTypePtr() ) ) {
	// 	if( auto cxxr = dyn_cast<CXXRecordDecl>( rt->getDecl() ) ) r &= is_bindable(cxxr);
	// }

	//outs() << " isIncompleteArrayType(): " << qt->isIncompleteArrayType() << " r:" << r << "\n";
	//qt->dump();

	//outs() << "Qt " << qt.getAsString() << " -> " << r /*<< " isInstantiationDependentType: " << qt->isInstantiationDependentType() */ << "\n";
	return r;
}



/// extract type info from QualType if any and bind relative type
void request_bindings(clang::QualType const &qt, Context &context)
{
	if( TagDecl *td = qt->getAsTagDecl() ) {
		//if( TypeDecl * type_decl = dyn_cast<TypeDecl>( b->named_decl() ) ) types[ typename_from_type_decl(type_decl) ] = b;
		if( td->isCompleteDefinition() ) context.request_bindings( typename_from_type_decl(td) );
		//else outs() << "void bind(clang::QualType, Context &): isCompleteDefinition is false for: "<< qt.getAsString() << "\n";
		//context.request_bindings( typename_from_type_decl(td) );
	}
}


/// check if given class/struct is builtin in Python and therefor should not be binded
bool is_python_builtin(NamedDecl const *C)
{
	//outs() << "Considering: " << C->getQualifiedNameAsString() << "\n";
	string name = C->getQualifiedNameAsString();
	//if( begins_with(name, "class ") ) name = name.substr(6); // len("class ")

	std::vector<string> known_builtin = {"std::pair", "std::tuple", "std::vector", "std::map", "std::list", "std::set", "std::basic_string",
										 "std::allocator", "std::initializer_list", "std::shared_ptr", "std::weak_ptr", "std::less",
										 "std::iterator", "std::reverse_iterator",
										 "std::_Rb_tree_iterator", "std::_Rb_tree_const_iterator", "__gnu_cxx::__normal_iterator",
 										 "std::_List_iterator", "std::_List_const_iterator",
 	};

	for(auto &k : known_builtin) {
		if( begins_with(name, k) ) return true;
	}

	return false;
}

/*

/// check if user requested binding for the given QualType
bool f_is_binding_requested::cxx_record_decl(clang::CXXRecordDecl const *C)
{
	bool bind = config.is_namespace_binding_requested( namespace_from_named_decl(C) );

	for(auto & t : get_type_dependencies(C) ) bind |= apply(t);

	return bind;
}



/// check if user requested skipping for the given QualType
bool f_is_skipping_requested::cxx_record_decl(clang::CXXRecordDecl const *C)
{
	bool skip = config.is_class_skipping_requested( C->getQualifiedNameAsString() ) or config.is_class_skipping_requested( class_qualified_name(C) ) or config.is_namespace_skipping_requested( namespace_from_named_decl(C) );

	for(auto & t : get_type_dependencies(C) ) skip |= apply(t);

	return skip;
}

  void f_add_relevant_includes::cxx_record_decl(clang::CXXRecordDecl const *C)
{
	if( stack.count(C) ) return; else stack.insert(C);

	add_relevant_include_for_decl(C, includes);

	for(auto & t : get_type_dependencies(C) ) apply(t);

	for(auto & d : get_decl_dependencies(C) ) add_relevant_include_for_decl(d, includes);
}

void f_add_relevant_includes::enum_decl(clang::EnumDecl const *E)
{
	add_relevant_include_for_decl(E, includes);
}

*/




} // namespace binder
