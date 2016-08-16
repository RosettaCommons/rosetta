// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/ResidueSelectorFactory.hh
/// @brief  Class for instantiating arbitrary ResidueSelectors from a string --> ResidueSelectorCreator map
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_select_residue_selector_ResidueSelectorFactory_HH
#define INCLUDED_core_select_residue_selector_ResidueSelectorFactory_HH

// Package headers
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/select/residue_selector/ResidueSelectorCreator.fwd.hh>

// Basic headers
#include <basic/datacache/DataMap.fwd.hh>

// Utility headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <map>
#include <string>

#ifdef MULTI_THREADED
#ifdef CXX11
// C++11 Headers
#include <atomic>
#include <mutex>
#endif
#endif

namespace core {
namespace select {
namespace residue_selector {

class ResidueSelectorFactory {
public:
	static ResidueSelectorFactory * get_instance();

	void factory_register( ResidueSelectorCreatorOP creator );
	bool has_type( std::string const & ) const;

	ResidueSelectorOP new_residue_selector(
		std::string const & selector_name,
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	) const;

	/// @brief Should the Factory throw an exception or call utility::exit when it encounters the
	/// second of two ResidueSelctorCreators with the same keyname?  It's default behavior is to
	/// call utility::exit, but this method allows you to set it so that it will throw an
	/// exception instead (which is unit testable).
	void set_throw_on_double_registration();

	/// @brief The %ResidueSelectorFactory is the point of entry for the definition of the XML Schemas
	/// for every ResidueSelector that may be instantiated from a file.  It is  responsible for defining
	/// an xs:group named "residue_selector" listing each of the residue-selector-complex types that may
	/// be initialized using the %ResidueSelectorFactory and to iterate across each of the
	/// ResidueSelectorCreators it contains asking them for the XML schema of the ResidueSelector they
	/// are responsible for creating.
	void define_residue_selector_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;

	static std::string residue_selector_xml_schema_group_name();

private:
	ResidueSelectorFactory();

	// Unimplemented -- uncopyable
	ResidueSelectorFactory( ResidueSelectorFactory const & );
	ResidueSelectorFactory const & operator = ( ResidueSelectorFactory const & );

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static ResidueSelectorFactory * create_singleton_instance();

#ifdef MULTI_THREADED
#ifdef CXX11
public:

	/// @brief This public method is meant to be used only by the
	/// utility::thread::safely_create_singleton function and not meant
	/// for any other purpose.  Do not use.
	static std::mutex & singleton_mutex();

private:
	static std::mutex singleton_mutex_;
#endif
#endif

private:
#if defined MULTI_THREADED && defined CXX11
	static std::atomic< ResidueSelectorFactory * > instance_;
#else
	static ResidueSelectorFactory * instance_;
#endif

private:
	std::map< std::string, ResidueSelectorCreatorOP > creator_map_;
	bool throw_on_double_registration_;
};


} //namespace residue_selector
} //namespace select
} //namespace core


#endif
