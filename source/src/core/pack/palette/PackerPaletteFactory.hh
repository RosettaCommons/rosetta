// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/palette/PackerPaletteFactory.hh
/// @brief  Headers for factory class for creating instances of PackerPalettes (e.g. for RosettaScripts).
/// @author Vikram K. Mulligan (vmullig@flatironinstitute.org).

#ifndef INCLUDED_core_pack_palette_PackerPaletteFactory_hh
#define INCLUDED_core_pack_palette_PackerPaletteFactory_hh

// Unit Headers
#include <core/pack/palette/PackerPaletteFactory.fwd.hh>

// Package Headers
#include <core/pack/palette/PackerPalette.fwd.hh>
#include <core/pack/palette/PackerPaletteCreator.fwd.hh>

// Basic headers
#include <basic/datacache/DataMap.fwd.hh>

// Utility Headers
#include <utility/SingletonBase.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <utility/thread/threadsafe_creation.hh>

// c++ headers
#include <string>
#include <map>

#include <utility/vector0.hh>

#ifdef MULTI_THREADED
// C++11 Headers
#ifdef OLDER_GXX_STDLIB
#include <atomic>
#endif //OLDER_GXX_STDLIB
#include <mutex>
#endif //MULTI_THREADED

namespace core {
namespace pack {
namespace palette {

// singleton class
class PackerPaletteFactory : public utility::SingletonBase< PackerPaletteFactory >
{
public:
	friend class utility::SingletonBase< PackerPaletteFactory >;

	typedef utility::vector1< PackerPaletteOP > PackerPaletteOPs;
	typedef std::map< std::string, PackerPaletteCreatorOP > PackerPaletteCreatorMap;
	typedef utility::tag::Tag Tag;
	typedef utility::tag::TagOP TagOP;
	typedef utility::tag::TagCOP TagCOP;

public:

	void factory_register( PackerPaletteCreatorOP );

	/// @brief add a prototype, using its default type name as the map key
	///
	void add_creator( PackerPaletteCreatorOP );

	/// @brief Does a type exist?
	///
	bool has_type( std::string const & ) const;

	/// @brief Return new PackerPalette by key lookup in packer_palette_creator_map_.
	/// @details New PackerPalette parses Tag if provided)
	PackerPaletteOP newPackerPalette(
		std::string const &,
		basic::datacache::DataMap & datamap,
		TagCOP = TagCOP( TagOP( new Tag() ) )
	) const;

	/// @brief Fills packer_palette_creator_map_ vector with new PackerPalettes from nested "PACKER_PALETTES" TagCOP
	///
	void newPackerPalettes( PackerPaletteOPs &, basic::datacache::DataMap & datamap, TagCOP ) const;

	/// @brief Fills packer_palette_creator_map_ vector with new PackerPalettes from xml-like tag file.
	///
	void newPackerPalettes( PackerPaletteOPs &, basic::datacache::DataMap & datamap, std::string const & ) const;

	static std::string packer_palette_xml_schema_group_name();

	/// @brief The PackerPaletteFactory is the point of entry for the definition of the XML Schemas
	/// for every PackerPalette that may be instantiated from an XML file.  It is responsible for defining
	/// an xs:group named "packer_palette" listing each of the packer-palette-complex types that may
	/// be initialized using the PackerPaletteFactory.  It will iterate across each of the
	/// PackerPaletteCreators that it contains, asking them for the XML schema of the PackerPalette
	/// that each one is responsible for creating.
	/// @details By convention, the name assigned to each of the complexTypes for PackerPalettes should be
	/// what is returned by the function "complex_type_name_for_packer_palette" (declared in
	/// core/pack/palette/xsd_util.hh) when given the argument returned by that PackerPalette's
	/// PackerPaletteCreator's keyname() function. So long as the writing of XML schema for your packer
	/// palette is accomplished by calling the functions in core/pack/palette/xsd_util.hh, then
	/// this should happen automatically.
	void define_packer_palette_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;

	/// @brief Get the XML schema for a given packer palette.
	/// @details Throws an error if the packer palette is unknown to Rosetta.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void provide_xml_schema( std::string const &packer_palette_name, utility::tag::XMLSchemaDefinition & xsd ) const;

	/// @brief Create a packer palette based on global defaults, and return an owning pointer to it.
	/// @details By default, makes a DefaultPackePalette.  If the user provides options for additional residues,
	/// makes a CustomBaseTypePackerPalette.
	PackerPaletteOP create_packer_palette_from_global_defaults() const;

private:
	/// @brief Private constructor.
	///
	PackerPaletteFactory();

	/// @brief Private destructor.
	///
	virtual ~PackerPaletteFactory();

	/// @brief Create the PackerPalette that we clone whenever create_packer_palette_from_global_defaults() is called.
	/// @details This function is intended ONLY to be called in a threadsafe manner from create_packer_palette_from_global_defaults().
	/// @note This function reads from the global options system, and can read from disk.
	static PackerPaletteOP initialize_packer_palette_from_global_defaults();

private:

	/// @brief Map of string->PackerPaletteOP.
	///
	PackerPaletteCreatorMap packer_palette_creator_map_;

	/// @brief The default PackerPalette, handed out by create_packer_palette_from_global_defaults().
	/// @details Lazily created.
	mutable PackerPaletteOP default_palette_;

#ifdef MULTI_THREADED
	mutable std::mutex default_palette_mutex_;
#ifdef OLDER_GXX_STDLIB
	mutable std::atomic_bool default_palette_bool_;
#endif //OLDER_GXX_STDLIB
#endif //MULTI_THREADED

};

} //namespace palette
} //namespace pack
} //namespace core

#endif
