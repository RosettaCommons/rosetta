// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/palette/PackerPaletteFactory.cc
/// @brief  Factory class for creating instances of PackerPalettes (e.g. for RosettaScripts).
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).

#include <core/pack/palette/PackerPaletteFactory.hh>

#include <core/pack/palette/DefaultPackerPalette.hh>
#include <core/pack/palette/NCAADefaultPackerPalette.hh>
#include <core/pack/palette/CustomBaseTypePackerPalette.hh>
#include <core/pack/palette/PackerPalette.hh>
#include <core/pack/palette/PackerPaletteCreator.hh>
#include <core/pack/palette/xsd_util.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <utility/exit.hh> // runtime_assert, utility_exit_with_message
#include <utility/io/izstream.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/xml_schema_group_initialization.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>

#include <utility/vector0.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

namespace core {
namespace pack {
namespace palette {

static basic::Tracer TR( "core.pack.palette.PackerPaletteFactory" );

PackerPaletteFactory::~PackerPaletteFactory(){}

/// @brief the default PackerPalettes are now initialized in core/init/init.cc via the registrator/creator scheme
PackerPaletteFactory::PackerPaletteFactory() :
	packer_palette_creator_map_(),
	default_palette_(nullptr)
#ifdef MULTI_THREADED
	,
	default_palette_mutex_()
#ifdef OLDER_GXX_STDLIB
	,
	default_palette_bool_(false)
#endif //OLDER_GXX_STDLIB
#endif //MULTI_THREADED
{}

void
PackerPaletteFactory::factory_register( PackerPaletteCreatorOP creator )
{
	if ( packer_palette_creator_map_.find( creator->keyname() ) != packer_palette_creator_map_.end() ) {
		utility_exit_with_message( "Factory Name Conflict: Two or more PackerPaletteCreators registered with the name " + creator->keyname() );
	}
	add_creator( creator );
}

/// @brief add a PackerPalette prototype creator
void
PackerPaletteFactory::add_creator( PackerPaletteCreatorOP creator )
{
	runtime_assert( creator != 0 );
	packer_palette_creator_map_[ creator->keyname() ] = creator;
}

bool PackerPaletteFactory::has_type( std::string const & type ) const
{
	return ( packer_palette_creator_map_.find( type ) != packer_palette_creator_map_.end() );
}


PackerPaletteOP
PackerPaletteFactory::newPackerPalette(
	std::string const & type,
	basic::datacache::DataMap & datamap,
	TagCOP tag /* = boost::shared_ptr< Tag >() */
) const
{
	PackerPaletteCreatorMap::const_iterator iter( packer_palette_creator_map_.find( type ) );
	if ( iter != packer_palette_creator_map_.end() ) {
		PackerPaletteOP packer_palette( iter->second->create_packer_palette() );
		// parse tag if tag pointer is pointing to one
		if ( tag.get() != nullptr ) packer_palette->parse_my_tag( tag, datamap );
		return packer_palette;
	} else {
		TR<<"Available options: ";
		for ( PackerPaletteCreatorMap::const_iterator to_iter = packer_palette_creator_map_.begin(); to_iter != packer_palette_creator_map_.end(); ++to_iter ) {
			TR<<to_iter->first<<", ";
		}
		TR<<std::endl;
		utility_exit_with_message( type + " is not known to the PackerPaletteFactory. Was its PackerPaletteCreator class registered at initialization?" );
		return nullptr;
	}
}

/// @brief recurse tag file to find PACKER_PALETTES definitions
void
PackerPaletteFactory::newPackerPalettes( PackerPaletteOPs & ppops, basic::datacache::DataMap & datamap, TagCOP tag ) const
{
	typedef utility::vector0< TagCOP > TagCOPs;
	TR.Trace << "Tag name " << tag->getName();
	if ( tag->getTags().empty() ) { TR.Trace << " (empty)" << std::endl; return; }
	else TR.Trace << std::endl;
	TagCOPs const subtags( tag->getTags() );
	if ( tag->getName() == "PACKER_PALETTES" ) {
		for ( TagCOPs::const_iterator tp( subtags.begin() ), tp_e( subtags.end() ); tp != tp_e; ++tp ) {
			std::string const type( (*tp)->getName() );
			PackerPaletteOP new_pp = newPackerPalette( type, datamap, *tp );
			runtime_assert( new_pp != 0 );
			ppops.push_back( new_pp );
			TR << "Created and parsed anonymous PackerPalette of type " << type << "." << std::endl;
		}
	}
	// recurse
	for ( TagCOPs::const_iterator tp( subtags.begin() ), tp_e( subtags.end() ); tp != tp_e; ++tp ) {
		newPackerPalettes( ppops, datamap, *tp );
	}
}

std::string
PackerPaletteFactory::packer_palette_xml_schema_group_name() {
	return "packer_palette";
}

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
void
PackerPaletteFactory::define_packer_palette_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	try {
		utility::tag::define_xml_schema_group(
			packer_palette_creator_map_,
			packer_palette_xml_schema_group_name(),
			& core::pack::palette::complex_type_name_for_packer_palette,
			xsd );
	} catch ( utility::excn::Exception const & e ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Could not generate an XML Schema for PackerPalettes from PackerPalettesFactory; offending class"
			" must call core::pack::palette::complex_type_name_for_packer_palette when defining"
			" its XML Schema\n" + e.msg() );
	}
}

/// @brief Get the XML schema for a given packer palette.
/// @details Throws an error if the packer palette is unknown to Rosetta.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
PackerPaletteFactory::provide_xml_schema(
	std::string const &packer_palette_name,
	utility::tag::XMLSchemaDefinition & xsd
) const {
	if ( ! has_type( packer_palette_name ) ) {
		std::string err_msg =  "No PackerPaletteCreator with the name '" + packer_palette_name + "' has been registered with the PackerPaletteFactory";
		throw CREATE_EXCEPTION(utility::excn::Exception,  err_msg );
	}
	auto iter = packer_palette_creator_map_.find( packer_palette_name );
	iter->second->provide_xml_schema( xsd );
}

void
PackerPaletteFactory::newPackerPalettes( PackerPaletteOPs & ppops, basic::datacache::DataMap & datamap, std::string const & tagfilename ) const
{
	utility::io::izstream fin;
	fin.open( tagfilename );
	runtime_assert( fin.good() );
	TagCOP tag = utility::tag::Tag::create(fin);
	fin.close();
	TR << "PackerPaletteFactory parsing " << tagfilename << " to create PackerPalettes:" << std::endl;
	TR << tag << std::endl;
	newPackerPalettes( ppops, datamap, tag );
}

/// @brief Create a packer palette based on global defaults, and return an owning pointer to it.
/// @details By default, makes a DefaultPackerPalette. (If the user provides -packer_palette:NCAA_expanded,
/// or equivalent local options, makes a NCAADefaultPackerPalette.)
/// If the user provides options for additional residues, makes a CustomBaseTypePackerPalette.
PackerPaletteOP
PackerPaletteFactory::create_packer_palette_from_global_defaults() const {
	boost::function< PackerPaletteOP () > creator( boost::bind( &PackerPaletteFactory::initialize_packer_palette_from_global_defaults ) );
	utility::thread::safely_create_load_once_object_by_OP( creator, default_palette_ , SAFELY_PASS_MUTEX(default_palette_mutex_), SAFELY_PASS_THREADSAFETY_BOOL(default_palette_bool_) ); //Creates this once in a threadsafe manner, iff it hasn't been created.  Otherwise, returns already-created object.
	return default_palette_->clone();
}

/// @brief Create the PackerPalette that we clone whenever create_packer_palette_from_global_defaults() is called.
/// @details This function is intended ONLY to be called in a threadsafe manner from create_packer_palette_from_global_defaults().
/// @note This function reads from the global options system, and can read from disk.
PackerPaletteOP
PackerPaletteFactory::initialize_packer_palette_from_global_defaults() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	static const std::string errmsg( "Error in core::pack::palette::PackerPaletteFactory::initialize_packer_palette_from_global_defaults(): " );

	// How will this interact with a NCAADefaultPackerPalette? I'd hope its universe would
	// be strictly bigger, right?
	runtime_assert_string_msg( !( option[ packing::packer_palette::extra_base_type_file ].user()
		&& option[ packing::packer_palette::NCAA_expanded ].user() ),
		errmsg + "The -extra_base_type_file and -NCAA_expanded flags are mutually incompatible for now." );

	if ( option[ packing::packer_palette::extra_base_type_file ].user() ) {
		std::string const filename( option[ packing::packer_palette::extra_base_type_file ].value() );
		runtime_assert_string_msg( !filename.empty(), errmsg + "The filename provided with the \"-packing:packer_palette:extra_base_type_file\" option cannot be empty!" );

		std::string file_contents;
		try {
			file_contents = utility::file_contents( filename );
		} catch( utility::excn::Exception & exception ) {
			utility_exit_with_message( errmsg + "Could not read file \"" + filename + "\".  Error:\n" + exception.msg() );
		}

		CustomBaseTypePackerPaletteOP new_palette( utility::pointer::make_shared< CustomBaseTypePackerPalette >() );
		try {
			new_palette->initialize_from_file_contents( file_contents, filename );
		} catch( utility::excn::Exception & exception ) {
			utility_exit_with_message( errmsg + "Could not parse file \"" + filename + "\".  Error:\n" + exception.msg() );
		}
		return new_palette;
	}

	if ( option[ packing::packer_palette::NCAA_expanded ].user() ) {
		return utility::pointer::make_shared< NCAADefaultPackerPalette >();
	}

	return utility::pointer::make_shared< DefaultPackerPalette >();
}

} //namespace palette
} //namespace pack
} //namespace core
