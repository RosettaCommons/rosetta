// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loophash/LoopHashLibraryLoader.cc
/// @brief Load the Loop Hash library using the resource manager
/// @author Tim Jacobs

//unit headers
#include <protocols/loophash/LoopHashLibraryLoader.hh>
#include <protocols/loophash/LoopHashLibraryLoaderCreator.hh>
#include <protocols/loophash/LoopHashLibrary.hh>

//package headers
#include <basic/resource_manager/loader_schemas.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

//C++ headers
#include <istream>

namespace protocols {
namespace loophash {


LoopHashLibraryLoader::LoopHashLibraryLoader() {}

LoopHashLibraryLoader::~LoopHashLibraryLoader() = default; // why are we doing this again?

/// @details Ensure the %ResourceOptions is a LoopHashLibraryOptions instance and construct a new LoopHashLibrary from
/// it.  The locator_id and istream are not used.
/// @throws utility::excn::EXCN_Msg_Exception
basic::resource_manager::ResourceCOP
LoopHashLibraryLoader::create_resource(
	basic::resource_manager::ResourceManager &,
	utility::tag::TagCOP resource_tag,
	std::string const &,
	std::istream &
) const {

	utility::vector1< core::Size > loop_sizes;
	if ( resource_tag->hasOption("loop_sizes") ) {
		utility::vector1< std::string > loop_sizes_strings =
			utility::string_split( resource_tag->getOption< std::string >("loop_sizes"), ',');
		for ( core::Size i=1; i<=loop_sizes_strings.size(); ++i ) {
			// Error handling needed? No: rely on the XSD validation to ensure that a positive integer
			// has been given for this attribute.
			loop_sizes.push_back( utility::string2int( loop_sizes_strings[i] ) );
		}
	} else {
		throw CREATE_EXCEPTION( utility::excn::Exception, "You must provide a 'loop_sizes' option to the LoopHashLibrary tag");
	}

	LoopHashLibraryOP lh_library( new LoopHashLibrary( loop_sizes ) );
	lh_library->load_mergeddb(); // why is the input file hard coded?

	// From here forward, the library may not be changed.
	return lh_library;
}

std::string
LoopHashLibraryLoader::classname()
{
	return "LoopHashLibrary";
}

void
LoopHashLibraryLoader::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute( "loop_sizes", xsct_nnegative_int_cslist, "The"
		" list of loop sizes that will be used to extract loops out of the loop-hash structure database" );
	basic::resource_manager::resource_loader_xsd_type_definition_w_attributes( xsd,
		classname(), "The LoopHashLibrary can be used to perform minimal-motion 'loop' conformational"
		" sampling. Indeed, this sampling need not be limited to only loops, but can help sample"
		" conformations of helices and strands. The library itself is quite large, and therefore"
		" should be managed with the ResourceManager", attlist );
}

/// @details Return an owning pointer to a newly constructed default instance of LoopHashLibraryLoader.
basic::resource_manager::ResourceLoaderOP
LoopHashLibraryLoaderCreator::create_resource_loader() const
{
	return basic::resource_manager::ResourceLoaderOP( new LoopHashLibraryLoader() );
}

/// @details Return a string specifying the type of %ResourceLoader to create (LoopHashLibrary).
std::string LoopHashLibraryLoaderCreator::loader_type() const
{
	return LoopHashLibraryLoader::classname();
}

void
LoopHashLibraryLoaderCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	LoopHashLibraryLoader::provide_xml_schema( xsd );
}


} // namespace
} // namespace
