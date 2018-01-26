// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/LoopsFileLoader.cc
/// @brief
/// @author

// unit headers
#include <protocols/loops/LoopsFileLoader.hh>
#include <protocols/loops/LoopsFileLoaderCreator.hh>

// package headers
#include <protocols/loops/LoopsFileIO.hh>

// Basic headers
#include <basic/resource_manager/loader_schemas.hh>

// utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

//C++ headers
#include <istream>

namespace protocols {
namespace loops {

LoopsFileLoader::LoopsFileLoader() = default;
LoopsFileLoader::~LoopsFileLoader() = default;

basic::resource_manager::ResourceCOP
LoopsFileLoader::create_resource(
	basic::resource_manager::ResourceManager &,
	utility::tag::TagCOP resource_tag,
	std::string const & input_id,
	std::istream & loopfstream
) const
{
	bool prohibit_single_residue_loops = resource_tag->getOption< bool >( "prohibit_single_residue_loops", true );
	LoopsFileIO lfio;
	LoopsFileDataOP lfd = lfio.read_loop_file_stream( loopfstream, input_id, prohibit_single_residue_loops );

	return lfd;
}

std::string
LoopsFileLoader::classname()
{
	return "LoopsFile";
}

void
LoopsFileLoader::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "prohibit_single_residue_loops", xsct_rosetta_bool,
		"Prevent loop definitions that describe only a single residue? Defaults to true", "true" );

	basic::resource_manager::resource_loader_xsd_type_definition_w_attributes( xsd,
		classname(), "If a different LoopsFile should be loaded for different jobs, then using the"
		" ResourceManager to specify which loops file goes with with job can be done; this may"
		" be strictly unnecessary because specifying the filename directly with an XML tag and"
		" the per-job script_vars could substitute in the filename you want for a given job.",
		attlist );
}


/// @details Return an owning pointer to a newly constructed default instance of LoopsFileLoader.
basic::resource_manager::ResourceLoaderOP LoopsFileLoaderCreator::create_resource_loader() const
{
	return basic::resource_manager::ResourceLoaderOP( new LoopsFileLoader() );
}

/// @details Return a string specifying the type of %ResourceLoader to create (LoopsFile).
std::string LoopsFileLoaderCreator::loader_type() const
{
	return LoopsFileLoader::classname();
}

void
LoopsFileLoaderCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	LoopsFileLoader::provide_xml_schema( xsd );
}

} // namespace loops
} // namespace protocols
