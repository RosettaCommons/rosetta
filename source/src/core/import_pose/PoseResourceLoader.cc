// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/import_pose/PoseResourceLoader.cc
/// @brief
/// @author Brian D. Weitzner brian.weitzner@gmail.com

//unit headers
#include <core/import_pose/PoseResourceLoader.hh>
#include <core/import_pose/PoseResourceLoaderCreator.hh>

//package headers
#include <core/import_pose/PoseResource.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/import_pose_options.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentStruct.hh>

#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/ResourceLoaderFactory.hh>
//#include <core/chemical/ResidueType.hh>
//#include <core/chemical/ResidueTypeSet.hh>
//#include <core/chemical/PoseResidueTypeSet.hh>
//#include <core/conformation/Conformation.hh>

// Basic headers
#include <basic/Tracer.hh>

//utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/pointer/memory.hh>

// numeric headers

static basic::Tracer TR( "core.import_pose.PoseResourceLoader" );

namespace core {
namespace import_pose {

PoseResourceLoader::PoseResourceLoader() {}
PoseResourceLoader::~PoseResourceLoader() = default;

basic::resource_manager::ResourceCOP
PoseResourceLoader::create_resource(
	basic::resource_manager::ResourceManager &,
	utility::tag::TagCOP resource_tag,
	std::string const & input_id,
	std::istream & input_stream
) const
{
	pose::PoseOP pose( new pose::Pose );
	using namespace utility::tag;
	runtime_assert( resource_tag->getTags().size() == 1 );
	TagCOP child( resource_tag->getTags()[ 0 ] );
	if ( child->getName() == "PDB" ) {
		ImportPoseOptions opts;
		if ( child->getTags().size() > 0 ) {
			TagCOP grandchild( child->getTags()[0] );
			opts.parse_my_tag( grandchild );
		}
		pose_from_pdb_stream( *pose, input_stream, input_id, opts );
	} else if ( child->getName() == "Silent" ) {
		using namespace core::io::silent;
		SilentFileOptions opts;
		runtime_assert( child->hasOption( "tag" ) );
		std::string silent_tag = child->getOption< std::string >( "tag" );
		utility::vector1< std::string > silent_tags( 1, silent_tag );
		if ( child->getTags().size() > 0 ) {
			TagCOP grandchild( child->getTags()[0] );
			opts.read_from_tag( grandchild );
		}
		SilentFileData sfd( opts );
		sfd.read_stream( input_stream, silent_tags, false, input_id );
		SilentStructOP silent_struct = sfd[ silent_tag ];
		silent_struct->fill_pose( *pose );
	}
	return utility::pointer::make_shared< PoseResource >( pose );
}

std::string
PoseResourceLoader::classname()
{
	return "Pose";
}

std::string
ctname_for_pose_resource_loader( std::string const & element_name )
{
	return "pose_resource_" + element_name + "_type";
}

void
PoseResourceLoader::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	// OK: two possible subtags: PDB and SilentFile.
	// PDB has no attributes and has one (optional) subtag PDBOptions
	// SilentFile has one attribute, tag, and one (optional) subtag SilentFileOptions

	using namespace utility::tag;

	AttributeList pdb_opts_attributes;
	ImportPoseOptions::append_schema_attributes( pdb_opts_attributes );

	XMLSchemaSimpleSubelementList pdb_subelements;
	pdb_subelements.add_simple_subelement( "PDBOptions", pdb_opts_attributes, "Configuration stating how to create a Pose from the indicated PDB/MMCIF file" );
	XMLSchemaComplexTypeGenerator pdb_ct;
	pdb_ct.element_name( "PDB" )
		.complex_type_naming_func( & ctname_for_pose_resource_loader )
		.description( "Initialize a Pose from a PDB- or MMCIF-formatted file" )
		.set_subelements_single_appearance_optional( pdb_subelements )
		.write_complex_type_to_schema( xsd );

	AttributeList silent_file_opts_attributes;
	core::io::silent::SilentFileOptions::append_attributes_for_tag_parsing( xsd, silent_file_opts_attributes );
	XMLSchemaSimpleSubelementList silent_file_subelements;
	silent_file_subelements.add_simple_subelement( "SilentFileOptions", silent_file_opts_attributes, "Configuration stating how to create a Pose from the indicated silent file" );
	AttributeList silent_file_attributes;
	silent_file_attributes
		+ XMLSchemaAttribute::required_attribute( "tag", xs_string, "The tag for the structure to be read from the silent file; the file name is indicated by the input_id, but the tag specifies which"
		" structure from that file should be read" );
	XMLSchemaComplexTypeGenerator silent_file_ct;
	silent_file_ct.element_name( "Silent" )
		.complex_type_naming_func( & ctname_for_pose_resource_loader )
		.description( "Initialize a Pose from a (rosetta-specific) silent-file-formatted file/input stream" )
		.set_subelements_single_appearance_optional( silent_file_subelements )
		.write_complex_type_to_schema( xsd );

	XMLSchemaSimpleSubelementList pose_resource_subelements;
	pose_resource_subelements
		.add_already_defined_subelement( "PDB", & ctname_for_pose_resource_loader )
		.add_already_defined_subelement( "Silent", & ctname_for_pose_resource_loader );
	XMLSchemaComplexTypeGenerator pose_resource_ct;
	pose_resource_ct.element_name( "Pose" )
		.complex_type_naming_func( & basic::resource_manager::ResourceLoaderFactory::complex_type_name_for_loader )
		.description( "Create a Pose as a resource so that it only needs to be loaded in from disk a single time; note that a Pose has mutable"
		" data members, and so the PoseResourceLoader creates a PoseResourceCOP instead of a PoseCOP. When a PoseResourceCOP is requested from"
		" the ResourceManager, the PoseResource only grants access to the Pose it contains by creating a *copy* of that Pose. This means that"
		" a Pose held by the ResourceManager will have multiple instances in memory. Truly bitwise-const classes (i.e. classes with no mutable"
		" data) can be guaranteed to have only a single instance in memory." )
		.set_subelements_pick_one( pose_resource_subelements )
		.write_complex_type_to_schema( xsd );
}

basic::resource_manager::ResourceLoaderOP PoseResourceLoaderCreator::create_resource_loader() const
{
	return basic::resource_manager::ResourceLoaderOP( new PoseResourceLoader() );
}

std::string PoseResourceLoaderCreator::loader_type() const
{
	return PoseResourceLoader::classname();
}

void
PoseResourceLoaderCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PoseResourceLoader::provide_xml_schema( xsd );
}

} // namespace import_pose
} // namespace core
