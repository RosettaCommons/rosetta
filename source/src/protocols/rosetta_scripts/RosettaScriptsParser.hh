// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/rosetta_scripts/RosettaScriptsParser.hh
/// @brief  header file for Parser class, part of August 2008 job distributor
/// @author Sarel Fleishman sarelf@u.washington.edu
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com). JD3/PyRosetta compatability, removal of JD2 code, general updates and docs.
///   Simplification - generate_mover_and_apply_to_pose functions.

#ifndef INCLUDED_protocols_rosetta_scripts_RosettaScriptsParser_hh
#define INCLUDED_protocols_rosetta_scripts_RosettaScriptsParser_hh

// Unit Headers
#include <protocols/rosetta_scripts/RosettaScriptsParser.fwd.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.fwd.hh>
#include <protocols/rosetta_scripts/XmlObjects.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Basic headers
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/resource_manager/ResourceManager.fwd.hh>

// Utility headers
#include <utility/options/StringVectorOption.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/tag/XMLSchemaValidation.fwd.hh>
#include <utility/VirtualBase.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/keys/OptionKeyList.fwd.hh>

// C++ headers
#include <iosfwd>
#include <set>

namespace protocols {
namespace rosetta_scripts {


/// @brief Reading the xml file and generating the mover
/// @details Uses the Tag interface to the xml reader library in boost to parse
/// an xml file that contains design protocol information.
///
/// SCOREFXNS provides a way to define scorefunctions as they are defined in the
/// rosetta database, using the weights/patch convenctions. Several default
/// scorefunctions are preset and can be used without defining them explicitly.
///
/// FILTERS defines a set of filters that can be used together with the
/// dockdesign movers to prune out poses that don't meet certain criteria
///
/// MOVERS defines the movers that will be used
///
/// PROTOCOLS is the only order-sensitive section where subsequent movers and
/// filters are expected to be defined. These movers and filters were defined
/// in the previous two sections. The order in which the protocols are specified
/// by the user will be maintained by the DockDesign mover.
///
/// APPLY_TO_POSE -- DEPRECATED/Removed. This section used to allow for applying certain
/// movers to the pose prior to protocol start. This functionality has been removed.
/// Anything other than an empty APPLY_TO_POSE section will raise an error.
///
/// xi:include statements can be placed anywhere to trigger read-in of another
/// RosettaScripts XML file.  Contents of the file replace the xi:include
/// statement.
class RosettaScriptsParser : public utility::VirtualBase
{
public:
	typedef utility::vector0< utility::tag::TagCOP > TagCOPs;
	typedef utility::tag::TagCOP TagCOP;
	typedef moves::MoverOP MoverOP;
	typedef std::pair<std::string, std::string> ImportTagName;

public:
	RosettaScriptsParser();

	~RosettaScriptsParser() override;

	/// @brief Convenience function to construct a ParsedProtocol from a file
	/// @details Reads in an RosettaScript from file, replaces script_vars if present,
	///   Then returns the full ParsedProtocol
	///   See also protocols::rosetta_scripts::XmlObjects
	ParsedProtocolOP
	generate_mover(
		std::string const & xml_fname,
		std::string const & current_input_name = "",
		std::string const & current_output_name = ""
	);

	/// @brief Convenience function to construct a ParsedProtocol from a file
	/// @details Reads in an RosettaScript from file, replaces script_vars if present,
	///   Then returns the full ParsedProtocol
	///   See also protocols::rosetta_scripts::XmlObjects
	ParsedProtocolOP
	generate_mover(
		std::string const & xml_fname,
		utility::vector1< std::string > const & script_vars,
		std::string const & current_input_name = "",
		std::string const & current_output_name = ""
	);

	/// @brief Convenience function to construct a ParsedProtocol from a file
	/// @details Reads in an RosettaScript from file, replaces script_vars if present,
	///   Then returns the full ParsedProtocol
	///   See also protocols::rosetta_scripts::XmlObjects.
	///   Any Resources requested in the RESOURCES block will be
	///   retrieved from the ResourceManager and placed in the DataMap.
	ParsedProtocolOP
	generate_mover(
		utility::options::OptionCollection const & options,
		std::string const & current_input_name,
		std::string const & current_output_name,
		basic::resource_manager::ResourceManagerOP resource_manager
	);

	/// @brief Convenience function to construct a ParsedProtocol from a file
	/// @details Reads in an RosettaScript from file, replaces script_vars if present,
	///   Then returns the full ParsedProtocol
	///   See also protocols::rosetta_scripts::XmlObjects
	ParsedProtocolOP
	generate_mover(
		utility::options::OptionCollection const & options,
		std::string const & xml_fname,
		std::string const & current_input_name = "",
		std::string const & current_output_name = "",
		XmlObjectsOP xml_objects = nullptr,
		basic::resource_manager::ResourceManagerOP resource_manager = nullptr
	);

	/// @brief Convenience function to construct a ParsedProtocol from a string
	/// @details Reads in an RosettaScript from file, replaces script_vars if present,
	///   Then returns the full ParsedProtocol
	///   See also protocols::rosetta_scripts::XmlObjects
	ParsedProtocolOP
	generate_mover_xml_string(
		utility::options::OptionCollection const & options,
		std::string const & xml_text,
		std::string const & current_input_name,
		std::string const & current_output_name,
		XmlObjectsOP xml_objects /* = nullptr */ );


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/// @brief Convenience function to construct a ParsedProtocol from a file
	/// @details Reads in an RosettaScript from file, replaces script_vars if present,
	///   Then returns the full ParsedProtocol
	///   See also protocols::rosetta_scripts::XmlObjects.
	///   Any Resources requested in the RESOURCES block will be
	///   retrieved from the ResourceManager and placed in the DataMap.
	///
	/// Places an OptionCollectionCOP into the datamap for local option parsing.
	///
	ParsedProtocolOP
	generate_mover(
		utility::options::OptionCollectionCOP options,
		std::string const & current_input_name,
		std::string const & current_output_name,
		basic::resource_manager::ResourceManagerOP resource_manager
	);

	/// @brief Convenience function to construct a ParsedProtocol from a file
	/// @details Reads in an RosettaScript from file, replaces script_vars if present,
	///   Then returns the full ParsedProtocol
	///   See also protocols::rosetta_scripts::XmlObjects
	///
	/// Places an OptionCollectionCOP into the datamap for local option parsing.
	///
	ParsedProtocolOP
	generate_mover(
		utility::options::OptionCollectionCOP options,
		std::string const & xml_fname,
		std::string const & current_input_name = "",
		std::string const & current_output_name = "",
		XmlObjectsOP xml_objects = nullptr,
		basic::resource_manager::ResourceManagerOP resource_manager = nullptr
	);

	/// @brief Convenience function to construct a ParsedProtocol from a string
	/// @details Reads in an RosettaScript from file, replaces script_vars if present,
	///   Then returns the full ParsedProtocol
	///   See also protocols::rosetta_scripts::XmlObjects
	///
	/// Places an OptionCollectionCOP into the datamap for local option parsing.
	///
	ParsedProtocolOP
	generate_mover_xml_string(
		utility::options::OptionCollectionCOP options,
		std::string const & xml_text,
		std::string const & current_input_name,
		std::string const & current_output_name,
		XmlObjectsOP xml_objects /* = nullptr */ );

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	///@brief
	///  Create a tag from an XML file.  Read from Script Vars options from a local options collection.
	TagCOP
	create_tag_from_xml( std::string const & xml_fname, utility::options::OptionCollection const & options);

	///@brief
	///  Create a tag from an XML file.  Use Script Vars bariable to do any substitutions.
	///  Script Var String: "xx=yy" for %xx% replacement in XML script.
	TagCOP
	create_tag_from_xml( std::string const & xml_fname, utility::vector1< std::string > const & script_vars );

	///@brief
	///  Create a tag from an XML string.  Read from Script Vars options from a local options collection.
	TagCOP
	create_tag_from_xml_string( std::string const & xml_text, utility::options::OptionCollection const & options);

	///@brief
	///  Create a tag from an XML string.  Use Script Vars bariable to do any substitutions.
	///  Script Var String: "xx=yy" for %xx% replacement in XML script.
	TagCOP
	create_tag_from_xml_string( std::string const & xml_text, utility::vector1< std::string > const & script_vars );

	///@brief
	///
	/// Create the ParsedProtocolMover for the protocol using a tag.
	///
	///@details
	///
	/// Creates mover via standard parsing machinery
	ParsedProtocolOP parse_protocol_tag( utility::tag::TagCOP protocol_tag, utility::options::OptionCollection const & options );

	///@brief
	///
	/// Create the ParsedProtocolMover for the protocol using a tag.
	///
	///@details
	///
	/// Load the datamap with the predefined object,
	/// then go through tags, generating object and the ParsedProtocol Mover.
	/// READS from LOCAL OptionCollection.
	ParsedProtocolOP
	generate_mover_for_protocol(
		utility::tag::TagCOP protocol_tag,
		utility::options::OptionCollection const & options,
		std::string const & current_input_name = "",
		std::string const & current_output_name = "",
		XmlObjectsOP xml_objects = nullptr,
		basic::resource_manager::ResourceManagerOP resource_manager = nullptr
	);

	///@brief
	///
	/// Create the ParsedProtocolMover for the protocol using a tag.
	///
	///@details
	///
	/// Uses the passed DataMap for passing to parse_my_tag.
	/// This function does the heavy lifting.
	///
	/// Does not add any predefined objects to the data map (though in-tag defined objects will be added.)
	ParsedProtocolOP
	generate_mover_for_protocol(
		utility::tag::TagCOP protocol_tag,
		basic::datacache::DataMap & data,
		utility::options::OptionCollection const & options,
		XmlObjectsOP xml_objects = nullptr,
		basic::resource_manager::ResourceManagerOP resource_manager = nullptr
	);

	///@brief
	///
	/// Create the ParsedProtocolMover for the protocol using a tag.
	///
	///@details
	///
	/// This interface is preserved for (PyRosetta) backwards compatibility only.
	/// The Pose and bool (modified_pose) parameters are now unused.
	ParsedProtocolOP
	generate_mover_for_protocol(
		core::pose::Pose &,
		bool &,
		utility::tag::TagCOP protocol_tag,
		utility::options::OptionCollection const & options,
		std::string const & current_input_name = "",
		std::string const & current_output_name = "",
		XmlObjectsOP xml_objects = nullptr,
		basic::resource_manager::ResourceManagerOP resource_manager = nullptr
	) { return generate_mover_for_protocol(protocol_tag, options, current_input_name, current_output_name, xml_objects, resource_manager ); }

	/// @brief Initialize the passed datamap with the standard set of pre-defined RosettaScripts objects
	void
	initialize_data_map(
		basic::datacache::DataMap & data,
		utility::options::OptionCollection const & options
	);


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/// @brief Actually read in the XML file.  Called recursively to read in XML files that
	/// this XML file includes.  At the end of this operation, substituted_contents contains the contents
	/// of the XML file, with all xi:includes replaced with the contents of included XML
	/// files.  Files are opened and closed here.
	/// @details Note that filenames_encountered is passed by copy rather than by reference
	/// DELIBERATELY.  This is so that it remains a list of parent files, so that only
	/// circular dependencies (attempts to include one's own parent, grandparent, etc.) are
	/// detected.
	/// If xml_text_for_top_level is set the filename passed will not be read and instead
	/// xml_text_for_top_level will be used as though it was the contents of the first file.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	/// @author Rocco Moretti (rmorettiase@gmail.com)
	void
	read_in_and_recursively_replace_includes(
		std::string const &filename,
		std::string & substituted_contents, // Return-by-reference
		utility::vector1 < std::string > filenames_encountered,
		core::Size const recursion_level,
		bool const do_not_recurse = false,
		std::string const & xml_text_for_top_level = ""
	) const;

	static void register_factory_prototypes();

	void
	parse_resources_tag(
		utility::tag::TagCOP resources_tag,
		basic::datacache::DataMap & data,
		basic::resource_manager::ResourceManagerOP resource_manager
	) const;

	void
	parse_apply_to_pose_tag(
		utility::tag::TagCOP apply_tag,
		basic::datacache::DataMap & data
	) const;

	void instantiate_filter(
		utility::tag::TagCOP const & tag_ptr,
		basic::datacache::DataMap & data
	) const;

	void instantiate_mover(
		utility::tag::TagCOP const & tag_ptr,
		basic::datacache::DataMap & data
	) const;

	/// @brief Instantiate a PackerPalette object.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void instantiate_packer_palette(
		utility::tag::TagCOP const & tag_ptr,
		basic::datacache::DataMap & data
	) const;

	void instantiate_taskoperation(
		utility::tag::TagCOP const & tag_ptr,
		basic::datacache::DataMap & data
	) const;

	utility::tag::TagCOP
	find_rosettascript_tag(
		utility::tag::TagCOP rootTag,
		const std::string & section_name,
		const std::string & option_name,
		const std::string & option_value
	);

	void
	parse_import_tag(
		utility::tag::TagCOP import_tag,
		basic::datacache::DataMap & data
	) const;

	void import_tags(
		std::set< ImportTagName > & import_tag_names,
		utility::tag::TagCOP my_tag,
		basic::datacache::DataMap & data
	) const;

	void
	validate_input_script_against_xsd(
		std::string const & fname,
		std::stringstream const & input_xml
	) const;

	static std::string xsd_for_rosetta_scripts();

	static void write_ROSETTASCRIPTS_complex_type( utility::tag::XMLSchemaDefinition & );
	static std::string rosetta_scripts_element_name(); // returns "ROSETTASCRIPTS"
	static std::string rosetta_scripts_complex_type_naming_func( std::string const & element_name );

	/// @brief List the command-line options used by this class.
	static
	void
	list_options_read( utility::options::OptionKeyList & read_options );

	void
	set_recursion_limit( utility::options::OptionCollection const & options );

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	///@brief
	/// Original, JD2-compatible mover generating function.
	///
	///@detals
	/// NOT accessible in PyRosetta.
	/// For general use and PyRosetta, use mover-returning generate_mover() functions.
	/// READS FROM GLOBAL OPTIONS.
	void
	generate_mover(
		MoverOP & mover,
		bool new_input,
		std::string const & xml_file,
		std::string const & current_input_name = "",
		std::string const & current_output_name = "",
		bool guarantee_new_mover = false
	);

	///@brief Checks the LOCAL options collection to see if a protocol has been set.
	bool
	protocol_option_set(utility::options::OptionCollection const & options) const;

	///@brief Substitute xx==yy using the XML syntax %%xx%%.
	//// script_vars is a list of strings of "xx==yy".
	static
	void
	substitute_variables_in_stream(
		std::istream & in,
		utility::vector1< std::string > const & script_vars,
		std::stringstream & out
	);

private: //Methods

	/// @brief Is the current tag one with a slash at the end?
	/// @details Starting from endpos, crawl backward until we hit a
	/// '/' or a non-whitespace character.  Return true if and only
	/// if it's a slash, false otherwise.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	bool
	has_slash_at_end(
		std::string const &str,
		core::Size const endpos
	) const;

	/// @brief Is the current tag one with something that will be replaced by variable replacement?
	/// @details Starting from endpos, crawl backward until we hit a '<' or '%'.  Return 'false' if
	/// the former is preceded by a percent sign and 'true' if we hit the latter.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	bool
	has_double_percent_signs(
		std::string const &str,
		core::Size const endpos
	) const;


	/// @brief To be called after the full xml tree has been loaded into a string from create_tag_from_xml*
	/// @details xml_fname only used for tracer output.
	TagCOP
	finish_creating_tag(
		std::string const & xml_fname,
		utility::vector1< std::string > const & script_vars,
		std::string const & substituted_contents ) const;

	/// @brief Load a native into this class, which will be placed into the datamap as a resource upon generating mover(s).
	void
	load_native( std::string const & fname );

private: //Varaibles

	/// @brief The object used to validate XML input files against the schema
	utility::tag::XMLSchemaValidatorOP validator_;

	/// @brief The depth of the rabbit hole one can find oneself in by including XML files that include XML files.
	/// Defaults to 8.
	core::Size recursion_limit_;

	core::pose::PoseOP native_pose_ = nullptr; //Loaded from options and added to datamap if present
	utility::options::OptionCollectionCOP local_options_ = nullptr; //Passed to datamap for local options access for RosettaScript elements.

}; // Parser

} // namespace rosetta_scripts
} // namespace protocols

#endif // INCLUDE GUARD
