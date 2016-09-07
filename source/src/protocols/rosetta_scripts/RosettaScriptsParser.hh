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


#ifndef INCLUDED_protocols_rosetta_scripts_RosettaScriptsParser_hh
#define INCLUDED_protocols_rosetta_scripts_RosettaScriptsParser_hh

// Unit Headers
#include <protocols/filters/Filter.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <protocols/jd2/Parser.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MoverFactory.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/options/StringVectorOption.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>

// C++ headers
#include <iostream>
#include <set>

namespace protocols {
namespace rosetta_scripts {

/// @brief Reading the xml file and generating the mover
class RosettaScriptsParser : public protocols::jd2::Parser
{
public:
	typedef utility::vector0< utility::tag::TagCOP > TagCOPs;
	typedef protocols::moves::MoverFactory MoverFactory;
	typedef protocols::moves::MoverFactoryOP MoverFactoryOP;
	typedef std::pair<std::string, std::string> ImportTagName;

public:
	RosettaScriptsParser();
	~RosettaScriptsParser() override;

	/// @brief Actually read in the XML file.  Called recursively to read in XML files that
	/// this XML file includes.  At the end of this operation, fin contains the contents
	/// of the XML file, with all xi:includes replaced with the contents of included XML
	/// files.  Files are opened and closed here.
	/// @details Note that filenames_encountered is passed by copy rather than by reference
	/// DELIBERATELY.  This is so that it remains a list of parent files, so that only
	/// circular dependencies (attempts to include one's own parent, grandparent, etc.) are
	/// detected.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void
	read_in_and_recursively_replace_includes(
		std::string const &filename,
		std::stringstream &fin,
		utility::vector1 < std::string > filenames_encountered
	) const;

	/// @brief Open the file given by xml_fname and construct a mover representing
	/// the script contained in that file. If new_input is true, run the APPLY_TO_POSE
	/// block on the input mover.  If both new_input and guarantee_new_mover are false,
	/// then the input mover is considered up-to-date and the file is not re-read.
	
	bool
	generate_mover_from_pose(
		protocols::jd2::JobCOP job,
		core::pose::Pose & pose,
		protocols::moves::MoverOP & mover,
		bool new_input,
		std::string const & xml_fname,
		bool guarantee_new_mover = false
	) override;

	MoverOP
	generate_mover_for_protocol(
		Pose & pose,
		bool & modified_pose,
		utility::tag::TagCOP protocol_tag
	);

	//@brief Temporary hook into parsing machinery with pose reference
	MoverOP parse_protocol_tag( Pose & pose, utility::tag::TagCOP protocol_tag );

	//@brief Temporary hook into parsing machinery w/o pose reference.
	MoverOP parse_protocol_tag( utility::tag::TagCOP protocol_tag );

	void register_factory_prototypes();

	void instantiate_filter(
		utility::tag::TagCOP const & tag_ptr,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map & filters,
		protocols::moves::Movers_map & movers,
		core::pose::Pose & pose
	);

	void instantiate_mover(
		utility::tag::TagCOP const & tag_ptr,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map & filters,
		protocols::moves::Movers_map & movers,
		core::pose::Pose & pose
	);

	void instantiate_taskoperation(
		utility::tag::TagCOP const & tag_ptr,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map & filters,
		protocols::moves::Movers_map & movers,
		core::pose::Pose & pose
	);

	utility::tag::TagCOP
	find_rosettascript_tag(
		utility::tag::TagCOP rootTag,
		const std::string & section_name,
		const std::string & option_name,
		const std::string & option_value
	);

	void import_tags(
		std::set< ImportTagName > & import_tag_names,
		utility::tag::TagCOP & my_tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map & filters,
		protocols::moves::Movers_map & movers,
		core::pose::Pose & pose
	);

private:

	static
	void
	substitute_variables_in_stream(
		std::istream & in,
		utility::options::StringVectorOption const& script_vars,
		std::stringstream & out
	);

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

}; // Parser

} // namespace rosetta_scripts
} // namespace protocols

#endif //INCLUDED_protocols_jd2_RosettaScriptsParser_HH
