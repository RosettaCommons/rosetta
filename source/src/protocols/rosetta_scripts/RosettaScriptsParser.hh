// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
	typedef protocols::moves::MoverFactory MoverFactory;
	typedef protocols::moves::MoverFactoryOP MoverFactoryOP;
	typedef std::pair<std::string, std::string> ImportTagName;

public:
	RosettaScriptsParser();
	virtual ~RosettaScriptsParser();

	virtual
	bool
	generate_mover_from_pose( protocols::jd2::JobCOP job, core::pose::Pose & pose, protocols::moves::MoverOP & mover, bool new_input, std::string const xml_fname );

	MoverOP generate_mover_for_protocol(Pose & pose, bool & modified_pose, utility::tag::TagCOP protocol_tag);

	//@brief Temporary hook into parsing machinery with pose reference
	MoverOP parse_protocol_tag(Pose & pose, utility::tag::TagCOP protocol_tag);

	//@brief Temporary hook into parsing machinery w/o pose reference.
	MoverOP parse_protocol_tag(utility::tag::TagCOP protocol_tag);

	void register_factory_prototypes();
	
	void instantiate_filter  (utility::tag::TagCOP const & tag_ptr, basic::datacache::DataMap & data, protocols::filters::Filters_map & filters, protocols::moves::Movers_map & movers, core::pose::Pose & pose);
	void instantiate_mover   (utility::tag::TagCOP const & tag_ptr, basic::datacache::DataMap & data, protocols::filters::Filters_map & filters, protocols::moves::Movers_map & movers, core::pose::Pose & pose);

	utility::tag::TagCOP find_rosettascript_tag(utility::tag::TagCOP rootTag, const std::string & section_name, const std::string & option_name, const std::string & option_value);
	
	void import_tags(std::set< ImportTagName > & import_tag_names, utility::tag::TagCOP & my_tag, basic::datacache::DataMap & data, protocols::filters::Filters_map & filters, protocols::moves::Movers_map & movers, core::pose::Pose & pose);
	
private:

	static
	void
	substitute_variables_in_stream( std::istream & in, utility::options::StringVectorOption const& script_vars, std::stringstream & out);

}; // Parser

} // namespace rosetta_scripts
} // namespace protocols

#endif //INCLUDED_protocols_jd2_RosettaScriptsParser_HH
