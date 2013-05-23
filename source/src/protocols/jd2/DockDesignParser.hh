// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/DockDesignParser.hh
/// @brief  header file for Parser class, part of August 2008 job distributor
/// @author Sarel Fleishman sarelf@u.washington.edu


#ifndef INCLUDED_protocols_jd2_DockDesignParser_hh
#define INCLUDED_protocols_jd2_DockDesignParser_hh

//unit headers
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <protocols/jd2/Parser.hh>

#include <protocols/moves/MoverFactory.fwd.hh>
//#include <protocols/protein_interface_design/DockDesignFilterFactory.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>

//utility headers
// AUTO-REMOVED #include <basic/options/option.hh>

#include <utility/vector1.hh>
#include <utility/options/StringVectorOption.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <iostream>


namespace protocols {
namespace jd2 {

/// @brief Reading the xml file and generating the mover
class DockDesignParser : public Parser
{
public:
	typedef protocols::moves::MoverFactory MoverFactory;
	typedef protocols::moves::MoverFactoryOP MoverFactoryOP;
	//typedef protocols::protein_interface_design::DockDesignFilterFactory DockDesignFilterFactory;
	//typedef protocols::protein_interface_design::DockDesignFilterFactoryOP DockDesignFilterFactoryOP;

public:
	DockDesignParser();
	virtual ~DockDesignParser();

	virtual
	bool
	generate_mover_from_pose( protocols::jd2::JobCOP job, core::pose::Pose & pose, protocols::moves::MoverOP & mover, bool new_input, std::string const xml_fname );

	MoverOP generate_mover_for_protocol(Pose & pose, bool & modified_pose, utility::tag::TagPtr protocol_tag);

	//@brief Temporary hook into parsing machinery with pose reference
	MoverOP parse_protocol_tag(Pose & pose, utility::tag::TagPtr protocol_tag);

	//@brief Temporary hook into parsing machinery w/o pose reference.
	MoverOP parse_protocol_tag(utility::tag::TagPtr protocol_tag);

	void register_factory_prototypes();

private:
	//MoverFactoryOP mover_factory_;

	//DockDesignFilterFactoryOP ddfilter_factory_;

	static
	void
	substitute_variables_in_stream( std::istream & in, utility::options::StringVectorOption const& script_vars, std::stringstream & out);

}; // Parser

} // namespace jd2
} // namespace protocols

#endif //INCLUDED_protocols_jd2_DockDesignParser_HH
