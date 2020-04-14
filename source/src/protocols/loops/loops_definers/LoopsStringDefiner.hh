// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/loops_definers/LoopsStringDefiner.hh
/// @brief  Creates a serialized loops list based on a string specification
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_loops_loops_definers_LoopsStringDefiner_HH
#define INCLUDED_protocols_loops_loops_definers_LoopsStringDefiner_HH

// Unit Headers
#include <protocols/loops/loops_definers/LoopsDefiner.hh>
#include <protocols/loops/loops_definers/LoopsStringDefiner.fwd.hh>

#ifdef WIN32
#include <protocols/loops/Loop.hh>
#else
#include <protocols/loops/Loop.fwd.hh>
#endif

// Platform Headers
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <string>


namespace protocols {
namespace loops {
namespace loops_definers {

/// @brief A LoopsDefiner which can be configured with a string in a
/// Start:End:Cut,Start:End:Cut... format. (Cut optional, PDB numbering acceptable.)
class LoopsStringDefiner : public LoopsDefiner {
public:

	LoopsStringDefiner(std::string const & in = "");

	~LoopsStringDefiner() override;

	LoopsStringDefiner(
		LoopsStringDefiner const & src);

	/// @brief Create another loops definer of the type matching the most-derived
	/// version of the class.
	LoopsDefinerOP
	clone() const override;

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap const & data
	) override;


	SerializedLoopList
	apply(
		core::pose::Pose const &) override;

	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	std::string loop_spec_;

};

} //namespace
} //namespace
} //namespace

#endif
