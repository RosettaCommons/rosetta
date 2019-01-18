// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/jd3/util.hh
/// @brief Utility functions for JD3.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_jd3_util_hh
#define INCLUDED_protocols_jd3_util_hh

#include <protocols/jd3/InnerLarvalJob.fwd.hh>

// Core headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/keys/OptionKeyList.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace jd3 {

///@brief Append Residue Selector elements for your JD file.
void
append_residue_selector_subelements(
	utility::tag::XMLSchemaDefinition & job_definition_xsd,
	utility::tag::XMLSchemaComplexTypeGenerator & job_ct_gen );

///@brief Append the ability to have a SINGLE ScoreFunction set.
///
///@details Used to have the scorefunction defined for the entire job.
void
append_single_scorefxn_subelement(
	utility::tag::XMLSchemaDefinition & job_definition_xsd,
	utility::tag::XMLSchemaComplexTypeGenerator & job_ct_gen );

///@brief Add all residue selectors to a given datamap from a tag.
void
add_residue_selectors_to_datamap(
	core::pose::Pose const & pose,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap);

///@brief Return a single residue selector from tag.
core::select::residue_selector::ResidueSelectorCOP
parse_single_residue_selector(
	core::pose::Pose const & pose,
	utility::tag::TagCOP tag,
	std::string const & selector_name );


///////////////////////////////////////////////////////////////////////

///@brief Get the nstruct from the job_tag, or fallback to the cmd-line
core::Size
nstruct_for_job( utility::tag::TagCOP job_tag );

///@brief Create a local options collection from an InnerLarvalJob.
utility::options::OptionCollectionOP
options_for_job(
	utility::options::OptionKeyList const & option_keys,
	InnerLarvalJob const & inner_job,
	utility::tag::TagCOP common_block_tags );

///@brief Create a local options collection from a tag.
utility::options::OptionCollectionOP
options_from_tag(
	utility::options::OptionKeyList const & option_keys,
	utility::tag::TagCOP job_options_tags,
	utility::tag::TagCOP common_block_tags );


///////////////////////////////////////////////////////////////////////
///@brief Prints the job template to a Tracer.
void
print_job_template();

///@brief Get the OptionType from a key
utility::options::OptionTypes
option_type_from_key(
	utility::options::OptionKey const & key );

///@brief Get the XMLSchemaType for an option key.
/// Used during SJQ's xsd parsing
utility::tag::XMLSchemaType
value_attribute_type_for_option(
	utility::options::OptionTypes const & key );



} //protocols
} //jd3


#endif //protocols/jd3_util_hh

