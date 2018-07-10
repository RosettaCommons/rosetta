// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/simple_metrics/PerResidueStringMetric.cc
///
/// @brief Main class for simple metrics.
/// @author Jared Adolf-Bryfogle ( jadolfbr@gmail.com )

// Unit Headers
#include <core/simple_metrics/PerResidueStringMetric.hh>
#include <core/simple_metrics/util.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/residue_selector/ResidueSelector.hh>

// Protocol Headers

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/extra_pose_info_util.hh>

#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/util.hh>

namespace core {
namespace simple_metrics {

PerResidueStringMetric::PerResidueStringMetric():
	SimpleMetric("PerResidueStringMetric")
{}


PerResidueStringMetric::~PerResidueStringMetric() = default;

PerResidueStringMetric::PerResidueStringMetric( PerResidueStringMetric const & src ):
	SimpleMetric( src )
{


}

utility::vector1< std::string >
PerResidueStringMetric::get_metric_names() const {
	utility::vector1< std::string > names;
	names.push_back( metric() );
	return names;
}

void
PerResidueStringMetric::set_residue_selector(select::residue_selector::ResidueSelectorCOP selector){
	selector_ = selector;
}

void
PerResidueStringMetric::set_output_as_pdb_nums(bool output_as_pdb_nums){
	output_as_pdb_nums_ = output_as_pdb_nums;
}

void
PerResidueStringMetric::parse_per_residue_tag(utility::tag::TagCOP tag, basic::datacache::DataMap & datamap){
	set_output_as_pdb_nums(tag->getOption< bool >("output_as_pdb_nums", false));
	if ( tag->hasOption("residue_selector") ) {
		set_residue_selector(select::residue_selector::parse_residue_selector( tag, datamap ));
	}
}

///@brief Add options to the schema from this base class.
void
PerResidueStringMetric::add_schema( utility::tag::XMLSchemaComplexTypeGeneratorOP ct_gen) {
	add_per_residue_simple_metric_schema( ct_gen );
}

void
PerResidueStringMetric::apply( pose::Pose & pose, std::string prefix, std::string suffix ) const {
	std::map< core::Size, std::string > const values = calculate( pose );

	std::string custom_type = get_custom_type();
	if ( custom_type != "" ) custom_type=custom_type+"_";

	for ( auto value_pair : values ) {

		std::string out_number = utility::to_string(value_pair.first);

		if ( output_as_pdb_nums_ ) {
			out_number = pose.pdb_info()->pose2pdb(value_pair.first);
		}

		std::string const out_tag = prefix + out_number + "_" + custom_type + metric() + suffix;
		//std::cout << out_tag << std::endl;
		core::pose::setPoseExtraScore( pose, out_tag, value_pair.second);
	}
}

} //namespace simple_metrics
} //namespace core
