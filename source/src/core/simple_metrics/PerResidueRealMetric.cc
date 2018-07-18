// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/simple_metrics/PerResidueRealMetric.cc
///
/// @brief Main class for simple metrics.
/// @author Jared Adolf-Bryfogle ( jadolfbr@gmail.com )

// Unit Headers
#include <core/simple_metrics/PerResidueRealMetric.hh>
#include <core/simple_metrics/util.hh>

// Protocol Headers

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>

#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/util.hh>


namespace core {
namespace simple_metrics {
using namespace core::select::residue_selector;

PerResidueRealMetric::PerResidueRealMetric():
	SimpleMetric("PerResidueRealMetric")
{
	selector_ = TrueResidueSelectorOP( new TrueResidueSelector());
}


PerResidueRealMetric::~PerResidueRealMetric() = default;

PerResidueRealMetric::PerResidueRealMetric( PerResidueRealMetric const & src ):
	SimpleMetric( src )
{


}

utility::vector1< std::string >
PerResidueRealMetric::get_metric_names() const {
	utility::vector1< std::string > names;
	names.push_back( metric() );
	return names;
}

void
PerResidueRealMetric::set_residue_selector(select::residue_selector::ResidueSelectorCOP selector){
	selector_ = selector;
}

void
PerResidueRealMetric::set_output_as_pdb_nums(bool output_as_pdb_nums){
	output_as_pdb_nums_ = output_as_pdb_nums;
}

select::residue_selector::ResidueSelectorCOP
PerResidueRealMetric::get_selector() const {
	return selector_;
}

void
PerResidueRealMetric::parse_per_residue_tag(utility::tag::TagCOP tag, basic::datacache::DataMap & datamap){
	set_output_as_pdb_nums(tag->getOption< bool >("output_as_pdb_nums", false));
	if ( tag->hasOption("residue_selector") ) {
		set_residue_selector(select::residue_selector::parse_residue_selector( tag, datamap ));
	}
}

///@brief Add options to the schema from this base class.
void
PerResidueRealMetric::add_schema( utility::tag::XMLSchemaComplexTypeGeneratorOP ct_gen) {
	add_per_residue_simple_metric_schema( ct_gen );
}

void
PerResidueRealMetric::apply( pose::Pose & pose, std::string prefix, std::string suffix ) const {
	std::map< core::Size, core::Real > const values = calculate( pose ); //Index to value map

	std::string custom_type = get_custom_type();
	if ( custom_type != "" ) custom_type=custom_type+"_";

	//Output sum as well.
	core::Real sum = 0;
	std::string strip=" ";
	for ( auto value_pair : values ) {

		std::string out_number = utility::to_string(value_pair.first);

		if ( output_as_pdb_nums_ ) {

			out_number = utility::remove_from_string(pose.pdb_info()->pose2pdb(value_pair.first), strip);
		}

		sum += value_pair.second;
		std::string const out_tag = prefix  + custom_type + metric() + "_"+ out_number + suffix;
		core::pose::setPoseExtraScore( pose, out_tag, value_pair.second);
	}

	std::string const out_tag = prefix + "SUM" + "_" + custom_type + metric() + suffix;
	core::pose::setPoseExtraScore(pose, out_tag, sum);
}





} //namespace simple_metrics
} //namespace core
